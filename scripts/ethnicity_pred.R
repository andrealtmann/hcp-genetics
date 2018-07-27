library(data.table)
library(glmnet)
library(ade4)
library(ROCR)

source("analyze_functions.R")

if (!exists("run_external")){
    nparc <- 300
    chip  <- "MEGA"
    nperm  <- 0
}

pstring <- "noperm"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")


ofname <- paste("../results/", paste(chip, nparc, "ethn",  pstring, sep="_"), ".RData", sep="")

rdata.fname <- paste("../data/netmats1_", nparc, ".RData", sep="")

#avoid re-loading data for multiple runs
if ( !exists("mydata.scale"))
  ddim <- 0
if ( nparc != ddim){
    message("(re-)loading .RData")

    load(rdata.fname)
    ddim <- nparc
}

##add sex info
sex <- subject.info[rownames(ethn.info),"Gender"]
ethn.info <- data.frame(ethn.info, sex)

##rm subjects not of > 0.90 CEU ethn.
#keep.ceu  <- rownames(subset(ethn.info, CEU>0.9))
#keep.ceu  <- rownames(subset(ethn.info, CEU>0.9 & sex=="M"))
#keep.ceu  <- rownames(subset(ethn.info, CEU>0.9 & sex=="F"))

#outcome
Y <- ethn.info[rownames(mydata.scale),"CEU"] > 0.9

excl <- is.na(Y)
mydata.scale2 <- mydata.scale[!excl,]
Y2 <- Y[!excl]


##dummy
message("10 fold CV")
abc <- cv.glmnet(mydata.scale2, Y2, type.measure = "auc", nfolds=10, family="binomial", alpha=0.5)
#plot(abc)

#nested CV
mycv      <- getCVfold(mydata.scale2, 10)
system.time(nested.cv <- doubleCV(Y2, mydata.scale2, myfold=mycv))

mauc <- mean(unlist(performance(prediction(nested.cv$prediction, nested.cv$labels), 'auc')@y.values))
oauc <- unlist(performance(prediction(unlist(nested.cv$prediction), unlist(nested.cv$labels)), 'auc')@y.values)

##permuation
rndm <- c()
if (nperm > 0){
  for(j in 1:nperm){
    message("rnd", j)
    Yrn <- sample(Y2)
    tmp.cv <- doubleCV(Yrn, mydata.scale2, myfold=mycv, nla=20, fam="binomial", measure="auc", lop="1se")

    rmauc <- mean(unlist(performance(prediction(tmp.cv$prediction, tmp.cv$labels), 'auc')@y.values))
    roauc <- unlist(performance(prediction(unlist(tmp.cv$prediction), unlist(tmp.cv$labels)), 'auc')@y.values)


    rndm <- rbind(rndm, c(rmauc, roauc))
    print(rndm)
  }
}
message("save results to ", ofname)
save(abc, mycv, nested.cv, mauc, oauc, rndm, file=ofname)



### END ###
