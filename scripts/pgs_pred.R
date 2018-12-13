#git version
library(data.table)
library(glmnet)
library(ade4)
library(ROCR)

source("analyze_functions.R")

if (!exists("run_external")){
    nparc <- 50
    chip  <- "MEGA"
    pgs   <- "education"
    do.binary <- F
    nperm  <- 0
}

#strings for output
bstring <- "regress"
if (do.binary)
  bstring <- "classify"

pstring <- "noperm"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")

ofname <- paste("../results/", paste(chip, nparc, pgs, bstring, pstring, sep="_"), ".RData", sep="")

rdata.fname <- paste("../data/netmats2_", nparc, ".RData", sep="")
pgs.fname   <- paste("../data/",pgs, "/", pgs, "_", chip, ".all.score",sep="")
fam.fname <- paste("../data/MEGA_Chip.fam")

message("loading ", pgs, " PGS")
pgs.info <- read.table(pgs.fname, head=T)
rownames(pgs.info) <- pgs.info$IID

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
keep.ceu  <- rownames(subset(ethn.info, CEU>0.9))
#keep.ceu  <- rownames(subset(ethn.info, CEU>0.9 & sex=="M"))
#keep.ceu  <- rownames(subset(ethn.info, CEU>0.9 & sex=="F"))
pgs2      <- pgs.info[keep.ceu, "X0.700000"]
names(pgs2) <- keep.ceu

#outcome
Y <- pgs2[rownames(mydata.scale)]

fam.info <- read.table(fam.fname)
myfam.map <- fam.info[,1:2]
rownames(myfam.map) <- myfam.map[,2]
myfam.use <- myfam.map[rownames(mydata.scale),]

if (do.binary){
    #binarize the outcome and do classification top 25% vs bottom 25%
    ccc <- quantile(Y,c(0.25,0.75),na.rm=T)
    #ccc <- quantile(Y,c(0.49,0.51),na.rm=T)
    Y3 <- rep(NA, length(Y))
    Y3[Y <= ccc[1]] <- F
    Y3[Y >= ccc[2]] <- T

    excl <- is.na(Y3)
    mydata.scale2 <- mydata.scale[!excl,]
    Y2 <- Y3[!excl]
    myfam.use2 <- myfam.use[!excl,]

    #rndm
    #Y2 <- sample(Y2, length(Y2), replace=F)

    ##dummy
    message("10 fold CV")
    abc <- cv.glmnet(mydata.scale2, Y2, type.measure = "auc", nfolds=10, family="binomial", alpha=0.5)
    plot(abc)

    #nested CV
    #mycv      <- getCVfold(mydata.scale2, 10)
    #nested CV with outer fold being family aware!
    mycv      <- getFamCVfold(mydata.scale2, 10, myfam.use2)
    system.time(nested.cv <- doubleCV(Y2, mydata.scale2, myfold=mycv))

    aucs <- unlist(performance(prediction(nested.cv$prediction, nested.cv$labels),'auc')@y.values)
} else {
    #do regression with the full dataset

    excl <- is.na(Y)
    mydata.scale2 <- mydata.scale[!excl,]
    Y2 <- Y[!excl]
    myfam.use2 <- myfam.use[!excl,]

    #mean 0 and sd 1
    Y2 <- scale(Y2)

    message("10 fold CV")
    abc <- cv.glmnet(mydata.scale2, Y2, type.measure = "mse", nfolds=10, family="gaussian", alpha=0.5)
    #plot(abc)

    #nested CV
    #mycv      <- getCVfold(mydata.scale2, 10)
    mycv      <- getFamCVfold(mydata.scale2, 10, myfam.use2)
    nested.cv <- doubleCV(Y2, mydata.scale2, myfold=mycv, fam="gaussian", measure="mse", lop="min")

    #rtr  <- cor(unlist(nested.cv$labels),  unlist(nested.cv$prediction))
    #mstr <- sqrt(mean((unlist(nested.cv$labels) -  unlist(nested.cv$prediction))^2))

    xperf <- cmp.CV.perf(nested.cv)

    rndm <- c()
    if (nperm > 0){
        for(j in 1:nperm){
          message("rnd", j)
          Yrn <- sample(Y2)
          tmp.cv <- doubleCV(Yrn, mydata.scale2, myfold=mycv, nla=20, fam="gaussian", measure="mse", lop="min")
          #rrn <- cor(unlist(tmp.cv$labels),  unlist(tmp.cv$prediction))
          #rms <- sqrt(mean((unlist(tmp.cv$labels) -  unlist(tmp.cv$prediction))^2))
          #rndm <- rbind(rndm, c(rrn, rms))
          rndperf <- cmp.CV.perf(tmp.cv)
          rndm <- rbind(rndm, rndperf)
          print(rndm)
        }
    }
    message("save results to ", ofname)
    save(abc, mycv, nested.cv, xperf, rndm, file=ofname)
}


### END ###
