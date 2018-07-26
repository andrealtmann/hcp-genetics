library(data.table)
library(glmnet)
library(ade4)
library(ROCR)

source("analyze_functions.R")

if (!exists("run_external")){
    nparc <- 300
    chip  <- "MEGA"
    pgs   <- "fiq"
    do.binary <- F
    perm.rnd <- F
}

data.path <- paste("/Volumes/My Passport/HCP/HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d", nparc, "_ts2/", sep="")
fname <- paste(data.path, "netmats1.txt", sep="")

message("loading data")
subjects <- read.table(paste(data.path,"/../../subjectIDs.txt",sep=""))
subject.info <- read.csv(paste(data.path,"/../../../unrestricted_altmann2_7_2_2018_17_45_27.csv",sep=""))
rownames(subject.info) <- subject.info[,"Subject"]

ethn.info <- read.table(paste(data.path, "/../../../snpweights/NA.predpc", sep=""))
rownames(ethn.info) <- ethn.info[,1]
colnames(ethn.info)[7:10] <- c("YRI","CEU","ASI","NAT")

pgs.info <- read.table(paste(data.path,"/../../../prsice/",pgs,"/",pgs,"_",chip,".all.score",sep=""), head=T)
rownames(pgs.info) <- pgs.info$IID

#avoid re-loading data for multiple runs
if ( !exists("mydata.scale"))
  ddim <- 0
if ( nparc != ddim){
    message("(re-)loading connectivity data")

    tmp <- fread(fname, data.table=F)
    rownames(tmp) <- subjects[,1]

    ddim <- sqrt(ncol(tmp))

    if (nparc != ddim)
      message("WRONG DIMENSION")


    #build upper triangle of matrix
    message("extract upper triangle of matrix")
    myidx <- c()
    for (i in 1:(nparc-1)){
      sstart <- (i-1) * nparc + (i+1)
      eend   <- i * nparc
      #message(sstart, ":", eend)
      myidx <- c(myidx, sstart:eend)
    }

    mydata <- tmp[,myidx]
    rm(tmp)

    ##standardize features
    message("scaling features")
    mydata.scale <- apply(mydata, 2, scale)
    rownames(mydata.scale) <- rownames(mydata)
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

if (do.binary){
    #binarize the outcome and do classification
    ccc <- quantile(Y,c(0.25,0.75),na.rm=T)
    #ccc <- quantile(Y,c(0.49,0.51),na.rm=T)
    Y3 <- rep(NA, length(Y))
    Y3[Y <= ccc[1]] <- F
    Y3[Y >= ccc[2]] <- T

    excl <- is.na(Y3)
    mydata.scale2 <- mydata.scale[!excl,]
    Y2 <- Y3[!excl]

    #rndm
    #Y2 <- sample(Y2, length(Y2), replace=F)

    ##dummy
    message("10 fold CV")
    abc <- cv.glmnet(mydata.scale2, Y2, type.measure = "auc", nfolds=10, family="binomial", alpha=0.5)
    plot(abc)

    #nested CV
    mycv      <- getCVfold(mydata.scale2, 10)
    system.time(nested.cv <- doubleCV(Y2, mydata.scale2, myfold=mycv))
} else {
    #do regression with the full dataset

    excl <- is.na(Y)
    mydata.scale2 <- mydata.scale[!excl,]
    Y2 <- Y[!excl]

    #mean 0 and sd 1
    Y2 <- scale(Y2)

    message("10 fold CV")
    abc <- cv.glmnet(mydata.scale2, Y2, type.measure = "mse", nfolds=10, family="gaussian", alpha=0.5)
    plot(abc)

    #nested CV
    mycv      <- getCVfold(mydata.scale2, 10)
    system.time(nested.cv <- doubleCV(Y2, mydata.scale2, myfold=mycv, fam="gaussian", measure="mse", lop="min"))

    rtr  <- cor(unlist(nested.cv$labels),  unlist(nested.cv$prediction))
    mstr <- sqrt(mean((unlist(nested.cv$labels) -  unlist(nested.cv$prediction))^2))

    #library(kernlab)
    #check out RVR
    #xxx <- rvm(x=mydata.scale2, y=Y2, kernel="rbfdot", cross=10)

    if (perm.rnd){
        rndm <- c()
        for(j in 1:100){
          message("rnd", j)
          Yrn <- sample(Y2)
          tmp.cv <- doubleCV(Yrn, mydata.scale2, myfold=mycv, nla=20, fam="gaussian", measure="mse", lop="min")
          #rrn <- mean(tmp.cv$metrics[,3])
          #rms <- mean(tmp.cv$metrics[,4])
          rrn <- cor(unlist(tmp.cv$labels),  unlist(tmp.cv$prediction))
          rms <- sqrt(mean((unlist(tmp.cv$labels) -  unlist(tmp.cv$prediction))^2))
          rndm <- rbind(rndm, c(rrn, rms))
          print(rndm)
        }
    }

}

if (0){
    #old and experimental stuff
    break.point

    #train test
    N <- nrow(mydata.scale2)
    trte <- sample(1:N, floor(N/10)*5)

    testset <- setdiff(1:N, trte)

    tr.cv <- cv.glmnet(mydata.scale2[trte,], Y2[trte], type.measure = "auc", nfolds=10, family="binomial", alpha=0.5)
    olam  <- tr.cv$lambda.1se
    #olam  <- tr.cv$lambda.min
    olam.idx <- which(tr.cv$lambda == olam)

    yhat <- predict(tr.cv$glmnet.fit, newx=mydata.scale2[testset,], s=olam, type="response")
    miss <- (yhat > 0.5) != Y2[testset]

    performance(prediction(yhat, Y2[testset]),'auc')
    wilcox.test(yhat[Y2[testset]],yhat[!Y2[testset]])



    ### randomForest ###
    library(randomForest)
    rf.mod <- randomForest(mydata.scale2[trte,], Y2[trte], ntree=1000, importance=T)
}

### END ###
