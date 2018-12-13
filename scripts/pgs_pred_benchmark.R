#git version
library(data.table)
library(glmnet)
library(ade4)
library(ROCR)

source("analyze_functions.R")

#leave these fixed
do.binary <- F
nperm  <- 100

if (!exists("run_external")){
    nparc <- 50
    chip  <- "MEGA"
    pgs   <- "education"
    cv.method <- "famaware"
    #cv.method <- "standard"
}

#strings for output
bstring <- "regress"
if (do.binary)
  bstring <- "classify"

pstring <- "noperm"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")

ofname <- paste("../results/", paste("benchmark", chip, nparc, pgs, bstring, pstring, sep="_"), ".RData", sep="")

rdata.fname <- list()
rdata.fname[[1]] <- paste("../data/netmats1_", nparc, ".RData", sep="")
rdata.fname[[2]] <- paste("../data/netmats2_", nparc, ".RData", sep="")
pgs.fname   <- paste("../data/",pgs, "/", pgs, "_", chip, ".all.score",sep="")
fam.fname <- paste("../data/MEGA_Chip.fam")

message("loading ", pgs, " PGS")
pgs.info <- read.table(pgs.fname, head=T)
rownames(pgs.info) <- pgs.info$IID

ddim <- 0
if ( nparc != ddim){
    message("(re-)loading .RData")
    load(rdata.fname[[1]])
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

excl <- is.na(Y)
mydata.scale2 <- mydata.scale[!excl,]
Y2 <- Y[!excl]
myfam.use2 <- myfam.use[!excl,]

#mean 0 and sd 1
Y2 <- scale(Y2)
#netmat1 data
mydata.scale2.nm1 <- mydata.scale2

#netmat2 data
load(rdata.fname[[2]])
mydata.scale2.nm2 <- mydata.scale[!excl,]

bench <- c()
for(rep in 1:nperm){
  message(rep, " of ", nperm)

  ##fix the cv for runs with nm1 and nm2
  ##this allows for a paired t-test/wilcoxon test
  if (cv.method == "standard"){
    mycv      <- getCVfold(mydata.scale2, 10)
  }
  if (cv.method == "famaware"){
    mycv      <- getFamCVfold(mydata.scale2, 10, myfam.use2)
  }

  #nested CV
  nested.cv.nm1 <- doubleCV(Y2, mydata.scale2.nm1, myfold=mycv, fam="gaussian", measure="mse", lop="min")
  nested.cv.nm2 <- doubleCV(Y2, mydata.scale2.nm2, myfold=mycv, fam="gaussian", measure="mse", lop="min")

  xperf <- c(cmp.CV.perf(nested.cv.nm1), cmp.CV.perf(nested.cv.nm2))
  bench <- rbind(bench, xperf)
  print(bench)
  save(bench, file=ofname)
}

### END ###
