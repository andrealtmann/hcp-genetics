#git version
library(data.table)
library(glmnet)
library(ade4)
library(ROCR)

source("analyze_functions.R")

#leave these fixed
#nperm  <- 100

if (!exists("run_external")){
    nperm  <- 100
    nparc <- 50
    fs_feat <- "FS_L_Hippo_Vol"
    #fs_feat <- "FS_InterCranial_Vol"
    #options: both, male, female
    keep.sex <- "both"
    #keep.sex <- "male"
    regress.icv <- T
    ceu_cut <- 0.0
}

#strings for output
bstring <- "regress"

icvstring <- "noICVreg"
if (regress.icv){
  icvstring <- "ICVreg"
}

pstring <- "noperm"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")

ofname <- paste("../results/", paste("FS", nparc, fs_feat, bstring, keep.sex, icvstring, paste("CEU",ceu_cut,sep=""), pstring, sep="_"), ".RData", sep="")
message(ofname)


rdata.fname <- list()
rdata.fname[[1]] <- paste("../data/netmats1_", nparc, ".RData", sep="")
rdata.fname[[2]] <- paste("../data/netmats2_", nparc, ".RData", sep="")

#FreeSurfer Data
fs.fname    <- "../data/unrestricted_hcp_freesurfer.csv"

#important for cross-validation
fam.fname <- paste("../data/MEGA_Chip.fam")

message("loading freeSurfer data")
fs.info <- fread(fs.fname, data.table=F)
rownames(fs.info) <- fs.info$Subject

ddim <- 0
if ( nparc != ddim){
    message("(re-)loading .RData")
    load(rdata.fname[[1]])
    ddim <- nparc
}

##add sex info
sex <- subject.info[rownames(ethn.info),"Gender"]
ethn.info <- data.frame(ethn.info, sex)

if (keep.sex != "both"){
  if (keep.sex == "male"){
    keep <- rownames(subset(ethn.info, sex=="M"))
  } else {
    keep <- rownames(subset(ethn.info, sex=="F"))
  }
} else {
  keep <- rownames(ethn.info)
}

##rm subjects not of > ceu_cut
#suggested settings either 0.0 to include all ethn
#or 0.9 to include only CEU subjects
keep.ceu <- rownames(subset(ethn.info[keep,], CEU>=ceu_cut))

fs2      <- fs.info[keep.ceu, fs_feat]
if (regress.icv){
  myicv <- fs.info[keep.ceu, "FS_InterCranial_Vol"]
  message("hello world")
  dummy <- lm(fs2 ~ myicv)
  fs2.bak <- fs2
  dummy2 <- dummy$res
  fs2 <- rep(NA, length(fs2.bak))
  fs2[as.integer(names(dummy2))] <- dummy2
}
names(fs2) <- keep.ceu

#outcome
Y <- fs2[rownames(mydata.scale)]

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

#it is n repetitions rather than n permutations. actually, misnomer
bench <- c()
for(rep in 1:nperm){
  message(rep, " of ", nperm)

  ##fix the cv for runs with nm1 and nm2
  ##this allows for a paired t-test/wilcoxon test
  if (cv.method == "standard"){
    message("sampling standard folds")
    mycv      <- getCVfold(mydata.scale2, 10)

  }
  if (cv.method == "famaware"){
    message("sampling family-aware folds")
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
