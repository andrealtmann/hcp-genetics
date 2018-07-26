library(data.table)

#currently works only for 50 and 300
nparc <- 300

data.path <- paste("/Volumes/My Passport/HCP/HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d", nparc, "_ts2/", sep="")
fname <- paste(data.path, "netmats1.txt", sep="")

message("loading data")
subjects <- read.table(paste(data.path,"/../../subjectIDs.txt",sep=""))
subject.info <- read.csv(paste(data.path,"/../../../unrestricted_altmann2_7_2_2018_17_45_27.csv",sep=""))
rownames(subject.info) <- subject.info[,"Subject"]

ethn.info <- read.table(paste(data.path, "/../../../snpweights/NA.predpc", sep=""))
rownames(ethn.info) <- ethn.info[,1]
colnames(ethn.info)[7:10] <- c("YRI","CEU","ASI","NAT")

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
  myidx <- c(myidx, sstart:eend)
}

mydata <- tmp[,myidx]
rm(tmp)

##standardize features
message("scaling features")
mydata.scale <- apply(mydata, 2, scale)
rownames(mydata.scale) <- rownames(mydata)

ofname <- paste("../data/netmats1_", nparc, ".RData",sep="")
message("storing in: ", ofname)
save(subject.info, ethn.info, nparc, mydata.scale, file=ofname)

### END ###
