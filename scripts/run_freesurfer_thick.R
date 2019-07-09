library(data.table)

run_external <- T

nparc <- 50
keep.sex <- "both"
regress.icv <- T
do.binary <- F
nperm <- 0
ceu_cut <- 0.0

icvstring <- "noICVreg"
if (regress.icv){
  icvstring <- "ICVreg"
}

fs.info <- fread("../data/unrestricted_hcp_freesurfer.csv", data.table=F)

fs.thick <- colnames(fs.info)[grep("Thck$",colnames(fs.info))]

ofnameX <- paste(paste("../results/FS_pred", nparc, keep.sex, ceu_cut, icvstring,sep="_" ), ".RData", sep="")

total.res <- data.frame()
for(fs_feat in fs.thick){
  message(fs_feat)
  source("./fs_feat_pred.R")
  tmp <- data.frame(name=fs_feat, cor=xperf[1], rmse=xperf[2])
  total.res <- rbind(total.res, tmp)
  print(total.res)
  save(total.res, file=ofnameX)
}
