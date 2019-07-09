library(data.table)

run_external <- T

nparc <- 50
keep.sex <- "male"
regress.icv <- T
do.binary <- F
nperm <- 0
ceu_cut <- 0.0

fs.info <- fread("../data/unrestricted_hcp_freesurfer.csv", data.table=F)

fs.thick <- colnames(fs.info)[grep("Thck",colnames(fs.info))]

for(fs_feat in fs.thick){
  message(fs_feat)
  source("./fs_feat_pred.R")
  print(c(fs_feat, xperf))
}



