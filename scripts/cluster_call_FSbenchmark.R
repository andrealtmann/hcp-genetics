
#default parameters
nparc <- 50
#fs_feat   <- ""
fs_feat_number <- 1
nperm  <- 100
keep.sex <- "both"
regress.icv <- T
run_external <- T
ceu_cut <- 0.0
cv.method <- ""

##reading command line parameters
args = commandArgs(trailingOnly=F)
print(args)
for(aa in args){
    tok = strsplit(aa, split="=")[[1]]
    if (length(tok) == 2){
        if (tok[1] == "--nparc")
          nparc <- as.integer(tok[2])
        #if (tok[1] == "--fsf")
        #    fs_feat <- tok[2]
        if (tok[1] == "--fsfn")
            fs_feat_number <- as.integer(tok[2])
        if (tok[1] == "--sex")
            keep.sex <- tok[2]
        if (tok[1] == "--ICVreg")
            regress.icv == as.logical(tok[2])
        if (tok[1]== "--CEU")
            ceu_cut <- as.double(tok[2])
        if (tok[1] == "--nperm")
            nperm <- as.integer(tok[2])
	    if (tok[1] == "--cv")
	       cv.method <- tok[2]
    }
}

#convert fs_feat_number into fs_feat
fs_map <- read.csv("../data/fs_feature_list.csv", row.names=1)
fs_feat <- paste(fs_map[fs_feat_number,1])

message("options:")
for(var in c("nparc","cv.method","keep.sex","fs_feat","regress.icv","ceu_cut","nperm"))
  message(var, ": ", get(var))

source("fs_feat_pred_benchmark.R")
