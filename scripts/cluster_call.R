
#default parameters
nparc <- 50
chip  <- "MEGA"
#pgs   <- "cad"
pgs   <- ""
do.binary <- F
nperm  <- 0
run_external <- T


##reading command line parameters
args = commandArgs(trailingOnly=F)
print(args)
for(aa in args){
    tok = strsplit(aa, split="=")[[1]]
    if (length(tok) == 2){
        if (tok[1] == "--nparc")
          nparc <- as.integer(tok[2])
        if (tok[1] == "--chip")
            chip <- tok[2]
        if (tok[1] == "--pgs")
            pgs <- tok[2]
        if (tok[1] == "--nperm")
            nperm <- as.integer(tok[2])
    }
}

message("options:")
for(var in c("nparc","chip","pgs","nperm"))
  message(var, ": ", get(var))

source("pgs_pred.R")
