

chip<-"MEGA"
if (!exists("nparc"))
  nparc <- 300 
bstring <- "regress"
nperm <- 10000
do.binary <- F

bstring <- "regress"
if (do.binary)
  bstring <- "classify"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")

all.pgs <- c("ad","bmi","blood_pressure","cad","education","education2","eps","fiq","height","mdd","sz")
#just to test the script
#all.pgs <- c("education")

res <- c()
for (mypgs in all.pgs){
  fname <- paste("../results/", paste(chip, nparc, mypgs, bstring, pstring, sep="_"), ".RData", sep="")
  load(fname)


  r.one.sided.p <- cor.test(unlist(nested.cv$labels), unlist(nested.cv$prediction), alternative="greater")$p.value
  r.perm.p  <- (sum(rndm[,"cor"] >= xperf["cor"])  + 1)/(nrow(rndm) + 1)
  ms.perm.p <- (sum(rndm[,"rmse"] <= xperf["rmse"]) + 1)/(nrow(rndm) + 1)
  tmp <- c(xperf["cor"], r.one.sided.p, r.perm.p, xperf["rmse"], ms.perm.p, nrow(rndm))
  res <- rbind(res, tmp)
}
#res <- data.frame(res)
rownames(res) <- all.pgs
colnames(res) <- c("r", "r_cor_test_p", "r_perm_p", "rmse", "rmse_perm_p", "nperm")
