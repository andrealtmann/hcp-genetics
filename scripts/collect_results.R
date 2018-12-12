

chip<-"MEGA"
nparc <- 300
bstring <- "regress"
nperm <- 250
do.binary <- F

bstring <- "regress"
if (do.binary)
  bstring <- "classify"
if (nperm > 0)
  pstring <- paste("perm",nperm, sep="")

all.pgs <- c("ad","blood_pressure","cad","education","eps","fiq","height","mdd","sz")

res <- c()
for (mypgs in all.pgs){
  fname <- paste("../results/", paste(chip, nparc, mypgs, bstring, pstring, sep="_"), ".RData", sep="")
  load(fname)

  
  r.one.sided.p <- cor.test(unlist(nested.cv$labels), unlist(nested.cv$prediction), alternative="greater")$p.value
  r.perm.p  <- (sum(rndm[,1] > rtr)  + 1)/(nrow(rndm) + 1)
  ms.perm.p <- (sum(rndm[,2] < mstr) + 1)/(nrow(rndm) + 1)
  tmp <- c(rtr, r.one.sided.p, r.perm.p, mstr, ms.perm.p, nrow(rndm))
  res <- rbind(res, tmp)
}
#res <- data.frame(res)
rownames(res) <- all.pgs
colnames(res) <- c("r", "r_cor_test_p", "r_perm_p", "rmse", "rmse_perm_p", "nperm")




