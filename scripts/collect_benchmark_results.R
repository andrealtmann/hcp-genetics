library(ggplot2)
library(gridExtra)

chip <- "MEGA"
nperm <- 100
#pgs <- "education2"
#pgs <- "ad"
mypgs <- c("ad","blood_pressure","bmi","cad","education2","eps","fiq","height","mdd","sz")
labmap <- c("AD", "BP", "BMI", "CAD", "EDU", "EPS", "FIQ", "HEI", "MDD", "SZ")
names(labmap) <- mypgs

nparcs <- c(15, 50, 100, 300)
cv.methods <- c("famaware","standard")

metr <- "cor"

perm.test <- T

myplots <- list()
all.median.perf <- c()
xwilcp <- c()

for(pgs in mypgs){
#loop over pgs
grp <- 0
result <- data.frame()
for (nparc in nparcs){
  for(cvme in cv.methods){
    fname <- paste("../results/benchmark_", cvme, "_", chip, "_", nparc, "_", pgs, "_regress_perm",nperm,".RData",sep="")
    load(fname)
    if (nrow(bench) != nperm){
      message("WARNING: ", cvme, " ", nparc, " ", pgs, ": ", nrow(bench))
    }
    grp <- grp+1
    #grp <- paste("nm1", cvme, sep="_")
    tmp.data <-  data.frame(cbind(bench[,1:2], ct="nm1", nd=nparc, cvm=cvme, grp=grp))
    grp <- grp+1
    #grp <- paste("nm2", cvme, sep="_")
    tmp2.data <- data.frame(cbind(bench[,3:4], ct="nm2", nd=nparc, cvm=cvme, grp=grp))
    tmp.data <- rbind(tmp.data, tmp2.data)
    result <- rbind(result, tmp.data)
  }
}
result$cor <- as.double(paste(result$cor))
result$rmse <- as.double(paste(result$rmse))

##

#if only one of the nm's
result <- subset(result, ct=="nm2")
median.perf <- sapply(nparcs, function(np){
  #tmp <- subset(result, nd==np & cvm=="standard")
  tmp <- subset(result, nd==np & cvm=="famaware")
  median(tmp[, metr])  
})
#message(pgs,": ", paste(round(median.perf,2), collapse="  "))
all.median.perf <- rbind(all.median.perf, median.perf)

  for(nd in nparcs){
    A <- subset(result, subset=(nd == nd) & (cvm=="standard"))
    B <- subset(result, subset=(nd == nd) & (cvm=="famaware"))
    xwilcp <- c(xwilcp, wilcox.test(A[,metr], B[,metr], alternative="greater")$p.value)
  }


if (metr == "cor"){
  p <- ggplot(result, aes(x=factor(nd), y=cor, fill=factor(cvm), group=factor(grp))) + geom_violin() + theme_bw()
  p <- p + ylab("correlation") + xlab("ICA") + ggtitle(labmap[pgs])
myplots[[pgs]] <- p + ylim(-0.2,0.4) + theme(legend.position="none", axis.title.x = element_blank())
}
if (metr == "rmse"){
  p <- ggplot(result, aes(x=factor(nd), y=rmse, fill=factor(cvm), group=factor(grp))) + geom_violin() + theme_bw()
  p <- p + ylab("RMSE") + xlab("ICA") + ggtitle(labmap[pgs])
myplots[[pgs]] <- p + ylim(0.95,1.05) + theme(legend.position="none", axis.title.x = element_blank())

}
}

rownames(all.median.perf) <- mypgs
colnames(all.median.perf) <- nparcs

grid.arrange(grobs=myplots)

all.pvalues.perf <- all.median.perf
if (perm.test){
  npermtest <- 1000
  for(pgs in mypgs){
    for(nd in nparcs){
      rndname <- paste("../results/",chip, "_", nd, "_", pgs, "_regress_FAMperm",npermtest, ".RData", sep="")
      load(rndname)
      if (nrow(rndm) != npermtest){
        message("WARNING: ", nd, " ", pgs, ": ", nrow(rndm))
	if (nrow(rndm) > npermtest){
	  #message("WARNING: using only first ", npermtest, " results")
	  rndm <- rndm[1:npermtest,]
	}
      }
      #true perm p
      ndem <- nrow(rndm)
      if (metr=="cor")
        ppval <- (sum(rndm[,metr] >= all.median.perf[pgs,paste(nd)]) + 1) / (ndem + 1)
      if (metr=="rmse")
        ppval <- (sum(rndm[,metr] <= all.median.perf[pgs,paste(nd)]) + 1) / (ndem + 1)      
      #parametric perm p
      #ppval <- pnorm(all.median.perf[pgs, paste(nd)], mean=mean(rndm[,metr]), sd=sd(rndm[,metr]), lower.tail=F)
      all.pvalues.perf[pgs, paste(nd)] <- ppval
    }
  }
}

##ugly hack, but works if now copy paste ;)
imgplot <- function(){
  image(-log10(all.pvalues.perf), zlim=c(0,3), axes=F)
  #indicate significance
  nom   <- all.pvalues.perf <= 0.05
  subs  <- all.pvalues.perf <= 0.01
  vsubs <- all.pvalues.perf <= 0.005
  collate <- nom + subs + vsubs
  
  xs <- seq(0,1,1/(length(mypgs)-1))
  ys <- seq(0,1,1/(length(nparcs)-1))

  yidx <- 0
  for(nd in nparcs){
    yidx <- yidx + 1
    xidx <- collate[, paste(nd)] == 1
    if (sum(xidx) > 0)
      text(xs[xidx], ys[yidx], labels=rep("*", sum(xidx)))
    xidx <- collate[, paste(nd)] == 2
    if (sum(xidx) > 0)
      text(xs[xidx], ys[yidx], labels=rep("**", sum(xidx)))
    xidx <- collate[, paste(nd)] == 3
    if (sum(xidx) > 0)
      text(xs[xidx], ys[yidx], labels=rep("***", sum(xidx)))
  }

  box()
  axis(1,at=xs, labels=labmap[mypgs], las=2)
  axis(2,at=ys, labels=nparcs, las=2)
}

imgplot2 <- function(){
  xdf <- data.frame()
  for(i in 1:ncol(all.pvalues.perf)){
    tmp <- data.frame(pv=all.pvalues.perf[,i], trait=labmap[rownames(all.pvalues.perf)], ica=nparcs[i])
    xdf <- rbind(xdf, tmp)
  }

  sigid <- rep("",nrow(xdf))
  sigid[xdf$pv <= 0.05] <- "*"
  sigid[xdf$pv <= 0.01] <- "**"
  sigid[xdf$pv <= 0.005] <- "***"
  xdf <- data.frame(xdf, xtxt=sigid)
  print(xdf)

  xdf$pv <- -log10(xdf$pv) 
  
  p <- ggplot(xdf, aes(trait, as.factor(ica))) + geom_raster(aes(fill=pv)) 
  p <- p + scale_fill_gradientn(colours=heat.colors(12)) + labs(fill="-log10(p)")
  p <- p + theme_minimal() +  theme(axis.title.y=element_blank(), axis.title.x = element_blank())
  p <- p + geom_text(aes(label=xtxt))
  return(p)

}