library(ggplot2)
library(gridExtra)

chip <- "MEGA"
nperm <- 100
#pgs <- "education2"
#pgs <- "ad"
mypgs <- c("ad","blood_pressure","bmi","cad","education2","eps","fiq","height","mdd","sz")
nparcs <- c(15,50, 100, 300)
cv.methods <- c("famaware","standard")

myplots <- list()

for(pgs in mypgs){
#loop over pgs
grp <- 0
result <- data.frame()
for (nparc in nparcs){
  for(cvme in cv.methods){
    fname <- paste("../results/benchmark_", cvme, "_", chip, "_", nparc, "_", pgs, "_regress_perm",nperm,".RData",sep="")
    load(fname)

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

#if only one of the nm's
result <- subset(result, ct=="nm2")
median.perf <- sapply(nparcs, function(np){
  #tmp <- subset(result, nd==np & cvm=="standard")
  tmp <- subset(result, nd==np & cvm=="famaware")
  median(tmp$cor)  
})
message(pgs,": ", paste(round(median.perf,2), collapse="  "))


#xlab <- paste(rep(c("nm1","nm2"), 2*length(nparcs)), rep(nparcs, each=4))

#p <- ggplot(result, aes(factor(grp), cor ))
#p <- p +  geom_violin(aes(fill=factor(cvm))) + geom_boxplot(width=0.1)
#geom_jitter(height=0, width=0.01, size=0.5)
#p <- p + scale_x_discrete(name="", labels=xlab) + theme(axis.text.x= element_text(angle=45)) + ylab("correlation")

#p +  scale_color_grey()+ theme_classic()

p <- ggplot(result, aes(x=factor(nd), y=cor, fill=factor(cvm), group=factor(grp))) + geom_violin() + theme_bw()
p <- p + ylab("correlation") + xlab("ICA") + ggtitle(pgs)
myplots[[pgs]] <- p + ylim(-0.2,0.4) + theme(legend.position="none", axis.title.x = element_blank())
}

grid.arrange(grobs=myplots)