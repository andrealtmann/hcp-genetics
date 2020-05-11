#classic stats to test for diffeences
#between FS features in different ancestries

library(data.table)

fs.data <- fread("../data/unrestricted_hcp_freesurfer.csv", data.table=F)

#subject info comes from:
load("../data/netmats1_15.RData")

#fam_info 
fam.info <- fread("../data/MEGA_Chip.fam", data.table=F)


mydat <- merge(merge(subject.info[,1:5], ethn.info, by.x="Subject", by.y="V1"), fs.data)

feat <- colnames(mydat)[16:949]
honk <- t(sapply(feat, function(ff){
  message(ff)
  bum <- mydat[,ff]
  form <- paste("bum", "~ Gender + Age + FS_InterCranial_Vol + eval(CEU >= 0.9)")
  tmp <- lm(as.formula(form), data=mydat)
  drop1(tmp, test="F")[,6]
}))

cthick <- rownames(honk)[grep("_Thck$",rownames(honk))]

honk[cthick,][order(honk[cthick,5]),]