library(data.table)
library(lme4)
library(glmnet)

source("analyze_functions.R")

##settings
ceu_cut <- 0.9
fs.feats <- c("Area","Thck")
#fs.feats <- c("Area")

#classic statistical analysis of FS (Area/Thck) wrt. ethnicity

#load this df to get subject information
load("../data/netmats1_15.RData")
colnames(ethn.info)[1] <- "Subject"

ttt <- colnames(ethn.info)[7:10]

MLETH <- ttt[apply(ethn.info[,ttt], 1, function(x){
  which.max(x)
})]
ethn.info <- data.frame(ethn.info, MLETH=MLETH)


sage <- read.table("../data/subject_age.txt", head=T)

#load FS data
fs_data <- fread("../data/unrestricted_hcp_freesurfer.csv", data.table=F)

#load FAM
fam <- fread("../data/MEGA_Chip.fam", data.table=F)
colnames(fam) <- c("FAM","Subject","P1","P2","SEX","PH")

confound <- subject.info[,c("Subject","Acquisition","Gender","Age")]
confound2 <- merge(confound, fam, by="Subject")
confound3 <- merge(confound2, sage, by="Subject")
basic_info <- merge(confound3, ethn.info, by="Subject")


#my_conf <- c("Gender","AGE","FS_InterCranial_Vol")
my_conf <- c("Gender","AGE")

my_result <- c()
for (fsfeat in fs.feats){
  fidx <- colnames(fs_data)[grep(paste("_",fsfeat,"$",sep=""),colnames(fs_data))]
  img_feat_tmp <- fs_data[,c("Subject","FS_InterCranial_Vol",fidx)]

  my.df <- merge(basic_info, img_feat_tmp, by="Subject")

  res <- sapply(fidx, function(my_feat){

    my_form <- paste(c(my_feat, "~ eval(CEU>",ceu_cut, ")", " + ", paste(my_conf, collapse = " + ") ), collapse=" ")
    #standard lm
    #tmp <- lm(as.formula(my_form), data=my.df)
    #summary(tmp)$coeff[2,]

    #lme
    my_form2 <- paste(c(my_feat, "~", paste(my_conf, collapse = " + ") ), collapse=" ")
    tmp2 <- lmer( as.formula(paste(my_form2, "+ (1 | FAM)")),  data=my.df, REML=F)
    tmp1 <- lmer( as.formula(paste(my_form, "+ (1 | FAM)")),  data=my.df, REML=F)
    return(anova(tmp1, tmp2)[2,8])
  })

  my_result <- cbind(my_result, res)
}
colnames(my_result) <- fs.feats


if (0){
##this is to run some predictions

#glmnet
rownames(fam) <- fam[,"Subject"]

for (fsfeat in fs.feats){
  fidx <- colnames(fs_data)[grep(paste("_",fsfeat,"$",sep=""),colnames(fs_data))]
  img_feat_tmp <- fs_data[,c("Subject","FS_InterCranial_Vol",fidx)]

  my.df <- merge(basic_info, img_feat_tmp, by="Subject")
  my.df2 <- sapply(fidx, function(feat){
    #my_res <- lm(as.formula(paste(feat, "~ FS_InterCranial_Vol + Gender + AGE")), data=my.df)$res
    my_res <- lm(as.formula(paste(feat, "~ Gender + AGE")), data=my.df)$res
    #return(feat)
  })
  #my.df2 <- my.df[,fidx]
  rownames(my.df2) <- my.df$Subject

  fam2 <- fam[rownames(my.df2),]
  my_cv <- getFamCVfold(my.df2,10,fam2)
  my_glm <- doubleCV( my.df$CEU > ceu_cut, my.df2, fam="binomial", measure="auc", alpha=0.5, lop="1se", nla = 100, myfold=my_cv)
}
}
