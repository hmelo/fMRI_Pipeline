library("lme4", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("lmerTest", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library(effects)

setwd("/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/")

mydata = read.csv("setup_headers.csv")

fmriDATA <- read.table("/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/Decision_SV_gene_k/Clust_02_t5_th3_m4_m4_m12.1D", header = FALSE)
fmri_seed <- read.table("/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/Decision_SV_gene_k/Clust_01_t5_th3_2_m44_38.1D", header = FALSE)

WM_COV = read.table("/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/all_mcdcdn_WM_confound_meants.txt", header = FALSE)
CSF_COV = read.table("/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/all_mcdcdn_CSF_confound_meants.txt", header = FALSE)


fullDAT <- cbind(mydata,fmriDATA, fmri_seed,WM_COV, CSF_COV)

fullDAT$decision_01 <- (fullDAT$decision + 1)/2
fullDAT$fMRI <- fmriDATA$V1/1000
fullDAT$fMRI_means <- fmri_seed$V1/1000
attach(fullDAT)

for (i in 1:length(fullDAT$COMT_rs4680_alleles)) {
  if (fullDAT$COMT_rs4680_alleles[i]==1){
    fullDAT$COMT[i]='Val/Val'} else if (fullDAT$COMT_rs4680_alleles[i]==2){
      fullDAT$COMT[i]='Met/Met'
    }else if (fullDAT$COMT_rs4680_alleles[i]==3){
      fullDAT$COMT[i]='Val/Met'
    }
  if (fullDAT$decision_01[i]==.5){
    fullDAT$decision_01[i]=NA}
  
  if (fullDAT$COMT_rs4680_bin[i]==0){
    fullDAT$COMT_bin[i]='Val/Val'
  }else if (fullDAT$COMT_rs4680_bin[i]==1){
    fullDAT$COMT_bin[i]='Met-'
  }else{
    fullDAT$COMT_bin[i]=NA
  }
  
  
}


thedata<-data.frame(Decision=fullDAT$decision_01,Genotype=fullDAT$COMT,COMT_bin=fullDAT$COMT_bin, k=fullDAT$k_estimate*1000, SV=fullDAT$subj_value_chosen)
thedata$fMRI <- fmriDATA$V1/1000
thedata$fMRI_seed <- fmri_seed$V1/1000
thedata$WM_COV <- WM_COV$V1/1000
thedata$CSF_COV <- CSF_COV$V1/1000
thedata$SubNo<-SubNo

#fMRI_data ~  (fMRI_seed * as.factor(COMT_rs4680_bin)) + WM_COV + CSF_COV + (1|subject)

summary(lmer_mri<- lmer(fMRI ~  (fMRI_seed * Genotype) + WM_COV + CSF_COV + (1|SubNo),data=thedata))
