library("lme4", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("lmerTest", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library(effects)

#### DARP32 ####

### Cluster 1
mydata = read.csv("setup_headers.csv")

fmriDATA <- read.table("/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/MPI_version/DARPP32_fmri_mean/Clust_02_mean_24m56_10_207voxels_OFC.1D", header = FALSE)
#fmriCOV <- read.table("/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/despiked_decision_bin_all_MEANS/Clust_05_meanCOV_m8_m4_42.1D", header = FALSE)
fullDAT <- cbind(mydata,fmriDATA)

fullDAT$explore <- (fullDAT$exploit -1)*-1 # so that exploration=1 , exploit=0
fullDAT$fMRI <- fmriDATA$V1/1000
fullDAT$fMRI_means <- fmriCOV$V1/1000

fullDAT$DARPP32 <- fullDAT$DARPP32_rs907094_bin

fullDAT$DARPP32[fullDAT$DARPP32==0]<-'C/C'
fullDAT$DARPP32[fullDAT$DARPP32==1]<-'T-carrier'
attach(fullDAT)


thedata<-data.frame(explore=fullDAT$explore,Genotype=fullDAT$DARPP32)
thedata$fMRI <- fmriDATA$V1/1000
thedata$fMRI_means <- fmriCOV$V1/1000
thedata$subject<-fullDAT$subject

summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|subject),data=thedata))

summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|subject),family=binomial,data=thedata))
summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))

# Random intercept and random slope for fMRI
summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1+fMRI|Genotype/subject),family=binomial,data=thedata, control=glmerControl(optimizer="bobyqa")))
##### USE MODEL BELOW
summary(glmer_mri<- glmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))
#####
plot(effect(term="fMRI:Genotype",mod=glmer_mri,default.levels=10579),multiline=TRUE,ci.style="bands",rescale.axis=FALSE)

summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))

## Calculating semi-partial R^2
# option 1
summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|subject),data=thedata))
anova.summary <- anova(lmer_mri); anova.summary
#effect.sizes <- with(anova.summary, NumDF/DenDF*F.value/(1+NumDF/DenDF*F.value))
#names(effect.sizes) <- row.names(anova.summary)
#effect.sizes

# option 2, get numDf and DenDf form above and compute semi-partial R^2
summary(glmer_mri<- glmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))
F_value <- anova(glmer_mri)$"F value"[3]; anova.summary  # use to get F value

effect.sizes <- with(anova.summary, NumDF/DenDF*F_value/(1+NumDF/DenDF*F_value))
names(effect.sizes) <- row.names(anova.summary)
effect.sizes

## Calculating Pseudo-R^2 at each level

#### Baseline Model to calculate Pseudo-R^2 and ICC ####
mixed.baseline <- glmer( explore ~ (1|Genotype/subject),family=binomial,data=thedata); summary(mixed.baseline) #Using the lmer() function

## Calculate ICC : determine whether data within schools are correlated with each other (by calculating ICC)
# variance of intercept = 0.6223   use pi^2/3 for variance of residual = 3.29 

icc = (0.6223)/ (0.6223 +3.29); icc

#The following syntax will calculate the harmonic mean of the number of legislators in the state, which you can use to calculate level 2 pseudo R^2
average.runs.state <- with( thedata, round(harmonic.mean(as.numeric(subject)))); average.runs.state
n2 <- average.runs.state;n2 #Harmonic mean

average.genotype.state <- with( thedata, round(harmonic.mean(as.numeric(Genotype)))); average.genotype.state
n3 <- average.genotype.state;n3 #Harmonic mean


#### Calculate Pseudo R^2 at levels 1 and 2, use lmer to get variance ?????? Probably NOT Right
# Effect size for the whole model
summary(glmer_mri)
summary(mixed.baseline)
# level 1 pseudo R^2
# R2McF = 1 – ln(LM) / ln(L0)
# R1^2 = 1 - (intercept.variance.model1 + residual.variance.model1)/(intercept.variance.baseline + residual.variance.baseline)
pseudo.r1.sq <- 1 - (0.6275 + 3.29) /(.6223 +3.29); pseudo.r1.sq 

pseudo.r1.sq <- 1 - (as.numeric(VarCorr(lmer_mri)[1,1])+as.numeric(VarCorr(lmer_mri)[2,1])) / (as.numeric(VarCorr(mixed.baseline)[1,1])+as.numeric(VarCorr(mixed.baseline)[2,1]))

#level 2 pseudo R^2
# R2^2 = 1 - (intercept.variance.model1 + residual.variance.model1/n)/(intercept.variance.baseline + residual.variance.baseline/n)
pseudo.r2.sq <- 1 - (0 +3.29) /(1.325e-10  +3.29); pseudo.r2.sq 
pseudo.r2.sq <- 1 - (0 +3.29/1345) /(1.325e-10  +3.29/1345); pseudo.r2.sq 

pseudo.r2.sq <- 1 - (as.numeric(VarCorr(mixed)[1,1])+as.numeric(VarCorr(mixed)[3,1])/n) / (as.numeric(VarCorr(mixed.baseline)[1,1])+as.numeric(VarCorr(mixed.baseline)[2,1])/n)

#level 3 pseudo R^2
pseudo.r3.sq <- 1 - (0 +3.29/2) /(1.325e-10  +3.29/2); pseudo.r3.sq 

### Cluster 2
mydata = read.csv("setup_headers.csv")
fmriDATA <- read.table("/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/MPI_version/DARPP32_fmri_mean/Clust_01_mean_m20m6210_214voxels_OFC.1D", header = FALSE)
fullDAT <- cbind(mydata,fmriDATA)

fullDAT$explore <- (fullDAT$exploit -1)*-1 # so that exploration=1 , exploit=0
fullDAT$fMRI <- fmriDATA$V1/1000
fullDAT$fMRI_means <- fmriCOV$V1/1000

fullDAT$DARPP32 <- fullDAT$DARPP32_rs907094_bin

fullDAT$DARPP32[fullDAT$DARPP32==0]<-'C/C'
fullDAT$DARPP32[fullDAT$DARPP32==1]<-'T-carrier'
attach(fullDAT)


thedata<-data.frame(explore=fullDAT$explore,Genotype=fullDAT$DARPP32)
thedata$fMRI <- fmriDATA$V1/1000
thedata$fMRI_means <- fmriCOV$V1/1000
thedata$subject<-fullDAT$subject

summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|subject),data=thedata))

summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|subject),family=binomial,data=thedata))
summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))
summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1+fMRI|Genotype/subject),family=binomial,data=thedata, control=glmerControl(optimizer="bobyqa")))
##### USE MODEL BELOW
summary(glmer_mri<- glmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))
#####
plot(effect(term="fMRI:Genotype",mod=glmer_mri,default.levels=10579),multiline=TRUE,ci.style="bands",rescale.axis=FALSE)

summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))

## Calculating semi-partial R^2
# option 1
summary(lmer_mri<- lmer(explore ~  fMRI * Genotype + (1|subject),data=thedata))
anova.summary <- anova(lmer_mri); anova.summary
#effect.sizes <- with(anova.summary, NumDF/DenDF*F.value/(1+NumDF/DenDF*F.value))
#names(effect.sizes) <- row.names(anova.summary)
#effect.sizes

# option 2, get numDf and DenDf form above and compute semi-partial R^2
summary(glmer_mri<- glmer(explore ~  fMRI * Genotype + (1|Genotype/subject),family=binomial,data=thedata))
F_value <- anova(glmer_mri)$"F value"[3]; anova.summary  # use to get F value

effect.sizes <- with(anova.summary, NumDF/DenDF*F_value/(1+NumDF/DenDF*F_value))
names(effect.sizes) <- row.names(anova.summary)
effect.sizes

## Calculating Pseudo-R^2 at each level

#### Baseline Model to calculate Pseudo-R^2 and ICC ####
mixed.baseline <- glmer( explore ~ (1|Genotype/subject),family=binomial,data=thedata); summary(mixed.baseline) #Using the lmer() function

## Calculate ICC : determine whether data within schools are correlated with each other (by calculating ICC)
# variance of intercept = 0.02546   variance of residual = 0.19295 

icc = (0.6293)/ (0.6275 +3.29); icc

#The following syntax will calculate the harmonic mean of the number of legislators in the state, which you can use to calculate level 2 pseudo R^2
average.runs.state <- with( thedata, round(harmonic.mean(as.numeric(subject)))); average.runs.state
n2 <- average.runs.state;n2 #Harmonic mean

average.genotype.state <- with( thedata, round(harmonic.mean(as.numeric(Genotype)))); average.genotype.state
n3 <- average.genotype.state;n3 #Harmonic mean

#### Calculate Pseudo R^2 at levels 1 and 2, use lmer to get variance ?????? NOT Right
# Effect size for the whole model
summary(glmer_mri)
summary(mixed.baseline)
# level 1 pseudo R^2
# R2McF = 1 – ln(LM) / ln(L0)
# R1^2 = 1 - (intercept.variance.model1 + residual.variance.model1)/(intercept.variance.baseline + residual.variance.baseline)
pseudo.r1.sq <- 1 - (0.6293 + 3.29) /(.6223 +3.29); pseudo.r1.sq 

pseudo.r1.sq <- 1 - (as.numeric(VarCorr(lmer_mri)[1,1])+as.numeric(VarCorr(lmer_mri)[2,1])) / (as.numeric(VarCorr(mixed.baseline)[1,1])+as.numeric(VarCorr(mixed.baseline)[2,1]))

#level 2 pseudo R^2
# R2^2 = 1 - (intercept.variance.model1 + residual.variance.model1/n)/(intercept.variance.baseline + residual.variance.baseline/n)
pseudo.r2.sq <- 1 - (0 +3.29) /(1.325e-10  +3.29); pseudo.r2.sq 
pseudo.r2.sq <- 1 - (0 +3.29/1345) /(1.325e-10  +3.29/1345); pseudo.r2.sq 

pseudo.r2.sq <- 1 - (as.numeric(VarCorr(mixed)[1,1])+as.numeric(VarCorr(mixed)[3,1])/n) / (as.numeric(VarCorr(mixed.baseline)[1,1])+as.numeric(VarCorr(mixed.baseline)[2,1])/n)

#level 3 pseudo R^2
pseudo.r3.sq <- 1 - (0 +3.29/2) /(1.325e-10  +3.29/2); pseudo.r3.sq 

## compare the two models
# find AIC for the models
AIC( glmer_mri )
AIC( lmer_mri )

# relative likelihood of random slope model (i.e., the one with the higher AIC)
exp( ( AIC( random.intercept ) - AIC( random.slope ) ) / 2 )
# value (.0168) is very small, so choose random intercept model

# find the BIC for the models
BIC( glmer_mri )
BIC( lmer_mri )
# calculate difference between the models' BIC values
abs( BIC( glmer_mri ) - BIC( lmer_mri ) )
# value is 8.18, indicating that the models "probably fit differently", therefore choose the model with the lower BIC (random intercept)



#### COMT ####

mydata = read.csv("setup_headers.csv")

fmriDATA <- read.table("/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/MPI_version/COMT_mean/Clust_01_mean_p4_26_4_46_2582voxels_cingulate.1D", header = FALSE)
#fmriCOV <- read.table("/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/despiked_decision_bin_all_MEANS/Clust_05_meanCOV_m8_m4_42.1D", header = FALSE)
fullDAT <- cbind(mydata,fmriDATA)

fullDAT$explore <- (fullDAT$exploit -1)*-1
fullDAT$fMRI <- fmriDATA$V1/1000
fullDAT$fMRI_means <- fmriCOV$V1/1000

fullDAT$COMT <- fullDAT$COMT_rs4680_bin

fullDAT$COMT[fullDAT$COMT==0]<-'Val/Val'
fullDAT$COMT[fullDAT$COMT==1]<-'Met-carrier'
#attach(fullDAT)


thedata<-data.frame(explore=fullDAT$explore,Genotype=fullDAT$COMT)
thedata$fMRI <- fmriDATA$V1/1000
#thedata$fMRI_means <- fmriCOV$V1/1000
thedata$subject<-fullDAT$subject

#summary(lmer_mri<- lmer(explore ~ (fMRI ) * Genotype + (1|subject),data=thedata))
summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|subject),family=binomial,data=thedata))

plot(effect(term="fMRI:Genotype",mod=lmer_mri,default.levels=10579),grid=TRUE, multiline=TRUE,ci.style="bands",rescale.axis=FALSE)

summary(lmer_mri<- glmer(explore ~ fMRI * Genotype + (1|subject),family=binomial,data=thedata))

#### Genetic effects on behavior
summary(lmer_mri<- lmer(expl_par ~ Genotype + (1|subject),data=thedata))

clean_expl_par<-expl_par[expl_par<50]
clean_COMT<- COMT[expl_par<50]
summary(clean_expl_par)

mean(clean_expl_par[clean_COMT=="Val/Val"],na.rm=TRUE)
mean(clean_expl_par[clean_COMT=="Met-carrier"],na.rm=TRUE)

anova(lm( clean_expl_par ~ clean_COMT)) # yes
boxplot( clean_expl_par ~ clean_COMT, ylab="expl_par", xlab="Genotype")




expl_par