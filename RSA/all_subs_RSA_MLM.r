

library("oro.nifti")
library("matlab")
source("faster.cor.r")
## Compute Model RDM : Behavioral data RDM

behav.data = read.csv("/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/setup_headers_rt3.csv")

subs <- unique(behav.data[,1])
subs <- subs[1:71] # last sub is NA so get rid of it

# create/select models

# All subjects RSA and MLM

## Arange in loop for Matrix of RDM correlations

for (i in 1:length(subs)){
  
#   # gets data fro each subject
#   Models_behav <- behav.data[behav.data[,1]==subs[i],c(1,8, 10,13,16,18)]
#   Models_behav <-behav.data.sub[1:149,]
  # get data for each subject
  Models_behav <- behav.data[behav.data[,1]==subs[i],c(1,3,5:31)]
  Models_behav <-Models_behav[1:149,]
  Models_behav$explore<-(Models_behav$exploit!=1)*1
  Models_behav$Switch_to_Explore<-c(zeros(149,1)*NA)
  #Models_behav$Switch_to_Explore<-array(0,c(1,149))*NA
  Models_behav$Switch_to_Explore_tp1<-c(zeros(149,1)*NA)
  #Models_behav$Switch_to_Explore_tp1<-array(0,c(1,149))*NA
  Models_behav$Explore_Same<-c(zeros(1,149)*NA)
  #Models_behav$Explore_Same<-array(0,c(1,149))*NA
    
  # create/select models
  for (trial in 2:149){
    # Switch to Explore
    if (Models_behav$explore[trial-1]==0 && Models_behav$explore[trial]==1){
      Models_behav$Switch_to_Explore[trial]<-1
    } else {Models_behav$Switch_to_Explore[trial]<-0}
    # Switch to explore t+1
    if (trial < 149){
    if (Models_behav$explore[trial]==0 && Models_behav$explore[trial+1]==1){
      Models_behav$Switch_to_Explore_tp1[trial]<-1
    } else {Models_behav$Switch_to_Explore_tp1[trial]<-0}}
    # Explore_Same/Explore_different
    if (!is.na(Models_behav$choice[trial]) &&  !is.na(Models_behav$choice[trial-1])){
      if (Models_behav$explore[trial]==1 && Models_behav$choice[trial]==Models_behav$choice[trial-1]){
        Models_behav$Explore_Same[trial]<-1} 
      if (Models_behav$explore[trial]==1 && Models_behav$choice[trial]!=Models_behav$choice[trial-1]){
        Models_behav$Explore_Same[trial]<-0
      }else {Models_behav$Explore_Same[trial]<-NA}
    }
  }
  
  # initialize Model_RSM
k<-1
#Model_sq_mat<-repmat(Models_behav[,k],size(Models_behav[,k])[2],1)
Model_sq_mat<-repmat(Models_behav[,k],size(Models_behav[,k])[2],1)

B=t(Model_sq_mat)
Model_RSM <- Model_sq_mat!=B

# initialize data.frame for RSM's
vModel_RSMs<- matrix(zeros(length(Model_RSM[upper.tri(Model_RSM)])*length(Models_behav),length(subs)),ncol=length(Models_behav)) # create matrix to store RSMs
vModel_RSMs<- data.frame(vModel_RSMs) # create data frame place holder
names(vModel_RSMs)<-names(Models_behav)

vlength <- size(vModel_RSMs[,1])[2]/length(subs)
sub_vec <- c(repmat(subs,vlength, 1))
vModel_RSMs$subject <- sub_vec

# Create Model (behavioral) RSM

for (j in 2:size(vModel_RSMs)[2]){
  #print(j)
  
  # create model RSMs
  Model_sq_mat<-repmat(Models_behav[,j],size(Models_behav[,j])[2],1)
  B=t(Model_sq_mat)
  Model_RSM <- Model_sq_mat!=B
#   heatmap(1-Model_RSM*1,Colv=NA,Rowv=NA) # create heatmap multiply by 1 to make boolean
#   image(Model_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

  vModel_RSMs[j]<- 1-Model_RSM[upper.tri(Model_RSM)]*1  # multiply by 1 to make boolean
  vModel_RSMs[vModel_RSMs[,1]==subs[i],j] <- 1-Model_RSM[upper.tri(Model_RSM)]*1

# create model RDMs t-1

}
}

# create matrix of all RMS's including fMRI

### Load fmri data

# Load Mask
tmp.mask.nii <- readNIfTI('/Volumes/HD1/Positivity_Project/masks/r_dlPFC_mask.nii.gz')
mask.vec <- as.vector(tmp.mask.nii@.Data) # graded dlpfc mask

bin.mask.vec <- (mask.vec>0)*1 # create binary mask

# Load frmi data .nii file

for( k in 1:length(subs)){
  
print(sprintf("Creating vectorized RDM for %s.", subs[k]))

tmp.nii <- readNIfTI(sprintf('mcdcdn_%s.nii.gz',subs[k]))
brain.dims=c(91,109,91,149)
brain.vec <- matrix(tmp.nii@.Data, nrow=prod(brain.dims[1:3]), ncol=brain.dims[4])

## RSA brain t0

# Apply mask to brain data
data.with.mask <- brain.vec[which(bin.mask.vec==1), ] # extract only values from mask FASTER
fMRI_RSM <- faster.cor(data.with.mask, method="pearson")
#heatmap(fMRI_RSM,Colv=NA,Rowv=NA) # create heatmap
fMRI_RDM <- 1-fMRI_RSM # get dissimilarity matrix

#data.with.mask <- brain.vec*bin.mask.vec # extract all brain values for mask (=0 outside mask) # SLOWER (takes full brain matrix)
#fMRI_RSM <- faster.cor(data.with.mask, method="pearson")
#heatmap(fMRI_RSM,Colv=NA,Rowv=NA) # create heatmap
#fMRI_RDM <- 1-fMRI_RSM # get dissimilarity matrix

#vec.length <- length(tmp.mask.nii@.Data)

# RSA brain t-2 and t-2
first_trial<- (0:902628)*149+1
second_trial<-first_trial+1
dlPFC_vec_t_1<-NaN
dlPFC_vec_t_1[2:134491721] <- brain.vec[1:134491720]
dlPFC_vec_t_1[first_trial]<- NaN
dlPFC_vec_t_1 <- matrix(dlPFC_vec_t_1, nrow=prod(brain.dims[1:3]), ncol=brain.dims[4])

dlPFC_vec_t_1_mask <- dlPFC_vec_t_1[which(bin.mask.vec==1), ] # extract only values from mask FASTER
fMRI_t_1_RSM <- cor(dlPFC_vec_t_1_mask, method="pearson", use = "pairwise.complete.obs") #### WORKS for NaNs! 
#heatmap(fMRI_t_1_RSM,Colv=NA,Rowv=NA) # create heatmap
fMRI_t_1_RDM <- 1-fMRI_t_1_RSM # get dissimilarity matrix

dlPFC_vec_t_2<-NaN
dlPFC_vec_t_2[2:134491721] <- dlPFC_vec_t_1[1:134491720]
dlPFC_vec_t_2[first_trial]<- NaN
dlPFC_vec_t_2 <- matrix(dlPFC_vec_t_2, nrow=prod(brain.dims[1:3]), ncol=brain.dims[4])

dlPFC_vec_t_2_mask <- dlPFC_vec_t_2[which(bin.mask.vec==1), ] # extract only values from mask FASTER
fMRI_t_2_RSM <- cor(dlPFC_vec_t_2_mask, method="pearson", use = "pairwise.complete.obs") #### WORKS for NaNs! 
#heatmap(fMRI_t_2_RSM,Colv=NA,Rowv=NA) # create heatmap
fMRI_t_2_RDM <- 1-fMRI_t_2_RSM # get dissimilarity matrix

# Ad fMRI to vModel_RSMs
vModel_RSMs$vfMRI_RSM[vModel_RSMs[,1]==subs[k]] <- fMRI_RSM[upper.tri(fMRI_RSM)]
vModel_RSMs$vfMRI_t_1_RSM[vModel_RSMs[,1]==subs[k]] <- fMRI_t_1_RSM[upper.tri(fMRI_RSM)]
vModel_RSMs$vfMRI_t_2_RSM[vModel_RSMs[,1]==subs[k]] <- fMRI_t_2_RSM[upper.tri(fMRI_RSM)]

}

# write vModel to csv

write.csv(vModel_RSMs, file = "vModel_RSMs_t_1_t_2.csv")

# brain RSMs

vBrain_RSMs <- vModel_RSMs$vfMRI_RSM

vModel_RSMs$vfMRI_RSM <- vBrain_RSMs
# create multiple correlation comparison matrix

# RDMs_cors <- cor(vModel_RSMs, method="kendall", use="pairwise")
# RDMs_cors <- cor(vModel_RSMs[,1:5], method="kendall", use="pairwise")

heatmap(1-RDMs_cors[1:5,1:5],Colv=NA,Rowv=NA)
heatmap(RDMs_cors[1:5,1:5])
image(1-RDMs_cors[1:5,1:5])

## Staritical inference


##
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy + switch_stay+ (1|subject),data=vModel_RSMs))

summary(lmer_mri<- lmer(switch_stay ~ vfMRI_RSM + vfMRI_t_1_RSM + vfMRI_t_2_RSM + (1|subject),data=vModel_RSMs))
summary(lmer_mri<- lmer(explore ~ vfMRI_RSM + vfMRI_t_1_RSM + (1|subject),data=vModel_RSMs))
summary(lmer_mri<- lmer(Switch_Policy ~ vfMRI_t_1_RSM + vfMRI_t_2_RSM + (1|subject),data=vModel_RSMs))


###

# Anslysis 4:

summary(lmer_mri1<- lmer(vfMRI_RSM ~ exploit + (1|subject),data=vModel_RSMs))# **
summary(lmer_mri2<- lmer(vfMRI_RSM ~ switch_stay + (1|subject),data=vModel_RSMs)) # ***
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy + (1|subject),data=vModel_RSMs))

summary(lmer_mri1<- lmer(vfMRI_RSM ~ exploit_tp1 + (1|subject),data=vModel_RSMs)) # ***
summary(lmer_mri2<- lmer(vfMRI_RSM ~ switch_stay_tp1 + (1|subject),data=vModel_RSMs)) # ***
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy_tp1 + (1|subject),data=vModel_RSMs))


summary(lmer_mri<- lmer(vfMRI_RSM ~ Choice_Cat + (1|subject),data=vModel_RSMs)) #**
summary(lmer_mri<- lmer(vfMRI_RSM ~ Choice_cat_2 + (1|subject),data=vModel_RSMs)) # ***
summary(lmer_mri<- lmer(vfMRI_RSM ~ newExplore + (1|subject),data=vModel_RSMs)) 
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_to_Explore + (1|subject),data=vModel_RSMs)) 
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy + (1|subject),data=vModel_RSMs)) 
summary(lmer_mri<- lmer(vfMRI_RSM ~ Explore_Same + (1|subject),data=vModel_RSMs)) # not enough DF

# Multiple predictors
summary(lmer_mri<- lmer(vfMRI_RSM ~ exploit + exploit_tp1 + (1|subject),data=vModel_RSMs)) # * , ***
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy + Switch_Policy_tp1 + (1|subject),data=vModel_RSMs)) 
summary(lmer_mri<- lmer(vfMRI_RSM ~ Switch_Policy + exploit + (1|subject),data=vModel_RSMs)) # , **
summary(lmer_mri<- lmer(vfMRI_RSM ~ switch_stay + exploit + (1|subject),data=vModel_RSMs)) # ***,

summary(lmer_mri<- lmer(vfMRI_RSM ~ exploit + exploit_tp1 + Switch_Policy + Switch_Policy_tp1 + (1|subject),data=vModel_RSMs)) # ***,

