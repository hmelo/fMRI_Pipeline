# RSA from MLM files
# June 8 2016
# Hans Melo

library("oro.nifti")
library("ff")

# Load Mask
tmp.mask.nii <- readNIfTI('/Volumes/HD1/Positivity_Project/masks/r_dlPFC_mask.nii.gz')
mask.vec <- as.vector(tmp.mask.nii@.Data) # graded mask

bin.mask.vec <- (mask.vec>0)*1 # create binary mask

# Load frmi data .nii file
tmp.nii <- readNIfTI('mcdcdn_1216.nii.gz')

brain.dims=c(91,109,91,149)
brain.vec <- matrix(tmp.nii@.Data, nrow=prod(brain.dims[1:3]), ncol=brain.dims[4])

i<-50
j<-50
k<-50
for (i in 1:91){
  for (j in 1:109){
    for (k in 1:91){
      brain_t_1[i,j,k,2:149]<-NaN      
      brain_t_1[i,j,k,2:149]<-tmp.nii@.Data[i,j,k,1:148]
}}}

# try moving vector instead ? double check, also check moved RDM ?

# Apply mask to brain data

data.with.mask <- brain.vec[which(bin.mask.vec==1), ] # extract only values from mask

data.with.mask <- brain.vec*bin.mask.vec # extract all brain values for mask (=0 outside mask)

vec.length <- length(tmp.mask.nii@.Data)

#

# save to .ff file
data.path=getwd()
file.ext=NULL
file.names=NULL
format="nifti"
bitmode="single"
save.path=paste(getwd(),"ffData",sep=.Platform$file.sep)
save.name="data"
mask=FALSE
vectorized=FALSE
output=FALSE
#

brain.ff <- ff(0, dim=c(vec.length, brain.dims[4], 1),
               filename=paste(save.path, "/", save.name, ".ff", sep=""),
               vmode=bitmode)
add(brain.ff, data.with.mask, , , 1)
save(hdr, file=paste(save.path, "/", save.name, "_hdr.RData", sep=""))
ffsave("brain.ff", file=paste(save.path, save.name, sep="/"), rootpath=save.path, compress=FALSE, move=TRUE)


# Compute fMRI data RDM - Representational Dissimilarity Matix

source("faster.cor.r")
fMRI_RSM <- faster.cor(data.with.mask, method="pearson")
heatmap(fMRI_RSM,Colv=NA,Rowv=NA)
fMRI_RDM <- 1-fMRI_RSM # get dissimilarity matrix
heatmap(fMRI_RDM,Colv=NA,Rowv=NA)



fMRI_RDM <- fMRI_RDM[upper.tri(fMRI_RDM)]

## Compute Model RDM : Behavioral data RDM

behav.data = read.csv("/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/setup_headers_rt3.csv")

#### EEG FFTs #####

eeg_ffts.data = read.csv("/Users/admin/Downloads/FFTs.csv")

eeg.data.401 <- eeg_ffts.data[1:224,]
fft.401 <- eeg.data.401[,15:38]

image(fft.401[1])

eeg.vec <- matrix(fft.401, nrow=24, ncol=225)

eeg_RSM <- faster.cor(eeg.vec, method="pearson")
heatmap(eeg_RSM,Colv=NA,Rowv=NA)


Model_Neg600toNeg300AlphaF3<-repmat(fft.401,size(fft.401)[2],1)
B=t(Model_Neg600toNeg300AlphaF3)
Model_Neg600toNeg300AlphaF3_RSM <- Model_Neg600toNeg300AlphaF3!=B
heatmap(fft.401*1,Colv=NA,Rowv=NA)

####################

#Exploit
behav.data.1216 <- behav.data[1:149,]
exploit <- behav.data.1216$exploit
Model_exploit<-repmat(exploit,size(exploit)[2],1)
B=t(Model_exploit)
Model_exploit_RSM <- Model_exploit!=B
heatmap(Model_exploit_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

# Explore (t-1)
exploit_t_1 <- behav.data.1216$exploit_t_1
Model_exploit_t_1<-repmat(exploit_t_1,size(exploit_t_1)[2],1)
B=t(Model_exploit_t_1)
Model_exploit_t_1_RSM <- Model_exploit_t_1!=B
heatmap(Model_exploit_t_1_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

# switch_stay Model
switch_stay <- behav.data.1216$switch_stay
Model_switch_stay<-repmat(switch_stay,size(switch_stay)[2],1)
B=t(Model_switch_stay)
Model_switch_stay_RSM <- Model_switch_stay!=B
heatmap(Model_switch_stay_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean
image(Model_switch_stay_RSM*1)
# Choice_Cat Model
Choice_Cat <- behav.data.1216$Choice_Cat
Model_Choice_Cat<-repmat(Choice_Cat,size(Choice_Cat)[2],1)
B=t(Model_Choice_Cat)
Model_Choice_Cat_RSM <- Model_Choice_Cat!=B
heatmap(1-Model_Choice_Cat_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

# Choice_cat2 Model
Choice_cat_2 <- behav.data.1216$Choice_cat_2
Model_Choice_cat_2<-repmat(Choice_cat_2,size(Choice_cat_2)[2],1)
B=t(Model_Choice_cat_2)
Model_Choice_cat_2_RSM <- Model_Choice_cat_2!=B
heatmap(1-Model_Choice_cat_2_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

# newExplore Model
newExplore <- behav.data.1216$newExplore
Model_newExplore<-repmat(newExplore,size(newExplore)[2],1)
B=t(Model_newExplore)
Model_newExplore_RSM <- Model_newExplore!=B
heatmap(1-Model_newExplore_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean

## correlate fMRI RDM and Model RDM's

# vectorize matrices for rank corrlelation
vfMRI_RSM <- fMRI_RSM[upper.tri(fMRI_RSM)]
vModel_exploit_RSM<- Model_exploit_RSM[upper.tri(Model_exploit_RSM)]

m <- cbind(vfMRI_RSM, vModel_exploit_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_exploit_RSM, method="kendall", use="pairwise") # runs faster

# exploit

vModel_exploit_RSM<- Model_exploit_RSM[upper.tri(Model_exploit_RSM)]

m <- cbind(vfMRI_RSM, vModel_exploit_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_exploit_RSM, method="kendall", use="pairwise") # runs faster


# explore(t-1)

vModel_exploit_t_1_RSM<- Model_exploit_t_1_RSM[upper.tri(Model_exploit_t_1_RSM)]

m <- cbind(vfMRI_RSM, vModel_exploit_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_exploit_t_1_RSM, method="kendall", use="pairwise") # runs faster

# switch_stay

vModel_switch_stay_RSM<- Model_switch_stay_RSM[upper.tri(Model_switch_stay_RSM)]

m <- cbind(vfMRI_RSM, vModel_switch_stay_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_switch_stay_RSM, method="kendall", use="pairwise") # runs faster

# Choice_Cat

vModel_Choice_Cat_RSM<- Model_Choice_Cat_RSM[upper.tri(Model_Choice_Cat_RSM)]

m <- cbind(vfMRI_RSM, vModel_Choice_Cat_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_Choice_Cat_RSM, method="kendall", use="pairwise") # runs faster

# Choice_cat_2

vModel_Choice_cat_2_RSM<- Model_Choice_cat_2_RSM[upper.tri(Model_Choice_cat_2_RSM)]

m <- cbind(vfMRI_RSM, vModel_Choice_Cat_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_Choice_cat_2_RSM, method="kendall", use="pairwise") # runs faster

# newExplore

vModel_newExplore_RSM<- Model_newExplore_RSM[upper.tri(Model_newExplore_RSM)]

m <- cbind(vfMRI_RSM, vModel_newExplore_RSM)
MDS <- cor(m, method="kendall", use="pairwise") 
cor(vfMRI_RSM, 1-vModel_newExplore_RSM, method="kendall", use="pairwise") # runs faster


## Arange in loop for Matrix of RDM correlations

#1. bring all behav matrices in single arrage

Models_behav <- behav.data.1216[,c(8, 10,13,16,18)]

i<-1
Model_sq_mat<-repmat(Models_behav[,i],size(Models_behav[,i])[2],1)
B=t(Model_sq_mat)
Model_RSM <- Model_sq_mat!=B

# initialize data.frame for RSM's
vModel_RSMs<- matrix(zeros(length(Model_RSM[upper.tri(Model_RSM)])*5,1),ncol=length(Models_behav)) # create matrix to store RSMs
vModel_RSMs<- data.frame(vModel_RSMs) # create data frame place holder
names(vModel_RSMs)<-names(Models_behav)

for (i in 1:length(Models_behav)){
  
  # create model RSM  
  Model_sq_mat<-repmat(Models_behav[,i],size(Models_behav[,i])[2],1)
  B=t(Model_sq_mat)
  Model_RSM <- Model_sq_mat!=B
  heatmap(1-Model_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean
  image(Model_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean
  
  vModel_RSMs[[i]]<- 1-Model_RSM[upper.tri(Model_RSM)]*1  # multiply by 1 to make boolean
  
}

# create matrix of all RMS's including fMRI

vModel_RSMs$vfMRI_RSM <- vfMRI_RSM

# create multiple correlation comparison matrix

RDMs_cors <- cor(vModel_RSMs, method="kendall", use="pairwise")
RDMs_cors <- cor(vModel_RSMs[,1:5], method="kendall", use="pairwise")

heatmap(1-RDMs_cors[1:5,1:5],Colv=NA,Rowv=NA)

heatmap(RDMs_cors[1:5,1:5])


image(1-RDMs_cors[1:5,1:5])

## Staritical inference

summary(y<-lm(vfMRI_RSM ~ vModel_exploit_RSM*1))
summary(y<-lm(vfMRI_RSM ~ vModel_Choice_cat_2_RSM*1 + vModel_switch_stay_RSM*1))

summary(y<-lm(vfMRI_RSM ~ vModel_Choice_cat_2_RSM*1))


