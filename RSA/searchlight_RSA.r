# Search light RSA Analysis


##Load Masks/ROIs
#Image Parameters
img <- list()
vxlsz <- c(2,2,2)
voxcount <- 0
searchradius <- 5   # in mm
alpha <- 0.05
fdr <- TRUE #Logical whether or not to perform FDR correction on results, will return uncorrected and FDR corrected p values

# Load Mask
tmp.mask.nii <- readNIfTI('/Volumes/HD1/Positivity_Project/masks/r_dlPFC_mask.nii.gz')
mask <- array(as.numeric(tmp.mask.nii@.Data > 0), dim=brain.dims[1:3])
mask.vec <- as.vector(tmp.mask.nii@.Data) # graded mask
bin.mask.vec <- (mask.vec>0)*1 # create binary mask
nvoxels <- sum(bin.mask.vec) # number of voxels in mask

#Analysis Parameters
analysis <- list()
analysis$preds.nointeractions <- c("exploit")
analysis$preds <- c("exploit") # Possible values: response, win, loss, value, expWin, expLoss, expVal
analysis$npreds <- length(analysis$preds)

##Specify Model 1
####-To be run on every voxel (within mask volume)
####-Initialize data containers
mdl1 <- list()
mdl1$formula <- as.formula("brain ~ exploit + (1|subject)")
mdl1$results <- list()
mdl1$results$voxel = matrix(0, nrow=nvoxels, ncol=3)
mdl1$results$betas = matrix(0, nrow=nvoxels, ncol=5)
mdl1$results$p     = matrix(0, nrow=nvoxels, ncol=5)
mdl1$results$r2    = vector("numeric", length=nvoxels)

##Specify Model 2
mdl2 <- list()
mdl2$formula <- as.formula("brain ~ switch_stay + (1|subject)")
mdl2$results <- list()
mdl2$results$voxel = matrix(0, nrow=nvoxels, ncol=3)
mdl2$results$betas = matrix(0, nrow=nvoxels, ncol=4)
mdl2$results$p     = matrix(0, nrow=nvoxels, ncol=4)
mdl2$results$r2    = vector("numeric", length=nvoxels)

mdlcomp <- list()
mdlcomp$comp1v2 <- list()
mdlcomp$comp1v2$dr2    <- matrix(0, nrow=nvoxels, ncol=1)
mdlcomp$comp1v2$AIC    <- matrix(0, nrow=nvoxels, ncol=2)
mdlcomp$comp1v2$BIC    <- matrix(0, nrow=nvoxels, ncol=2)
mdlcomp$comp1v2$loglik <- matrix(0, nrow=nvoxels, ncol=2)
mdlcomp$comp1v2$chi2   <- matrix(0, nrow=nvoxels, ncol=1)
mdlcomp$comp1v2$p      <- matrix(0, nrow=nvoxels, ncol=1)



##Load Brain Data
brainff <- list()
brain.dims <- c(91,109,91,149)
for( k in 1:length(subs)){
  
  tmp.nii <- readNIfTI(sprintf('mcdcdn_%s.nii.gz',subs[k]))
  
  
  brain.vec <- matrix(tmp.nii@.Data, nrow=prod(brain.dims[1:3]), ncol=brain.dims[4])
  
#   # Apply mask to brain data
#   data.with.mask <- brain.vec[which(bin.mask.vec==1), ] # extract only values from mask
#   data.with.mask <- brain.vec*bin.mask.vec # extract all brain values for mask (=0 outside mask)
#   vec.length <- length(tmp.mask.nii@.Data)
#     
  brainff[[k]] <- ff(brain.vec, dim=c(brain.dims[1]*brain.dims[2]*brain.dims[3], ncol=brain.dims[4]), filename=sprintf("%s.ff", subs[k]))
}

rm("brain.vec")

##Load Predictors

for (i in 1:length(subs)){
  
  #   # gets data fro each subject
  #   Models_behav <- behav.data[behav.data[,1]==subs[i],c(1,8, 10,13,16,18)]
  #   Models_behav <-behav.data.sub[1:149,]
  # get data for each subject
  Models_behav <- behav.data[behav.data[,1]==subs[i],c(1,3,5:31)]
  Models_behav <-Models_behav[1:149,]
  Models_behav$explore<-(Models_behav$exploit!=1)*1
  Models_behav$Switch_to_Explore<-c(zeros(149,1)*NA)
  Models_behav$Switch_to_Explore_tp1<-c(zeros(149,1)*NA)
  Models_behav$Explore_Same<-c(zeros(1,149)*NA)
  
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
  
  for (j in 2:size(vModel_RSMs)[2]){
    print(j)
    
    # create model RSM  
    Model_sq_mat<-repmat(Models_behav[,j],size(Models_behav[,j])[2],1)
    B=t(Model_sq_mat)
    Model_RSM <- Model_sq_mat!=B
    #   heatmap(1-Model_RSM*1,Colv=NA,Rowv=NA) # create heatmap multiply by 1 to make boolean
    #   image(Model_RSM*1,Colv=NA,Rowv=NA) # multiply by 1 to make boolean
    
    vModel_RDMs[j]<- 1-Model_RSM[upper.tri(Model_RSM)]*1  # multiply by 1 to make boolean
    vModel_RDMs[vModel_RSMs[,1]==subs[i],j] <- 1-Model_RSM[upper.tri(Model_RSM)]*1
    
  }
}

##Run Analysis
source("getVoxelSphere.r")
for (iter in 454032:prod(brain.dims[1:3])) {
  t0 = Sys.time()
  # Get Voxel Index
  idx <- vectorIndex2arrayIndex(iter, dim=brain.dims[1:3], dimorder=1:3)
  x <- idx[1]
  y <- idx[2]
  z <- idx[3]
  
  if (mask[x,y,z] == 1) { # Check if brain in voxel, if YES analyze
    roi <- getVoxelSphere(searchradius, vxlsz, c(x,y,z), brain.dims) # Get Searchlight Voxels
    searchlight <- array(0, dim=dim(mask)) # Initialize searchlight mask
    searchlight[roi] <- 1                  # Set voxels within searchlight to TRUE
    searchlight <- searchlight * mask      # Remove non-brain voxels from searchlight
    searchlight <- as.vector(searchlight)  # Vectorize searchlight
    
    if (sum(searchlight) > 4) {
      voxcount <- voxcount + 1
      mdl1$results$voxel[voxcount, ] <- c(x, y, z) # Save Voxel Coordinates
      mdl2$results$voxel[voxcount, ] <- c(x, y, z) # Save Voxel Coordinates
      
      brainrdm <- vector("numeric", length=0)     # Initialize RDM vector
      
      for (i in 1:length(subs)) {                 # Calculate RDMs, looping through each subject
        mx <- brainff[[i]][which(searchlight==1), ] # Retrieve brain data within searchlight
        # xy
        # x^2
        # y^2
        # colSums x, y, xy, x^2, y^2
        # r = (n(Exy)-ExEy) / sqrt((nEx^2 - (Ex)^2)*(nEy^2 - (Ey)^2))
        rtemp <- cor(mx)
        rtemp <- 1-rtemp                     # Convert to representational disimilarity matrix
        brainrdm <- c(brainrdm, rtemp[upper.tri(rtemp)])    # contains all brainrdms for all subjects
      }
      
      # Primary Model
      temp.df <- data.frame(brain=brainrdm, vModel_RDMs)
#       temp.df <- na.omit(temp.df) #remove unsampled voxels ### CAUSES MLM TO CRASH
      lme.mdl1.summary <- summary(lme.mdl1 <- lmer(mdl1$formula, data=temp.df)) # Summary of Linear Mixed-effects Model

      # Add data from LME to output variables
      mdl1$results$betas[voxcount, ] <- lme.mdl1.summary$coefficients[(2:(1+npreds)),1]
      mdl1$results$p[voxcount, ]     <- lme.mdl1.summary$coefficients[(2:(1+npreds)),5]
#       mdl1$results$r2[voxcount]      <- summary(lm(brainrdm ~ fitted(lme.mdl1)))$r.squared # crashes! ?
      
      # Model 2
      
      lme.mdl2.summary <- summary(lme.mdl2 <- lmer(mdl2$formula, data=temp.df)) # Summary of Linear Mixed-effects Model
      
      # Add data from LME to output variables
      mdl2$results$betas[voxcount, ] <- lme.mdl2.summary$coefficients[(2:(1+npreds)),1]
      mdl2$results$p[voxcount, ]     <- lme.mdl2.summary$coefficients[(2:(1+npreds)),5]
#       mdl2$results$r2[voxcount]      <- summary(lm(brainrdm ~ fitted(lme.mdl2)))$r.squared # crashes!
      

      # Calculate correlations btween RDMs
      
      
      #Compare LME Models
      mdl.aov <- anova(lme.mdl1, lme.mdl2)
#       mdlcomp$comp1v2$dr2[voxcount]      <- mdl2$results$r2[voxcount] - mdl1$results$r2[voxcount] # crashes
      mdlcomp$comp1v2$AIC[voxcount, ]    <- mdl.aov$AIC
      mdlcomp$comp1v2$BIC[voxcount, ]    <- mdl.aov$BIC
      mdlcomp$comp1v2$loglik[voxcount, ] <- mdl.aov$logLik
      mdlcomp$comp1v2$chi2[voxcount]     <- mdl.aov$Chisq[2]
      mdlcomp$comp1v2$p[voxcount]        <- mdl.aov$P[2]
   
    }
  }
  t2 <- Sys.time() - t0
  print(iter)
#   print(sprintf("x: %0.0f, y: %0.0f, z: %0.0f, Complete: %0.4f%%, Elapsed Time: %0.4f", x, y, z, voxcount/nvoxels*100, t2))
}