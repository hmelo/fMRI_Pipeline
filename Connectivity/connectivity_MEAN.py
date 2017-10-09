### functional connectivity ###

# does functional connectivity for the betas i.e. task specific!
# need beta time series, and seed (mask) mean time series 
# also needs setup.csv, and white matter (WM) and CSF covariate mean time series.

import numpy 
import nibabel as nib
import sys

from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

#from numpy import genfromtxt
#from numpy import tile

r = robjects.r
r.library('lme4')
stats = importr('stats')
base = importr('base')

print 'Loading mask...'
maskImage = nib.load('/Users/Home/Google Drive/Positivity Project/MLM test/hlmflexfmri/mask.nii.gz')
mask = maskImage.get_data()

# load beta fmri (deconvolved fmri); USE MEAN_CENTERED BETAS!
print 'Loading centered fmri data ...'
pp01 = nib.load('/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/mcdn_functional_merged.nii')
s01 = pp01.get_data() 

# load ave_seed_beta fmri merged, for all subjects; USE MEAN_CENTERED SEED!
print 'Loading mean seed fmri data ...'
seed = numpy.genfromtxt('/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/mPFC_seed.txt', delimiter=",")
seed_fmri=FloatVector(seed[:,0])
robjects.globalenv["seed_fmri"] = seed_fmri

# load CSF and WM covariates; USE MEAN_CENTERED WM and CSF COVARIATES!
print 'Loading CSF and WM covariates betas ...'
wm = numpy.genfromtxt('/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/all_WM_cov.txt', delimiter=",")
WM_cov=FloatVector(wm[:,0])
csf = numpy.genfromtxt('/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/all_CSF_cov.txt', delimiter=",")
CSF_cov=FloatVector(csf[:,0])

robjects.globalenv["CSF_cov"] = CSF_cov
robjects.globalenv["WM_cov"] = WM_cov

#laods csv file takes in missing data and replaces it with nan
print 'loading csv file ...'
mlm = numpy.genfromtxt('/Users/Home/Google Drive/Positivity Project/MLM test/Single_Beta_Analysis/setup_no_NA.csv', delimiter=",")

subject = FloatVector(mlm[:,0])
trial = FloatVector(mlm[:,2])
trial_type = FloatVector(mlm[:,3])
decision = FloatVector(mlm[:,4])
correct = FloatVector(mlm[:,5])
flex_search = FloatVector(mlm[:,6])
reg_ctl = FloatVector(mlm[:,7])
PA = FloatVector(mlm[:,8])
NA = FloatVector(mlm[:,9])
PA_NA = FloatVector(mlm[:,10])
SWLS = FloatVector(mlm[:,11])
COMT_rs4680_alleles = FloatVector(mlm[:,12])
COMT_rs4680_bin = FloatVector(mlm[:,13])
DARPP32_rs907094_alleles = FloatVector(mlm[:,14])
DARPP32_rs907094_bin = FloatVector(mlm[:,15])

robjects.globalenv["subject"] = subject
robjects.globalenv["trial"] = trial
robjects.globalenv["trial_type"] = trial_type
robjects.globalenv["decision"] = decision
robjects.globalenv["correct"] = correct
robjects.globalenv["flex_search"] = flex_search
robjects.globalenv["reg_ctl"] = reg_ctl
robjects.globalenv["PA"] = PA
robjects.globalenv["NA"] = NA
robjects.globalenv["PA_NA"] = PA_NA
robjects.globalenv["SWLS"] = SWLS
robjects.globalenv["COMT_rs4680_alleles"] = COMT_rs4680_alleles
robjects.globalenv["COMT_rs4680_bin"] = COMT_rs4680_bin
robjects.globalenv["DARPP32_rs907094_alleles"] = DARPP32_rs907094_alleles
robjects.globalenv["DARPP32_rs907094_bin"] = DARPP32_rs907094_bin


T0_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T1_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T2_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T3_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T4_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T5_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T6_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T7_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
T8_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)

beta0_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta1_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta2_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta3_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta4_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta5_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta6_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta7_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)
beta8_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], 1), dtype=numpy.float32)


r.library("nlme")
r.library("lmerTest")
from rpy2.robjects import Formula
noConverge = []


#got r
for x in xrange(int(sys.argv[2]), int(sys.argv[3])):
    print [x]
    for y in xrange(0, 108):
		for z in xrange(0, 90):
			if mask[x,y,z] > .5:
				#print [x, y, z]
				try:
					fMRI_data = FloatVector((s01[x,y,z,:]))
					robjects.globalenv["fMRI_data"] = fMRI_data
                    
					lmer_mri = r.lmer(Formula(" fMRI_data ~ seed_fmri + WM_cov + CSF_cov + (1|subject)"))
					lmer_mri = r.lmer(Formula(" fMRI_data ~ (seed_fmri + WM_cov + CSF_cov) * DARPP32_rs907094_bin + (1|subject)"))
                    
					lmer_mri = r.lmer(Formula(" fMRI_data ~ (fMRI_data + fMRI_seed) * DARPP32_rs907094_bin * correct + WM_cov + CSF_cov (1|subject)"))
					
                    betas = r.fixef(lmer_mri) #still the betas in lmer
					ses = r.sqrt(r.diag(r.vcov(lmer_mri))) #still works for lmer as well

					T0_img[x,y,z] = betas[0] / ses[0]
					T1_img[x,y,z] = betas[1] / ses[1]
					T2_img[x,y,z] = betas[2] / ses[2]
					T3_img[x,y,z] = betas[3] / ses[3]
					T4_img[x,y,z] = betas[4] / ses[4]
					T5_img[x,y,z] = betas[5] / ses[5]
					T6_img[x,y,z] = betas[6] / ses[6]
					T7_img[x,y,z] = betas[7] / ses[7]
					T8_img[x,y,z] = betas[8] / ses[8]

					beta0_img[x,y,z] = betas[0] 
					beta1_img[x,y,z] = betas[1] 
					beta2_img[x,y,z] = betas[2] 
					beta3_img[x,y,z] = betas[3]
					beta4_img[x,y,z] = betas[4] 
					beta5_img[x,y,z] = betas[5] 
					beta6_img[x,y,z] = betas[6]
					beta7_img[x,y,z] = betas[7]
					beta8_img[x,y,z] = betas[8]		


				except:
					print [x, y, z, ' did not converge']
					noConverge.append([x, y, z])

for i in xrange(0, 9):
    

    T_img="T{}_img".format(i)
    tempImg = nib.Nifti1Image(eval(T_img), pp01.get_affine())
    T_name='t{}_darpp32_{}.nii.gz'.format(i,sys.argv[1])
    print T_name
    tempImg.to_filename(T_name)
    
    beta_img="beta{}_img".format(i)
    tempImg = nib.Nifti1Image(eval(beta_img), pp01.get_affine())
    beta_name='beta{}_darpp32_{}.nii.gz'.format(i,sys.argv[1])
    tempImg.to_filename(beta_name)



