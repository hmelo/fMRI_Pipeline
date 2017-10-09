####################################################################################################
# Single Beta fmri analysis 
# Cunningham Lab Oct. 2014
# put together by Hans Melo

# This script does single beta analysis of fMRI data, using R.
# it also controls for mean activation for every subject (fmri+fmri_mean)
# This version uses mpi4py for parallel processing
# must have mpi4py installed, as well as all other required libraries

# to run from terminal type:
# mpiexec -n 4 python MPI_from_currentVersion_tmp_Hans.py
# where 4 is the number of threads desired (normally 2x number of processors)
####################################################################################################


#################################
# Set up model and files below! #
#################################
# specify R model for lmer
R_model = "fMRI_data ~  (fMRI_seed * COMT_rs4680_bin) + WM_COV + CSF_COV + (1|subject)"

 
# specify required files to load
maskFile = '/Users/Home/Google Drive/Positivity Project/MLM test/hlmflexfmri/mask.nii.gz' # mask i.e. whole brain 
centeredFile = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/mcdcdn_functional_merged.nii' # mean centered fmri betas
meanFile = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/mean_functional_merged.nii' # mean images
FC_seed = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/Decision_SV_gene_k/Clust_02_t5_th3_m4_m4_m12.1D'
WM_COV_file = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/all_mcdcdn_WM_confound_meants.txt'
CSF_COV_file = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/all_mcdcdn_CSF_confound_meants.txt'

csvFile = '/Users/Home/Google Drive/Positivity Project/TD_Single_Beta_Analysis/setup.csv' #csv file, it's ok to have a header line

####################################################################################################

# regarding memory https://github.com/nipy/nibabel/pull/211
# http://nipy.bic.berkeley.edu/nightly/nibabel/doc/devel/image_design.html


#load required libraries
import numpy
import nibabel as nib

from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import Formula
r = robjects.r
r.library('lme4')
stats = importr('stats')
base = importr('base')

# set up mpi4py for parallel processing
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank() # process number
nproc = comm.Get_size() # number of processors required


# load files
print("Loading mask image...")
maskImage = nib.load(maskFile)
mask = maskImage.get_data()

print("Loading centered data...")
pp01 = nib.load(centeredFile)
s01 = pp01.get_data()

print("Loading mean images...")
meanImageFile = nib.load(meanFile)
meanImages = meanImageFile.get_data()

print("loading seed file")
seed = numpy.genfromtxt(FC_seed, delimiter=',', skip_header=0) # it's ok to have a header line
WM_COV = numpy.genfromtxt(WM_COV_file, delimiter=',', skip_header=0) # it's ok to have a header line
CSF_COV = numpy.genfromtxt(CSF_COV_file, delimiter=',', skip_header=0) # it's ok to have a header line

print("loading MLM csv file")
mlm = numpy.genfromtxt(csvFile, delimiter=',')#, skip_header=0) # it's ok to have a header line

# need to specify variable names below
# could improve to read directly from csv header
fMRI_seed = FloatVector(seed)
WM_COV = FloatVector(WM_COV)
CSF_COV = FloatVector(CSF_COV)
subject = FloatVector(mlm[:,0])
trial = FloatVector(mlm[:,1])
decision = FloatVector(mlm[:,3])
SV = FloatVector(mlm[:,4])
rwd_diff = FloatVector(mlm[:,5])
time_diff = FloatVector(mlm[:,6])
k = FloatVector(mlm[:,7])
PA = FloatVector(mlm[:,8])
NA = FloatVector(mlm[:,9])
PA_NA = FloatVector(mlm[:,10])
SWLS = FloatVector(mlm[:,11])
COMT_rs4680_alleles = FloatVector(mlm[:,12])
COMT_rs4680_bin = FloatVector(mlm[:,13])
DARPP32_rs907094_alleles = FloatVector(mlm[:,14])
DARPP32_rs907094_bin = FloatVector(mlm[:,15])

robjects.globalenv["fMRI_seed"] = fMRI_seed
robjects.globalenv["WM_COV"] = WM_COV
robjects.globalenv["CSF_COV"] = CSF_COV
robjects.globalenv["subject"] = subject
robjects.globalenv["trial"] = trial
robjects.globalenv["decision"] = decision
robjects.globalenv["SV"] = SV
robjects.globalenv["rwd_diff"] = rwd_diff
robjects.globalenv["time_diff"] = time_diff
robjects.globalenv["k"] = k
robjects.globalenv["PA"] = PA
robjects.globalenv["NA"] = NA
robjects.globalenv["PA_NA"] = PA_NA
robjects.globalenv["SWLS"] = SWLS
robjects.globalenv["COMT_rs4680_alleles"] = COMT_rs4680_alleles
robjects.globalenv["COMT_rs4680_bin"] = COMT_rs4680_bin
robjects.globalenv["DARPP32_rs907094_alleles"] = DARPP32_rs907094_alleles
robjects.globalenv["DARPP32_rs907094_bin"] = DARPP32_rs907094_bin

###### mpi4py #######
comm.Barrier() # wait for all processes to finish before proceding
## SPLIT 90 INTO PIECES FOR ALL PROCESSORS
# for the loop
npts = 90    # number of loops to split into processors
nmin = npts/nproc # number of elements per processor 
nextra = npts%nproc # number of extra elements
k=0 # displacements

nsendcounts=numpy.empty([nproc,1]) #length of each process's portion of the original vector
ndispls=numpy.empty([nproc,1]) # displacement of each portion

for i in xrange(0, nproc):
    if i<nextra:
        nsendcounts[i] = nmin+1 # add 1 if extra
    else:
        nsendcounts[i] = nmin # don't add if not extra
    
    ndispls[i] = k; # update displacement
    k = k+nsendcounts[i]; # update interval
    
# assign nsendcounts for each process
local_start = ndispls[rank] #start loop at
local_end = ndispls[rank]+nsendcounts[rank] # end loop at     
print "rank", rank, "nsendcounts",nsendcounts[rank], "ndispls", ndispls[rank]

########

# runs sample R model to extract number of betas

#numBetas =example_funct(s01,meanImages, ...);

x = 33; y = 62; z = 23

fMRI_data = FloatVector((s01[x,y,z,:])) # takes sample points
robjects.globalenv["fMRI_data"] = fMRI_data
trial_num=s01.shape[3]/meanImages.shape[3] # GET NUMBER OF TRIALS (total number of trials devided by number of participants, assuming equal number of trials for all)
i = 0
meanfMRI = numpy.zeros(len(mlm), dtype=numpy.float32)
for subjCount in xrange(0, meanImages.shape[3]):
#    print subjCount # printing subject ok, no problem!
    for trialCount in xrange(0, trial_num): # improved to automatically use number of trials
	    meanfMRI[i] = meanImages[x,y,z,subjCount]
	    i = i+1
    fMRI_meanIs = FloatVector(meanfMRI)
    robjects.globalenv["fMRI_meanIs"] = fMRI_meanIs
    lme_mri = r.lmer(Formula(R_model))
    betas = r.fixef(lme_mri) 

numBetas = len(betas) 

#Creates arrays for the T, beta, and p scores
T_img_local = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)
beta_img_local = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)
p_img_local = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)


r.library("nlme")
r.library("lmerTest")
from rpy2.robjects import Formula
noConverge = []


#LOCAL PROCESSING
for x in xrange(local_start, local_end):
    print [x]
    for y in xrange(0, 108):
		for z in xrange(0, 90):
			if mask[x,y,z] > .5:
				#print [x, y, z]
				try:
					fMRI_data = FloatVector((s01[x,y,z,:]))
					robjects.globalenv["fMRI_data"] = fMRI_data
					lmer_mri = r.lmer(Formula(R_model))
					#modsum = base.summary(fit)
     					betas = r.fixef(lmer_mri)
					ses = r.sqrt(r.diag(r.vcov(lmer_mri)))
					#b = modsum.rx2('coefficients')[0:numBetas]
					#t = modsum.rx2('coefficients')[2*numBetas:3*numBetas]
					#p = modsum.rx2('coefficients')[3*numBetas:4*numBetas] 
					for BetaCount in xrange(0,numBetas):
					    T_img_local[x,y,z,BetaCount] = betas[BetaCount]/ses[BetaCount]
					    beta_img_local[x,y,z,BetaCount] = betas[BetaCount]
					    #p_img_local[x,y,z,BetaCount] = (p[BetaCount]-1)*-1
				except:
					print [x, y, z, ' did not converge']
					noConverge.append([x, y, z])

comm.Barrier() # wait for all processes to finish before proceding

# pre-allocate images for merged results
T_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)
beta_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)
p_img = numpy.zeros((s01.shape[0], s01.shape[1], s01.shape[2], numBetas), dtype=numpy.float32)

# mpi4py # Merge results from all processes

comm.Reduce(T_img_local, T_img, op = MPI.SUM)
comm.Reduce(beta_img_local, beta_img, op = MPI.SUM)
comm.Reduce(p_img_local, p_img, op = MPI.SUM)
## write stats images
# Images have the format such that each slice represents the T/b/p val for the corresponding regression parameter
# in this example, t1 = fMRI_data, t2 = fMRI_meanIs, etc. (try code in R if you're unsure what the order things come out in)
if rank == 0: # only one processor needs to do this
    print "Writing stats images ..."
    tempImg = nib.Nifti1Image(T_img, pp01.get_affine())
    tempImg.to_filename('T_img.nii.gz')
    tempImg = nib.Nifti1Image(beta_img, pp01.get_affine())
    tempImg.to_filename('beta_img.nii.gz')
    tempImg = nib.Nifti1Image(p_img, pp01.get_affine())
    tempImg.to_filename('p_img.nii.gz')
    print numBetas
#    print modsum # print stats table



