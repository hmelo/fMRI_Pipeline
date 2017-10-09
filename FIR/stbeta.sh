#!/bin/bash

# use for preparing data for single beta analysis

# 1. Creates registered normalized nii from preprocessed fsl nii's
# 2. Single-trial deconvolve, bring in trial onset info
# 3. Concatenate AFNI sub-bricks, create write individual concatenated deconvolved normalized (cdn_) .nii files
# 4. Despiking each subject deconvolved data
# 5. Mean centering
# 6. Create massive concatenated .nii
# 7. Create cocatenated mean

# to obtain all subject numbers in the folder
# subs=`for subject in cdn*; do echo ${subject:4:4};done`


######
#home=/Volumes/MacintoshHD3/YaleClinicalData/fMRI/single_trials_analysis
#cd $home

## Start AFNI
#echo 'export PATH=$PATH:$HOME/abin' >> .profile
#echo 'export DYLD_FALLBACK_LIBRARY_PATH=$HOME/abin' >> .profile
#. .profile

##################################################################################
## full path to subject_list
#subject_list=${home}/subjects.txt
## location of onset files
#onset_data=${home}/st_onsets.txt

##################################################################################

#d_dir=/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/deconvolved
#c_dir=/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/concatenated
#fsl_dir=/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/FSL_Analysis/Subj_Value
#m_dir=/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/mean


fsl_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR
m_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR/Single_Beta_Analysis/mean
d_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR/Single_Beta_Analysis/deconvolved
c_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR/Single_Beta_Analysis/concatenated

#fsl_dirs=/scratch/s/scslab/hansmelo/Positivity_Project/Positivity_fMRI_Analysis/Flexibility_fMRI_Analysis/Flexibility_fMRI_Analysis/FSL_Analysis
sb_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR/Single_Beta_Analisis
p_dir=/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/FIR/Single_Beta_Analisis/preprocessed

mkdir $p_dir; mkdir ${d_dir};mkdir ${c_dir}; mkdir ${m_dir}

for item in $@;do
	subject=${item:0:4} 

echo "Processing: ${subject}"

cd $fsl_dir
echo ''
echo 'Constructing registered preprocessed .nii for' $subject
echo ''

## 1. Creates registered normalized nii from preprocessed fsl nii's
featregapply ./$subject.feat -l filtered_func_data
##
mv ./${subject}.feat/reg_standard/filtered_func_data.nii.gz ./${subject}.feat/reg_standard/${subject}_filtered_func_data.nii.gz

#done
## Single-trial Deconvolve 
#for subject in $@
#do
cd ${d_dir}
mkdir tmp;
mkdir tmp/${subject}
mkdir tmp/${subject}/st_output

## 2. Single-trial deconvolve
# NB: The input file has not been cleaned of WM and CSF confounds!
onset=(`awk '{print $1}' ${fsl_dir}/onsets/${subject}/onset.txt`)
onsets=${onset[@]:0:149} # take only the first 150 trials
echo $onsets> tmp/${subject}/onsets_only.txt
#awk '{print $1}' ${fsl_dir}/${subject}/onset.txt > tmp/${subject}/onsets_only.txt # make onsets only file

echo ''
echo 'Processing single-trial Deconvolve for' ${subject} 
echo ''
3dDeconvolve -input ${fsl_dir}/${subject}.feat/reg_standard/${subject}_filtered_func_data.nii.gz \
-num_stimts 1 \
-stim_times_IM 1 tmp/${subject}/onsets_only.txt 'SPMG1' -stim_label 1 onsets \
-bucket ${d_dir}/tmp/${subject}/st_output/st_stats -nofullf_atall

3dAFNItoNIFTI -prefix ${d_dir}/${subject}_st_statsout ${d_dir}/tmp/${subject}/st_output/st_stats+tlrc.BRIK # makes .nii from deconvolved file

## 3. Concatenate AFNI sub-bricks and extract mean ts

cd ${c_dir}

echo 'Concatenating .nii for' ${subject} 

3dTcat -prefix ${c_dir}/cdn_${subject}.nii ${d_dir}/${subject}_st_statsout.nii[0..$] -tr 2  # 

## 4. Despiking
echo 'Despiking' ${subject} 
3dDespike -ignore 0 -nomask -prefix despiked_cdn_${subject}.nii.gz cdn_${subject}.nii
## 5. Mean centering
echo 'Mean centering' ${subject}

Fslmaths ${c_dir}/despiked_cdn_${subject}.nii.gz -Tmean ${m_dir}/mean_despiked_${subject}.nii.gz # CHECK THAT IT DOES IT RIGHT!
Fslmaths ${c_dir}/despiked_cdn_${subject}.nii -sub ${m_dir}/mean_despiked_${subject}.nii.gz ${m_dir}/mcdcdn_${subject}.nii.gz


cd $home
done

printf 'Finished!' 

### Run separately in Terminal after above processes are done.

## 6. Merge all .nii files into massive nii
# fslmerge -t mcdcdn_functional_merged.nii.gz mcdcdn*

## 7. Merge all means .nii files into massive nii to use in model
# fslmerge -t mean_functional_merged.nii.gz mean*
