#!/bin/bash

# Prepares WM and CSF comfound vectors for connectivity analysis

#find Structural T1 Bravo for each subject

## Paths for Mac Pro
#conn_path='/Volumes/MacintoshHD3/Hans/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Connectivity'
#beta_path='/Volumes/MacintoshHD3/Hans/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_analysis/concatenated'

# Paths for Hans' Macbook Pro
conn_path='/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/Connectivity/mcdcdn_tmp'
beta_path='/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/mean'


for SubNo in $@
do
## Paths for Mac Pro
#struct_path='/Volumes/MacintoshHD3/Hans/Positivity_Project/Positivity_fMRI_DATA/T1s'
#struct_file=`ls "$struct_path/${SubNo}"_*`

# hans macbook path
fmri_path=`echo "/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_DATA/MOOD1MR/"*_$SubNo`

info_file=`ls "${fmri_path}"/Name_Info.txt`
struct_path=`grep BRAVO "$info_file" | awk '{print $1}'` 
struct_file=`ls "${fmri_path}/${struct_path}/"*_brain*`

## Segmentation using FAST on each subjects T1

cd $conn_path
mkdir $SubNo; cd $SubNo
mkdir WM_CSFconfounds
echo 'Segmentation for' $SubNo
fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g --nopve -o $conn_path/$SubNo/WM_CSFconfounds/${SubNo} $struct_file
# Output: 0=CSF; 1=GM; 2=WM 

## FLIRT transform segmentation to standard space: WM and CSF

echo 'Transform WM map for' $SubNo # In *.feat/reg/
flirt -in $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_seg_2.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain -out $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_WM_mask -omat $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_WM_mask.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear

echo 'Transform CSF map for' $SubNo
flirt -in $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_seg_0.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain -out $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_CSF_mask -omat $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_CSF_mask.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear

## Create binary masks out of standard space WM and CSF maps

echo 'Binarize WM map to create mask for' $SubNo
fslmaths $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_WM_mask -thr 0.9 -bin $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_WM_mask90_bin

echo 'Binarize CSF map to create mask for' $SubNo
fslmaths $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_CSF_mask -thr 0.9 -bin $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_CSF_mask90_bin


## Extracting mean ts for masked, standard WM and CSF

echo 'Extracting mean TS in WM and CSF for' $SubNo
# input subjects betas ts ; use mean_centered betas
fslmeants -i $beta_path/mcdcdn_${SubNo} -m $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_WM_mask90_bin -o $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_WM_confound_meants.txt

fslmeants -i $beta_path/mcdcdn_${SubNo} -m $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_standard_CSF_mask90_bin -o $conn_path/$SubNo/WM_CSFconfounds/${SubNo}_CSF_confound_meants.txt

cd $home
done

# need to concatenate ts for all subjects to include in model!!!
# cat file1.txt file2.txt file3.txt >>file.txt # could use cat *_WM* >> all_WM_confound_meants.txt if in same folder

#for SubNo in $subs; do echo $SubNo; cat "$conn_path/$SubNo/WM_CSFconfounds/${SubNo}_WM_confound_meants.txt" >> all_mcdcdn_WM_confound_meants.txt; cat "$conn_path/$SubNo/WM_CSFconfounds/${SubNo}_CSF_confound_meants.txt" >> all_mcdcdn_CSF_confound_meants.txt	;done



printf 'FINISHED!'


