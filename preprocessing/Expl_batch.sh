#!/bin/bash
# Exploration first level batch script
# for each subject generates and runs new .fsf file
 
#loads the fsl program
#export FSLDIR=/usr/local/packages/fsl
#.  ${FSLDIR}/etc/fslconf/fsl.sh
  
# set subjet variables
#v=(`ls /Volumes/My_Passport/Positivity_Project/Positivity_fMRI_DATA/MOOD1MR`)
# ${v[@]:1:37}
# ${v[@]:37:75}1395 1327 1426 1346

for item in $@;do
#	SubNo=${item:16:4}; 
	SubNo=$item
#	echo $SubNo
#done

#for item in 1341;do
#	SubNo=${item}; 
echo "Processing: ${SubNo}"

./Expl_onsets.sh ${SubNo}
output_path="/Volumes/HD5/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/preprocessing/${SubNo}"
expl_fmri_path=`echo "/Volumes/HD5/Positivity_Project/Positivity_fMRI_DATA/MOOD1MR/"*_${SubNo}`
#struct_file=`echo "${td_fmri_path}/Ex01367Se16/s016a1001_brain"`
fMRI_file=`echo "${expl_fmri_path}/EXP/sprl"`

info_file=`echo "${expl_fmri_path}/Name_Info.txt"`
struct_path=`grep BRAVO "$info_file" | awk '{print $1}'` 
struct_file=(`ls -1 "${expl_fmri_path}/${struct_path}/"s*.nii* | sed -e 's/\..*$//'`)
 
bet_filename="${struct_file[0]}_brain.nii.gz"

if [ -f $bet_filename ] # test if bet has been done
then
	echo "BET file found."
else
	echo "Processing BET for ${SubNo} ..."
	bet "${struct_file}" "$bet_filename"	
	echo "BET for ${SubNo} completed."
fi

 #makes the fsf files from the template fsf file
 for i in 'design_subj2.fsf'; do
  sed -e 's@output_path@'${output_path}'@g' \
   -e 's@struct_file@'${struct_file}'@g' \
   -e 's@fMRI_file@'${fMRI_file}'@g' <$i> ${output_path}/FEAT_${SubNo}.fsf
 done
 #runs the analysis using the newly created fsf file
feat ${output_path}/FEAT_${SubNo}.fsf
done