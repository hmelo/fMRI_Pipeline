#!/bin/sh 

# use code for setting up file to input into hlm, contains all trial information for all participants

# /Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/

# use in terminal to get all subjects
# subjects=`for i in mcdcdn*;do echo ${i:7:4}; done`

# use to check how many


for SubNo in $@;do
	
	echo $SubNo

path="/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/preprocessing/${SubNo}"
#mkdir $path
#path2="/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_DATA/fMRI Behav/Dec. 3. 2013/fMRI Tasks/EFT/results"
#path2="/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_DATA/fMRI\ Behav/Dec.\ 3.\ 2013/fMRI\ Tasks/Attention/Results/"
sb_path="/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis"
modelling_file="/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/modelling/indirectActor/modelled_data/data_${SubNo}_PE.txt"
self_report_file="/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/Single_Beta_Analisis/google_setup.tsv"

# echo's participant number
#logfilename=`ls "$path2/"*${SubNo}*.log`
#resfilename=`ls "$path2/"*${SubNo}*.res`

cd $path

## Extract Behavioral and modelling data
#SubNo trial choice reward reward_cent exploit exp_par prob prob_cent action_value action_value_cent PE PE_cent
choice=(`grep $SubNo "${modelling_file}" | awk '{print $3}'`) # soon= -1, delayed = 1, no response=0
reward=(`grep $SubNo "${modelling_file}" | awk '{print $4}'`) 
reward_c=(`grep $SubNo "${modelling_file}" | awk '{print $5}'`) 
exploit=(`grep $SubNo "${modelling_file}" | awk '{print $6}'`)
expl_par=(`grep $SubNo "${modelling_file}" | awk '{print $7}'`) 
prob=(`grep $SubNo "${modelling_file}" | awk '{print $8}'`) 
prob_c=(`grep $SubNo "${modelling_file}" | awk '{print $9}'`) 
action_value=(`grep $SubNo "${modelling_file}" | awk '{print $10}'`) 
action_value_c=(`grep $SubNo "${modelling_file}" | awk '{print $11}'`)  
PE=(`grep $SubNo "${modelling_file}" | awk '{print $12}'`) 
PE_c=(`grep $SubNo "${modelling_file}" | awk '{print $13}'`) 


# Self-report
# load all available measures then compute mean and subtract mean for each subject - easier to excel it then save file.
# put all in one file and read from there!

PA=`grep $SubNo "$self_report_file" | awk '{print $3}'`
NA=`grep $SubNo "$self_report_file" | awk '{print $5}'`
PA_NA=`grep $SubNo "$self_report_file" | awk '{print $7}'`
SWLS=`grep $SubNo "$self_report_file" | awk '{print $9}'`

#echo $SubNo $PA $NA $PA_NA $SWLS
#echo ""
# Genotypes
COMT_rs4680_alleles=`grep $SubNo "$self_report_file" | awk '{print $11}'`;
COMT_rs4680_bin=`grep $SubNo "$self_report_file" | awk '{print $12}'`;
DARPP32_rs907094_alleles=`grep $SubNo "$self_report_file" | awk '{print $13}'`;
DARPP32_rs907094_bin=`grep $SubNo "$self_report_file" | awk '{print $14}'`; #echo $SubNo $DARPP32_rs907094_bin

# create/add to file
cd $sb_path
in=0
for i in ${choice[@]};do
	trial=`echo "$in+1"| bc -l`;	
	echo $SubNo $trial ${choice[in]} ${reward[in]} ${reward_c[in]} ${exploit[in]} ${expl_par[in]} ${prob[in]} ${prob_c[in]} ${action_value[in]} ${action_value_c[in]} ${PE[in]} ${PE_c[in]} $PA $NA $PA_NA $SWLS $COMT_rs4680_alleles $COMT_rs4680_bin $DARPP32_rs907094_alleles ${DARPP32_rs907094_bin}  >> setup.txt
	((in++))
done
done
echo 'finished!'
