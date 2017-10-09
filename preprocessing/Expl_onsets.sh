#!/bin/sh 

# use code for retrieving Exploration data files for all participants
# create folders and subfolders for each participant
# create onset .txt files

#/Documents/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis
# run TD fMRI analysis for all participants

#cd "Documents/Positivity_fMRI_DATA/fMRI Behav/P1331 Jul23.2013/fMRI Tasks/Temporal Discounting/results"

#for item in $v;do
#	f=${item:15:4}; 
#	echo $f
#done

#for SubNo in 1282; do
SubNo=$1
echo $SubNo

path="/Volumes/HD5/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/preprocessing/${SubNo}"
mkdir $path
path2="/Volumes/HD5/Positivity_Project/Positivity_fMRI_DATA/fMRI Behav/Dec. 3. 2013/fMRI Tasks/Exploration/results"
path3="/Volumes/HD5/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/modelling/indirectActor/modelled_data"

# echo's participant number
logfilename=`ls "$path2/"*${SubNo}*.log`
resfilename=`ls "$path3/"*${SubNo}*.txt`
#k_filename="/Volumes/My_Passport/Positivity_Project/Positivity_fMRI_Analysis/TD_fMRI_Analysis/k_estimate.txt"

#for item in $v;do
#	echo $item;
#done

cd $path
# extract Stimulus onset from log file
t_scanner=`grep Start\ Slice "$logfilename" | awk '{print $7}'` # scanner starts
t_scanner=${t_scanner:0:5}
t_onset=`grep Stimulus\ Onset "$logfilename" | awk '{print $1}'` # stimulus onsets

## ONSET: create 3-col onset.txt file
for t in $t_onset ; do
  tsec=`echo "scale=3;($t-$t_scanner)/1000.0" | bc`
  echo "$tsec 3.0 1.0" >> onset.txt
done

## CHOICE: extract behavioural data to create .txt 3 column onset files 
choice=`grep $SubNo "$resfilename" | awk '{print $3}'` # choice
#choice=`awk '{print $9}' "$filename"`
#paste <(awk '{print $1, $2}' onset.txt ) <(awk '{print $1}' tmpchoice.txt) > choice_2.txt # doesn't compile :(
#choice=`awk '{print $1}' tmpchoice.txt` # choice
onset=(`awk '{print $1}' onset.txt`) # onset

in=0
for c in $choice; do
  echo "${onset[in]} 3.0 $c" >> choice.txt;
  ((in++))
done

## REWARD: CREATE 3-COLUMN centered_REWARD .TXT

tmprwd=`grep $SubNo "$resfilename" | awk '{print $5}'` # choice reward centered

in=0
for c in $tmprwd; do
  echo "${onset[in]} 3.0 $c" >> reward.txt;
  ((in++))
done

# EXPLOIT: CREATE 3-COLUMN exploration/exploitation .TXT FILE

exploit=`grep $SubNo "$resfilename" | awk '{print $6}'` # exploitation=1, exploration=0

in=0
for c in $exploit; do
  echo "${onset[in]} 3.0 $c" >> exploit.txt;
  ((in++))
done

#----