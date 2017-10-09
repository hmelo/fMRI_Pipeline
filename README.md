# fMRI_Pipeline

fMRI Analysis

Step 1: Preprocessing and first level (can also do higher level by setting up the fsl design)

Package: FSL via shell scripts

1.	Must set up model ONCE for one participant to get .fsf file. (you’ll need onset files)
2.	Then change elements in .fsf file to make it generic so that a batch script can be used. See example design ‘design_subj2.fsf' file 
3.	Use batch script to run for all subjects. Example script: ‘Expl_batch.sh’
4.	The above will generate ‘.feat’ folders for all subjects

Step 2: Single Beta Analysis: this is where things get interesting

1.	you’ll need to make a set up .csv file with all the information for all subjects and all trials. See example ‘setup_headers.csv’
2.	you’ll need to create a single .nii file that contains all fmri data for all subjects and all trials. To do this you’ll use a shell script ‘stbeta.sh’
3.	stbeta.sh takes all the preprocessed fmri data from Step 1 and runs a series of subroutines to get it ready for merging. 
4.	Merging: I recommend that you run the actual merging (steps 6 and 7 in the stbeta.sh file) separately in terminal. Of course you can just run it directly on the script, but in practice you’ll need to check every subject’s data and make sure it’s good. I usually get a couple of bad data which I choose not to include in the merged file, which is why I usually end up running this final steps separately in terminal. The resulting file will be MASSIVE because it contains all fmri trials for all participants (for 70 subjects I got ~30GB)

Step 3: Run actual single trial beta model: this is the part that involves the python script.

1.	Install mpi4py : I HIGHLY recommend using mpi4py. It can be a bit tedious to install it, but once it’s working. It allows parallel processing and can save you a lot of time, days maybe even weeks. 
2.	Setup model in python script. See ‘MLM_decision.py’ for an example. All you really need to do here is set up your model at the start of the file, and specify your paths i.e. where your data is. *I spent a lot of time figuring on this script, and it’s quite efficient in terms of processing and practicality.
3.	Wait for your model to run. If you set up your files correctly and followed all steps above this should run smoothly.
4.	Relish on your results.
