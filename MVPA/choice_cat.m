% Exploration MVPA
% May 28. 2016
% by Hans Melo


% Step 1: Set up structure
subj = init_subj('explore','1216'); % initinalize structure

% Step 2: Load Mask
subj = load_afni_mask(subj,'DLPFC_category-selective','/Volumes/HD1/Positivity_Project/masks/r_dlPFC_mask+tlrc'); % load mask

% for i=[1216 1256 1257 1258 1268 1269 1272 1273 1280 1281]
%   raw_filenames{i} = sprintf('mcdcdn_%i',i);
% end

% Step 3: Load fMRI data
filename = '/Volumes/HD1/Positivity_Project/Positivity_fMRI_Analysis/Exploration_fMRI_Analysis/Single_Beta_Analysis/mean/mcdcdn_1216+tlrc.BRIK';
subj = load_afni_pattern(subj,'fmri','DLPFC_category-selective',filename);

% Step 4: Load condition regressors

subj = init_object(subj,'regressors','conds_choice_cat');

% conditions loaded directly from csv file
% load('tutorial_regs');

% conditions loaded directly from csv file
% load('tutorial_regs');
choice1=ChoiceCat==1;
choice2=ChoiceCat==2;
choice3=ChoiceCat==3;
choice4=ChoiceCat==4;
choice5=ChoiceCat==5;
regs_ChoiceCat=[choice1 choice2 choice3 choice4 choice5]';
regs_ChoiceCat=double(regs_ChoiceCat);

subj = set_mat(subj,'regressors','conds_choice_cat',regs_ChoiceCat);
condnames_ChoiceCat = {'choice1','choice2','choice3','choice4','choice5'};
subj = set_objfield(subj,'regressors','conds_choice_cat','condnames',condnames_ChoiceCat);


% Step 5: Load run information
subj = init_object(subj,'selector','trials');
% load('tutorial_runs');
subj = set_mat(subj,'selector','trials',trial');

% Step 6: Z-Scoring

subj = zscore_runs(subj,'fmri','trials');

% Step 7: create cross-validation indices

subj = create_xvalid_indices(subj,'trials');

% Step 8: Feature Selection
conds_choice_cat
[subj] = feature_select(subj,'fmri_z','conds_choice_cat','trials_xval');

% Step 9: Classification

% set some basic arguments for a backprop classifier
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;

% now, run the classification multiple times, training and testing
% on different subsets of the data on each iteration

[subj results] = cross_validation(subj,'fmri_z','conds_choice_cat','trials_xval','fmri_z_thresh0.05',class_args);

% Step 10: Results

results

% Visualize Results

subj = load_afni_mask(subj,'wholebrain2','MNI152_T1_2mm_brain+tlrc');

% run in terminal AFNI:

3dresample \
  -master mcdcdn_1216+tlrc \
  -prefix MNI152_T1_2mm_brain_down \
  -inset MNI152_T1_2mm_brain+tlrc

subj = load_afni_pattern(subj,'anat','wholebrain2','MNI152_T1_2mm_brain_down+tlrc');

view_montage(subj,'pattern','anat','mask','wholebrain2')

view_montage(subj,'mask','wholebrain2','mask','DLPFC_category-selective')

view_montage(subj,'mask','wholebrain2','mask','fmri_z_thresh0.05_1');

plot_stability(subj,'fmri','trials')
