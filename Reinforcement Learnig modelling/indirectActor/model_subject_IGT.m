% Modelling of Exploration data
% uses softmax RL modelling for n-arm bandits
% Hans Melo, March 2015
%
% Takes a list of behavioral data files
% extracts data and runs, modelling indirectActorAuto for each separately
%
% calculates all relevant parameters including:
% NLL, learning rate, exploration parameter, prediction error,
% probabilities, action values.
%
% Outputs a file with 
% [SubNo trial choice reward reward_centered exploitation(1=exploit, 0=explore) exploration_parameter chosen_prob chosen_prob_centered action_value action_value_centered prediction_error] prediction_error_centered

clear all
clc

% Specify path to behavioral data
path='/Volumes/HD1/Positivity_Project/Positivity_fMRI_DATA/fMRI Behav/Dec. 3. 2013/fMRI Tasks/Exploration/EEG_results/';


%% Load data 

num_trials=length(choice);
data=NaN(num_trials,4);

subs=unique(subject);
sessions=unique(session);

for j=1:length(subs)

sub=subs(j);
sub_trials=find(subject==sub);
sub_task_trials=find(taskcontrol(sub_trials)==1);

rename_choice=choice;
rename_choice(rename_choice==3)=1;
rename_choice(rename_choice==4)=2;
rename_choice(rename_choice==8)=3;
rename_choice(rename_choice==9)=4;

for current_session=1:4,
    
sub_choice=[];sub_choice_rwd=[];good_trials=[];
for current_trial=sub_trials(1):sub_trials(end)
if taskcontrol(current_trial)==1 && strcmp(session(current_trial),sessions(current_session)),
    sub_choice=[sub_choice rename_choice(current_trial)];
    sub_choice_rwd=[sub_choice_rwd rwd(current_trial)];
    good_trials=[good_trials current_trial];
end
end

numChoice=4; % 4-arm bandir

%% Run models

% simple model - indirect actor 
[negLogLike, learningRate, exploration, td, exploitation, probs, mRec,action_value, chosen_prob] = indirectActorAuto(sub_choice_rwd, sub_choice, numChoice,1);

% indirect actor with decay
%[negLogLike, learningRate, exploration, td, exploitation] = indirectActor2Auto(choice_rwd_c, keys_c, numChoice);

% Kalman filter
%[negLogLike, learningRate, exploration, td, exploitation] = kalmanAuto(choice_rwd_c, keys_c, numChoice);


%% back to full set for fmri, with missing values as NaN's

% record modelling results

for i=1:length(good_trials)
    real_trial=good_trials(i);
    data(real_trial,1)=good_trials(i);
    data(real_trial,2)=exploitation(i);
    data(real_trial,3)=td(i)-mean(td);
    data(real_trial,4)=chosen_prob(i);
end

%data(800:900,:)

% cd 'modelled_data'
%headers='SubNo  trial   choice  reward  reward_cent exploit exp_par prob    prob_cent   action_value    action_value_cent   PE  PE_cent rt prob1 prob2 prob3 prob4 av1 av2 av3 av4';
% dlmwrite(sprintf('data_%s_PE_rt.txt',SubNo),headers,'-append','delimiter','')
% dlmwrite(sprintf('data_%s_PE_rt.txt',SubNo),data(1:149,:),'-append','delimiter','\t')
% cd ..

% all_subs_data=[all_subs_data;data(1:75,:)];

% explore_rt = mean(rt(data(b,6)==0)); % rt for explore trials
% exploit_rt = mean(rt(data(b,6)==1)); % rt for exploit trials
% per_explore =  sum(exploitation==0)/length(exploitation); % percent of explore trials
% sum_data=[sum_data; str2double(SubNo) explore_rt exploit_rt per_explore exploration learningRate];
end % session
end

dlmwrite('413_subs_EEG_data.txt',headers,'-append','delimiter','')
dlmwrite('413_subs_EEG_data.txt',all_subs_data,'-append','delimiter','\t')
