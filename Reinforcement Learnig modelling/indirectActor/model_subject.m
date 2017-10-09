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
path='/Volumes/HD1/Positivity_Project/Positivity_fMRI_DATA/fMRI Behav/Dec. 3. 2013/fMRI Tasks/Exploration/results/';


%% Load data 

fileID = fopen('filelist.txt');
C = textscan(fileID,'%s',...
'Delimiter','\n');
fclose(fileID);
sum_data=[];all_subs_data=[];

for i=1:length(C{1})

file_name=sprintf('%s/%s',path,cell2mat(C{1}(i))); % get file name with full path
SubNo=cell2mat(C{1}(i)); SubNo=SubNo(20:23); % Get SubNo
[choice_rwd, choice,b,trial,rt] = import_data(file_name); % import relevant data
numChoice=4; % 4-arm bandir

%% Run models

% simple model - indirect actor 
[negLogLike, learningRate, exploration, td, exploitation, probs, mRec,action_value, chosen_prob] = indirectActorAuto(choice_rwd, choice, numChoice,1);

% indirect actor with decay
%[negLogLike, learningRate, exploration, td, exploitation] = indirectActor2Auto(choice_rwd_c, keys_c, numChoice);

% Kalman filter
%[negLogLike, learningRate, exploration, td, exploitation] = kalmanAuto(choice_rwd_c, keys_c, numChoice);


%% back to full set for fmri, with missing values as NaN's

data=NaN(150,4);
data(:,1)=repmat(str2double(SubNo), 150,1);
data(:,2)=trial;
data(b,3)=choice; % get all choices back indo right order with no answers as NaN's
data(b,4)=choice_rwd; % reward of chosen option
data(b,5)=choice_rwd-mean(choice_rwd); % mean centered reward
data(b,6)=exploitation; % exploit=1, explore=0
data(:,7)=repmat(exploration, 150,1); % exploration parameter for each subject
data(b,8)=chosen_prob; % probability of each chosen option
data(b,9)=chosen_prob-mean(chosen_prob); % mean centered probability of choice
data(b,10)=action_value; % action value of chosen option
data(b,11)=action_value-mean(action_value); % mean centered action calues
data(b,12)=td; %prediction error
data(b,13)=td-mean(td); %mean centered prediction error
data(b,14)=rt; %reaction time 
data(b,15:18)=probs; %probabilities
data(b,19:22)=mRec; %action values


% cd 'modelled_data'
headers='SubNo  trial   choice  reward  reward_cent exploit exp_par prob    prob_cent   action_value    action_value_cent   PE  PE_cent rt prob1 prob2 prob3 prob4 av1 av2 av3 av4';
% dlmwrite(sprintf('data_%s_PE_rt.txt',SubNo),headers,'-append','delimiter','')
% dlmwrite(sprintf('data_%s_PE_rt.txt',SubNo),data(1:149,:),'-append','delimiter','\t')
% cd ..

all_subs_data=[all_subs_data;data(1:149,:)];

explore_rt = mean(rt(data(b,6)==0)); % rt for explore trials
exploit_rt = mean(rt(data(b,6)==1)); % rt for exploit trials
per_explore =  sum(exploitation==0)/length(exploitation); % percent of explore trials
sum_data=[sum_data; str2double(SubNo) explore_rt exploit_rt per_explore exploration learningRate];

end

dlmwrite('all_subs_data.txt',headers,'-append','delimiter','')
dlmwrite('all_subs_data.txt',all_subs_data,'-append','delimiter','\t')
