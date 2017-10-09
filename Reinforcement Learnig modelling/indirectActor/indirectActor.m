function [negLogLike, exploration] = indirectActor(learningRate, reward, choice, numChoice, nodisplay)
%indirectActor  Negative log likelihood of the indirect actor model on the
%   basis of behavioral data from a group of participants. It assume that all
%   parameters are the same across pariticipants except for explaration
%   parameter of Softmax.
%
%   NLL = indirectActor(parameter, reward, choice, numChoice, nodisplay)
%   returns negative log likelihood NLL given behavioral data and parameters.
%
%   Behavioral data: reward, choice, numChoice
%       reward : a matrix of received rewards of participants
%       choice : a matrix of choices of participants
%           In the matrix, row represents trial and column is participant.
%       numChoice : the number of choices in the experiment
%
%   learningRate : learning rate of the indirect actor
%
%   optional parameter : nodisplay
%       nodisplay = 0, parameters being searched is displayed. (default)
%       nodisplay = 1, it is not displayed.
%       
%   [NLL, E] = indirectActor(...)
%   returns an array of exploration parameters E.
%       each value in the array is the optimum value of each participant.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

if nargin == 4
    nodisplay = 0; % set to default
end

[trials numOfData] = size(choice);
% get the number of trials and participants

initValues = [0.2 0.5 0.8];
initResult = zeros(1, 3);
% use 3 different exploration parameters to find a starting point

negLogLike = 0;
for i = 1:numOfData
    for j = 1:3
        initResult(j) = indirectActorIndv(initValues(j), reward, choice, numChoice, learningRate);
    end
    [minV minI]                     = min(initResult);
    % initial exploration parameter is set to the one with minimum
    % negative log likelihood. 
    
    [exploration(i) indvResult]     = fminsearch(@indirectActorIndv, initValues(minI), [], reward(:, i), choice(:, i), numChoice, learningRate);
    % get the optimum value of exploration parameter 
    % and its negative log likelihood
    
    negLogLike                      = negLogLike + indvResult;
end

if (nodisplay == 0)
    disp(['Learning rate = ' num2str(learningRate) ', NLL = ' num2str(negLogLike)]);
end
