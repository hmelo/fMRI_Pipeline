function [choice_rwd_c, keys_c,b, trial, rt_c] = import_data(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [SUBNO,TRIAL,REWARD1,REWARD2,REWARD3,REWARD4,KEY,CHOICE_RWD,RT,JITTER]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [SUBNO,TRIAL,REWARD1,REWARD2,REWARD3,REWARD4,KEY,CHOICE_RWD,RT,JITTER]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [SubNo,trial,reward1,reward2,reward3,reward4,key,choice_rwd,rt,jitter]
%   = importfile('expl_daw_fmri_SubNo1216_[2013 6 11 15 11 55.765].res',2,
%   151);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/03/10 14:56:40

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
SubNo = dataArray{:, 1};
trial = dataArray{:, 2};
reward1 = dataArray{:, 3};
reward2 = dataArray{:, 4};
reward3 = dataArray{:, 5};
reward4 = dataArray{:, 6};
key = dataArray{:, 7}-27;
choice_rwd = dataArray{:, 8};
rt = dataArray{:, 9};
jitter = dataArray{:, 10};

% Get rid of no answer for modelling
b=find(dataArray{:, 7});
choice_rwd_c=choice_rwd(key~=-27);
keys_c = key(key~=-27);
rt_c = rt(key~=-27);
