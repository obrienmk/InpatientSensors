function CVAPareticSide = get_side(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   CVAPARETICSIDE = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   CVAPARETICSIDE = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   CVAPareticSide = importfile('CVA_Paretic_Side.csv', 1, 55);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/03/24 12:51:46

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%s%[^\n\r]';

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

%% Create output variable
CVAPareticSide = table(dataArray{1:end-1}, 'VariableNames', {'ID','Side'});

