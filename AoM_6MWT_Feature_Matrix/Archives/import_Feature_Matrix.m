function MWT6AoMFeatureMatrix360 = import_Feature_Matrix(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MWT6AOMFEATUREMATRIX360 = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   MWT6AOMFEATUREMATRIX360 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads
%   data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   MWT6AoMFeatureMatrix360 = importfile('MWT6_AoM_Feature_Matrix360.csv', 2, 141);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/07/08 16:21:33

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]);
rawStringColumns = string(raw(:, 2));


%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% Create output variable
MWT6AoMFeatureMatrix360 = table;
MWT6AoMFeatureMatrix360.Sub_ID = cell2mat(rawNumericColumns(:, 1));
MWT6AoMFeatureMatrix360.Sub_Type = categorical(rawStringColumns(:, 1));
MWT6AoMFeatureMatrix360.FIM = cell2mat(rawNumericColumns(:, 2));
MWT6AoMFeatureMatrix360.BBS = cell2mat(rawNumericColumns(:, 3));
MWT6AoMFeatureMatrix360.MWT10_SSV = cell2mat(rawNumericColumns(:, 4));
MWT6AoMFeatureMatrix360.MWT10_FV = cell2mat(rawNumericColumns(:, 5));
MWT6AoMFeatureMatrix360.MWT6 = cell2mat(rawNumericColumns(:, 6));
MWT6AoMFeatureMatrix360.TUG = cell2mat(rawNumericColumns(:, 7));
MWT6AoMFeatureMatrix360.AoM_Pel_tilt = cell2mat(rawNumericColumns(:, 8));
MWT6AoMFeatureMatrix360.AoM_Pel_ro = cell2mat(rawNumericColumns(:, 9));
MWT6AoMFeatureMatrix360.AoM_Pel_oblq = cell2mat(rawNumericColumns(:, 10));
MWT6AoMFeatureMatrix360.AoM_Ankle_US_x = cell2mat(rawNumericColumns(:, 11));
MWT6AoMFeatureMatrix360.AoM_Ankle_US_y = cell2mat(rawNumericColumns(:, 12));
MWT6AoMFeatureMatrix360.AoM_Ankle_US_z = cell2mat(rawNumericColumns(:, 13));
MWT6AoMFeatureMatrix360.AoM_Ankle_AS_x = cell2mat(rawNumericColumns(:, 14));
MWT6AoMFeatureMatrix360.AoM_Ankle_AS_y = cell2mat(rawNumericColumns(:, 15));
MWT6AoMFeatureMatrix360.AoM_Ankle_AS_z = cell2mat(rawNumericColumns(:, 16));
MWT6AoMFeatureMatrix360.AoM_Pel_Norm = cell2mat(rawNumericColumns(:, 17));
MWT6AoMFeatureMatrix360.AoM_Ankle_US_Norm = cell2mat(rawNumericColumns(:, 18));
MWT6AoMFeatureMatrix360.AoM_Ankle_AS_Norm = cell2mat(rawNumericColumns(:, 19));
MWT6AoMFeatureMatrix360.Steps = cell2mat(rawNumericColumns(:, 20));

