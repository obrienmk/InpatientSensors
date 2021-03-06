function FeatureMatrixAoMAdmission6MWT = import_AoM(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   FEATUREMATRIXAOMADMISSION6MWT = IMPORTFILE(FILENAME) Reads data from
%   text file FILENAME for the default selection.
%
%   FEATUREMATRIXAOMADMISSION6MWT = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   FeatureMatrixAoMAdmission6MWT = importfile('Feature_Matrix_AoM_Admission_6MWT.csv', 2, 56);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/06/04 22:45:35

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: text (%q)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: double (%f)
%	column46: double (%f)
%   column47: double (%f)
%	column48: double (%f)
%   column49: double (%f)
%	column50: double (%f)
%   column51: double (%f)
%	column52: double (%f)
%   column53: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%q%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
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
FeatureMatrixAoMAdmission6MWT = table(dataArray{1:end-1}, 'VariableNames', {'ID','AS','Time','Gyr_Pel_tilt_BW','Gyr_Pel_tilt_FW','Gyr_Pel_Oblq_US','Gyr_Pel_Oblq_AS','Gyr_Pel_Ro_US','Gyr_Pel_Ro_AS','Gyr_Pel_norm','Gyr_Ankle_US_x','Gyr_Ankle_AS_x','Gyr_Ankle_US_y','Gyr_Ankle_AS_y','Gyr_Ankle_US_z','Gyr_Ankle_AS_z','Gyr_Ankle_US_norm','Gyr_Ankle_AS_norm','Acc_Pel_x_US','Acc_Pel_x_AS','Acc_Pel_y_P','Acc_Pel_y_N','Acc_Pel_z_P','Acc_Pel_z_N','Acc_Pel_norm','Acc_Ankle_US_x','Acc_Ankle_AS_x','Acc_Ankle_US_y','Acc_Ankle_AS_y','Acc_Ankle_US_z','Acc_Ankle_AS_z','Acc_Ankle_US_norm','Acc_Ankle_AS_norm','Gyr_SI_Pel_Tilt','Gyr_SI_Pel_Oblq','Gyr_SI_Pel_Ro','Acc_SI_Pel_x','Acc_SI_Pel_y','Acc_SI_Pel_z','Gyr_SI_Ankle_x','Gyr_SI_Ankle_y','Gyr_SI_Ankle_z','Gyr_SI_Ankle_norm','Acc_SI_Ankle_x','Acc_SI_Ankle_y','Acc_SI_Ankle_z','Acc_SI_Ankle_norm','FIM','BBS','MWT10_SSV','MWT10_FV','MWT6','TUG'});

