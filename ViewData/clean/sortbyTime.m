function [ DataSorted ] = sortbyTime( Data )
% *************************************************************************
% clean_sortbyTime.m
%
% This sorts all fields of a structure by time.
%
%     Input:   Data  = structure with field time (Data.time) and other
%                      unknown fields
%     Output:  DataSorted   = structure with all fields sorted by time
%
% Megan O'Brien, 2018
% *************************************************************************

% Initialize
DataSorted = Data;

% Get all fieldnames in structure
fnames = fieldnames(Data);

if ismember('time',fnames)
    % Get indices of sorted time
    [DataSorted.time, itime] = sort(Data.time);
    
    % Loop through all fields and arrange them by sorted time index
    for i = 1:length(fnames)
        
        field = fnames{i};
        fieldData = Data.(field);
        fieldDataSorted = fieldData(itime,:);
        
        DataSorted.(field) = fieldDataSorted;
    end
    
else
    error(['Fieldname ''time'' must be in data structure ' structname]);
end
    
end

