function [ DataNoOutliers ] = removeOutliers( Data, StdGain )
% *************************************************************************
% clean_removeOutliers.m
%
% This function removes outlier points
%
%     Input:   Data  = structure with fields lat, lng, time
%              StdGain = multiplier on standard deviation. Values outside
%              of this range (+/-) from the mean are considered outliers
%              and will be removed.
%     Output:  DataNoOutliers   = structure with outliers removed
%
% Megan O'Brien, 2018
% *************************************************************************

% Initialize
DataNoOutliers = Data;

% Get all fieldnames in structure
fnames = fieldnames(Data);
    
% Get location of time field and other fields in Data structure
ifields_time = strfind(fnames,'time');
ifields_other = find(cellfun('isempty', ifields_time));


% ------------------------- OUTLIER CRITERIA --------------------------
ioutlier = find(Data.lat > mean(Data.lat) + StdGain*std(Data.lat) | ...
    Data.lat < mean(Data.lat) - StdGain*std(Data.lat) | ...
    Data.lng > mean(Data.lng) + StdGain*std(Data.lng) | ...
    Data.lng < mean(Data.lng) - StdGain*std(Data.lng) );
% ---------------------------------------------------------------------


% Loop through all fields and remove indices where outlier was detected
for i = 1:length(fnames)
    field = fnames{i};
    DataNoOutliers.(field)(ioutlier,:) = [];
end

fprintf('-->%i OUTLIER POINTS REMOVED \n',length(Data.time) - length(DataNoOutliers.time));


end

