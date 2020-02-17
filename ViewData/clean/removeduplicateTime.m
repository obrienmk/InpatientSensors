function [ DataNoDuplicates ] = removeduplicateTime( Data )
% *************************************************************************
% clean_removeduplicateTime.m
%
% This sorts removes any duplicate times from data.
%
%     Input:   Data  = structure with field time (Data.time) and other
%                      unknown fields
%     Output:  DataNoDuplicates   = structure without duplicate times
%
% Megan O'Brien, 2018
% *************************************************************************

% Initialize
DataNoDuplicates = Data;

% Get all fieldnames in structure
fnames = fieldnames(DataNoDuplicates);

if ismember('time',fnames)    
    % Get indices of duplicated time
    irep=find(diff(DataNoDuplicates.time)==0);  % index of 1st replicate
    
    % Loop through fields and average values where time is duplicated
    elapsedTime = 0;
    while ~isempty(irep)
        tic           
        for ifield = 1:length(fnames)
            field = fnames{ifield};
            DataNoDuplicates.(field)(irep,:)=(DataNoDuplicates.(field)(irep,:)+DataNoDuplicates.(field)(irep+1,:))./2;
            DataNoDuplicates.(field)(irep+1,:)=[];  % Remove duplicated rows
        end
        elapsedTime = toc;

        irep=find(diff(DataNoDuplicates.time)==0);
    end
    

    %fprintf('   -->%i DUPLICATE POINTS REMOVED (%3.3f s) \n',length(Data.time) - length(DataNoDuplicates.time),elapsedTime);
    
else
    error(['Fieldname ''time'' must be in data structure ' structname]);
end
    
end

