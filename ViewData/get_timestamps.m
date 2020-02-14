function [ DataTime ] = get_timestamps( Labels )
%**************************************************************************
% get_timestamps.m
%
% Written by Megan O'Brien
% version 9.18.2018
%
% Convert epoch timestamps from tablet labels into a readable timestamp
%**************************************************************************

% Get date (yyyy-mm-dd) that each label was collected
LabelTimeUTC = cellfun(@str2double,Labels(:,1));
N=size(LabelTimeUTC);

DataTime=datetime([repmat(1970,N), ones(N), ones(N), zeros(N), zeros(N), LabelTimeUTC*0.001],'Format','yyyy-MM-dd'); % yyyy-mm-dd only
%DataTime=datetime([repmat(1970,N), ones(N), ones(N), zeros(N), zeros(N),LabelTimeUTC*0.001]);  % full timestamp
    
end

