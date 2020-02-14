function [ AnnotationsOut ] = get_annotationind( Task, Activity, Labels )
%**************************************************************************
% get_annotationind.m
%
% Written by Megan O'Brien
% version 6/1/2018
%
% Get row index from annotation file for task of interest.
%
% Task = list of tasks within a Biostamp activity label (e.g. Hip Flexion, 
% Hip Extension, and Knee Extension in Modified Ashworth test)
%
% indIn = list of all indices for a Biostamp activity label (e.g. all
% instances of Modified Ashworth)
%
% Needs to find Right and Left for relevant tasks (MAS, MMT)
%**************************************************************************

AllActivities = Labels(:,3);    % get names of all Activities
indAct = strcmp(AllActivities,Activity);  % get all indices of desired Activity (e.g Clinical - BBS)

TimeEnd = Labels(indAct,6);  % End time is same for Activity and all Survey Questions

% Full block of annotations related to one activity label (e.g. side, type, etc.)
AllActivityInfo = cell(length(TimeEnd),3);
for i=1:length(TimeEnd)
    indBlock = find(strcmp(Labels(:,6),TimeEnd(i)));  
    
    AllActivityInfo{i,1} = indBlock; % indices for this block
    
    % Get additional annotation info depending on type of test
    % ***** MMT *****
    if strcmp(Activity,'Clinical - MMT')
        indType = find(strcmp(Labels(indBlock,3),'Type of Activity'));
        type = Labels(indBlock(indType(end)),7);
        
        indSide = find(strcmp(Labels(indBlock,3),'Side'));
        side = Labels(indBlock(indSide(end)),7);
        
        AllActivityInfo{i,2} = type; % Hip Flexion, Hip Extension, etc.
        AllActivityInfo{i,3} = side; % Right or Left
    end
    
    % ***** MAS *****
    if strcmp(Activity,'Clinical - MAS')
        indType = find(strcmp(Labels(indBlock,3),'Movement Type'));
        type = Labels(indBlock(indType(end)),7);
        
        indSide = find(strcmp(Labels(indBlock,3),'Side'));
        side = Labels(indBlock(indSide(end)),7);
        
        AllActivityInfo{i,2} = type; % Hip Flexion, Hip Extension, etc.
        AllActivityInfo{i,3} = side; % Right or Left
    end
    
    % ***** BBS *****
    if strcmp(Activity,'Clinical - BBS')
        indType = find(strcmp(Labels(indBlock,3),'Type of Assessment'));
        type = Labels(indBlock(indType(end)),7);
        
        AllActivityInfo{i,2} = type; % Sit to Stand, Transfers, etc.
    end
    
    % ***** Activity Recognition *****
    if strcmp(Activity,'Activity Recognition')
        indType = find(strcmp(Labels(indBlock,3),'Activity type'));
        type = Labels(indBlock(indType(end)),7);
        
        AllActivityInfo{i,2} = type; % Sitting, Standing, etc.
    end
    
    % ***** Physical Therapy *****
    if strcmp(Activity,'Physical Therapy')
        indType = find(strcmp(Labels(indBlock,3),'Type of Physical Therapy'));
        type = Labels(indBlock(indType(end)),7);
        
        AllActivityInfo{i,2} = type; % Walking on a Treadmill, etc.
    end
    
end

% Find annotations relating to task of interest (e.g. Hip Flexion only) for
% activities with sub-surveys. Otherwise, all annotations (i.e. walking
% tasks)
AnnotationsOut = [];
if strcmp(Activity,'Clinical - 10MWT SSV') || strcmp(Activity,'Clinical - 10MWT FV') || ...
        strcmp(Activity,'Clinical - 6MWT') || strcmp(Activity,'Clinical - TUG')
    AnnotationsOut = AllActivityInfo;
else
    for i=1:length(Task)
        indTask = strcmp([AllActivityInfo{:,2}],Task(i));
        AnnotationsOut = [AnnotationsOut; AllActivityInfo(indTask,:)];
    end
end

end