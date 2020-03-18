clc
clear all
close all
% *************************************************************************
% main_computeFeatures_Sway.m
%
% Load extracted data from main_Data_Extract.m. 
%
% Run after main_Data_Extract.m. 
%
% Megan O'Brien
% Version date 03-18-2020
% *************************************************************************

%% ------------------------------- INPUT ----------------------------------
% Select data to load

Subject_type = 'cva';    % cva or controls
Activity = 'BBS';     % Select one from list of nicknames above
Subj2Load = [1:55];   % Subject numbers to load (1:55 cva, 1:51 controls)

%% ---------------------------- INITIALIZE --------------------------------
% Directories
dirMfile = pwd;
dirData = ['Z:\Inpatient Sensors -Stroke\MC10 Study\Data analysis\2_Clean_Data_Extracted\' Activity filesep];  % Segmented sensor data

%********** TO DO ********** 
% Pair Score with Sensor Data in main_Data_Extract
%***************************

%% ---------------------------- MAIN LOOP ---------------------------------

% SUBJECTS
for indSub = 1:length(Subj2Load)
    
    filenameData = [upper(Subject_type) '_' Activity '_ID' sprintf('%02d',Subj2Load(indSub)) '.mat'];
    
    Subject_num = Subj2Load(indSub);
    Subject = get_subjectId(Subject_num,Subject_type);
    
    fprintf('* %s \n',Subject)
    
    % 1. ----------------- LOAD CLEAN, EXTRACTED DATA ---------------------
    load([dirData filenameData]);
    
    numDates = numel(data.Session_date);

    for indDate = 1:numDates
        numTrials = numel(data.Session_trials{indDate});
        
        %figure; 
        for indTrial = 1:numTrials
            [~,trialname] = strtok(data.Session_trials{indDate}(indTrial),'_');
            trialname = trialname{1}(2:end);
            trialname(regexp(trialname,'_'))=[];
            
            t = data.Session{indDate}.Motion.Time{indTrial};
            AccData_sacrum = data.Session{indDate}.Motion.SC.Acc{indTrial};
            Acc.x = AccData_sacrum(:,1);
            Acc.y = AccData_sacrum(:,2);
            Acc.z = AccData_sacrum(:,3);
            
            GyrData_sacrum = data.Session{indDate}.Motion.SC.Gyr{indTrial};
            Gyr.x = GyrData_sacrum(:,1);
            Gyr.y = GyrData_sacrum(:,2);
            Gyr.z = GyrData_sacrum(:,3);
            
            %subplot(numTrials,1,indTrial); hold on;
            %plot(t,AccData_sacrum); 
            
            temp.Acc = Acc;
            temp.Gyr = Gyr;
            
            dataTransformed = acctransformation(temp,Subject,indDate);
            
            %plot(t,dataTransformed.AP,'k-','LineWidth',1.5);
            %plot(t,dataTransformed.ML,'m-','LineWidth',1.5);
            %plot(t,dataTransformed.V,'b-','LineWidth',1);
            
            %plot(t,dataTransformed.angX,'k-','LineWidth',1.5); % Pelvic Tilt (forward/backward)
            %plot(t,dataTransformed.angY,'m-','LineWidth',1.5); % Pelvic Rotation
            %plot(t,dataTransformed.angZ,'b-','LineWidth',1);   % Pelvic Obliquity
            
            %title(trialname)
            
            AllFeatures.(trialname){indSub,indDate} = staticBalance(dataTransformed,t);
            
        end %indTrial
    end %indDate
            
end %indSub

% ********** TO DO ************
% Add Gyr features
% 
% Save features to data struct. Or add Session_Date and metadata to feature struct
% *****************************

%% ANALYSIS

% load('features_sway_BBS.mat')

C = [AllFeatures.STANDUNSUPPORTEDS1];
for indSub=1:length(Subj2Load)
    % Find first available datapoint (Admission time-point if they could
    % complete test)
    col_adm = find(~cellfun(@isempty,C(indSub,:)),1,'first');
    if isempty(col_adm) % No time-points available, so admission is NaN
        Feat_adm(indSub) = NaN;
    else
        Feat_adm(indSub) = C{indSub,col_adm}.data.length_swayAPAcc;
    end

    
    % Find last available datapoint (Discharge time-point)
    col_dis = find(~cellfun(@isempty,C(indSub,:)),1,'last');
    if isempty(col_dis) || col_dis == 1 % No time-points available, or Single time-point only, so discharge is NaN
        Feat_dis(indSub) = NaN;
    else
        Feat_dis(indSub) = C{indSub,col_dis}.data.length_swayAPAcc;
    end
end

figure;
boxplot([Feat_adm' Feat_dis'],'Labels',{'Adm','Dis'})
