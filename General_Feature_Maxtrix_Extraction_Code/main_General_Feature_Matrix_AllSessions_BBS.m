%-------------------------------------------------------
%
% BBS Feature Extraction for Inpatient Sensor Data
% Adam P. Horin
% June 23, 2020
%
% Functions needed:
%   FFeatures.m
%   FFT_data.m
%   acctransformation.m
%   ellipsoid.m
%   sampen.m
%   staticBalance.m
%
% Description:
%   1. Loops through clean extracted data of BBS inpatient sensor data for
%   CVA and HC and extracts general and balance features
%   2. Extracts features for all 14 activities of the BBS at 15 second
%   cutoffs (15-120, and 200 seconds)
%   3. Extracts BBS scores for all activities and the total score for each
%   participant
%
%------------------------------------------------------------------

clc
clear all
close all

%% Define variables
Type_of_Subject = {'CVA', 'CONTROLS'}; % enter CONTROLS and/or CVA
Activity = 'BBS'; %Activity Selection

Hz = 31.25;
dt = 1/Hz;
k = 1;
BBS_table = [];
TrialNames = {'N1_SIT_TO_STAND_S1', ...
    'N2_STAND_UNSUPPORTED_S1', ...
    'N3_SIT_W__BACK_S1', ...
    'N4_STAND_TO_SIT_S1', ...
    'N5_TRANSFERS_S1', ...
    'N6_STAND_W__EYES_CLOSED_S1', ...
    'N7_STAND_W__FEET_TOGETHE_Tr_S1', ...
    'N8_REACH_FORWARD_S1', ...
    'N9_PICK_UP_OBJECT_S1', ...
    'N10_TURN_S1', ...
    'N11_TURN_360_DEGS_S1', ...
    'N12_FOOT_ON_STEP_STOOL_S1', ...
    'N13_STAND_W__ONE_FOOT_IN_Tr_S1', ...
    'N14_STAND_ON_ONE_LEG_S1'};

TrialNumbers = {'N01', ...
    'N02', ...
    'N03', ...
    'N04', ...
    'N05', ...
    'N06', ...
    'N07', ...
    'N08', ...
    'N09', ...
    'N10', ...
    'N11', ...
    'N12', ...
    'N13', ...
    'N14'};

% cuttoff times for processing; 200 is used to get all of the data from the
% trial; no activity should exceed more than 122 seconds
TrialTimes = [15, 30, 45, 60, 75, 90, 105, 120, 200];

%directory
path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Data analysis\2_Clean_Data_Extracted\';

%% Main Loop: identifies group, subject ID, trial names, and then trial times

for h = 1:1:length(Type_of_Subject) % Identify group to loop through
    if strcmp(Type_of_Subject{h}, 'CVA') == 1
        ID = [1:55]; %[1:55]
        Type_of_Subject_Group = 'CVA';
        BBS_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_BBS.xlsx';
    elseif strcmp(Type_of_Subject{h}, 'CONTROLS') == 1
        ID = [1:51]; %[1:51]
        Type_of_Subject_Group = 'CONTROLS';
        BBS_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_BBS.xlsx';
    end
    
    for n = 1:1:length(ID) % Loop through IDs for the specified group
        file_input = [path Activity '\' Type_of_Subject_Group '_' Activity '_ID' sprintf('%.2d',ID(n)) '.mat'];
        load(file_input);
        data
        
        %Get Subject ID
        if isnumeric(ID(n))
            str = num2str(ID(n));
            if length(str) == 1
                if strcmp(Type_of_Subject_Group,'CVA')
                    SubjectID=['CVA0' str];
                else strcmp(Type_of_Subject_Group,'CONTROLS')
                    SubjectID=['HC0' str];
                end
            else
                if strcmp(Type_of_Subject_Group,'CVA')
                    SubjectID=['CVA' str];
                else strcmp(Type_of_Subject_Group,'CONTROLS')
                    SubjectID=['HC' str];
                end
            end
        end
        
        
%-------------------------------------------------------------------------%
% Check Dates for Session Numbers
%
% This section gets the index of the session date based on the number of
% sessions. If there are only 2 or 3 sessions the last session is changed
% to session 4. The dates orders are then checked and the index is switched
% if necessary.
%
%-------------------------------------------------------------------------%
        % first get the value of each date and define the index
        SN1_start=1;
        for k = 1:1:length(data.Session_date)
            if length(data.Session_date) == 1
                SN1 = 1;
                SN2 = [];
                SN3 = [];
                SN4 = [];
            elseif length(data.Session_date) == 2
                SN1 = 1;
                SN2 = 2;
                SN3 = [];
                SN4 = [];
            elseif length(data.Session_date) == 3
                SN1 = 1;
                SN2 = 2;
                SN3 = 3;
                SN4 = [];
            elseif length(data.Session_date) == 4
                SN1 = 1;
                SN2 = 2;
                SN3 = 3;
                SN4 = 4;
            end
            
        end
        
        SN = {SN1, SN2, SN3, SN4};

        % check to make sure the the dates are in the correct order
        % make the last session SN4
        
        if isempty(SN2) == 1 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(data.Session_date(1))];
            order = 1;
        elseif isempty(SN2) == 0 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(data.Session_date(1)), datenum(data.Session_date(2))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                SN4 = SN2; %makes the discharge session SN4
                SN2 = [];
                order = 1;
            else
                %dates are not in correct order
                %date_index = find(contains(dates(),dates_sorted(1))))
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                SN1 = SN{date_index1};
                SN2 = [];
                SN3 = [];
                SN4 = SN{date_index2};
                order = 0;
            end
        elseif isempty(SN2) == 0 && isempty(SN3) == 0 && isempty(SN4) == 1
            
            dates = [datenum(data.Session_date(1)), datenum(data.Session_date(2)), datenum(data.Session_date(3))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                SN4 = SN3;
                SN3 = [];
                order = 1;
            else
                %dates are not in correct order
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                date_index3 = find(ismember(dates, dates_sorted(3)));
                SN1 = SN{date_index1};
                SN2 = SN{date_index2};
                SN3 = [];
                SN4 = SN{date_index3};
                order = 0;
                
            end
        elseif isempty(SN2) == 0 && isempty(SN3) == 0 && isempty(SN4) == 0
            dates = [datenum(data.Session_date(1)), datenum(data.Session_date(2)), datenum(data.Session_date(3)), datenum(data.Session_date(4))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                order = 1;
            else
                %dates are not in correct order
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                date_index3 = find(ismember(dates, dates_sorted(3)));
                date_index4 = find(ismember(dates, dates_sorted(4)));
                SN1 = SN{date_index1};
                SN2 = SN{date_index2};
                SN3 =  SN{date_index3};
                SN4 =  SN{date_index4};
                order = 0;
            end
        end
        
        % Redefine SN in case dates were out of order
        SN = {SN1, SN2, SN3, SN4}
        
     
% Extract features for each session and trial
        for s = 1:1:length(SN)  
            for j = 1:1:length(TrialNames) % Loop through all BBS activities
                %Trial N
                if isfield(data, 'Session_trials') == 0 % check for data; if empty participant will get NaN
                    N_index = [];
                else
                    N_index = find(strcmp(data.Session_trials{SN{1}}, TrialNames(j)));
                end
                
                for g = 1:1:length(TrialTimes) % Export features at different time cutoffs
                    if isempty(N_index) == 1 || isempty(SN{s}) == 1
                        %Trial Information
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        activity(g,:) = Activity;
                        session(g,:) = s;
                        trial_No(g,:) = TrialNumbers(j);
                        cutoff(g,:) = TrialTimes(g);
                        trialtime(g,:) = nan;
                        %                         BBS_TotalScore(g,:) = sum(BBS(1:14));
                        %                         BBS_Subscore(g,:) = BBS(j);
                        
                        
                        %General Features
                        %Gyro mean
                        SC_Gyr_x_mean(g,:) = nan;
                        SC_Gyr_y_mean(g,:) = nan;
                        SC_Gyr_z_mean(g,:) = nan;
                        SC_Gyr_norm_mean(g,:) = nan;
                        
                        % Range
                        SC_Gyr_x_range(g,:) = nan;
                        SC_Gyr_y_range(g,:) = nan;
                        SC_Gyr_z_range(g,:) = nan;
                        SC_Gyr_norm_range(g,:) = nan;
                        
                        % RMS
                        SC_Gyr_x_rms(g,:) = nan;
                        SC_Gyr_y_rms(g,:) = nan;
                        SC_Gyr_z_rms(g,:) = nan;
                        SC_Gyr_norm_rms(g,:) = nan;
                        
                        % Standard Deviation
                        SC_Gyr_x_std(g,:) = nan;
                        SC_Gyr_y_std(g,:) = nan;
                        SC_Gyr_z_std(g,:) = nan;
                        SC_Gyr_norm_std(g,:) = nan;
                        
                        % Skew
                        SC_Gyr_x_skew(g,:) = nan;
                        SC_Gyr_y_skew(g,:) = nan;
                        SC_Gyr_z_skew(g,:) = nan;
                        SC_Gyr_norm_skew(g,:) = nan;
                        
                        % Kurtosis
                        SC_Gyr_x_kurtosis(g,:) = nan;
                        SC_Gyr_y_kurtosis(g,:) = nan;
                        SC_Gyr_z_kurtosis(g,:) = nan;
                        SC_Gyr_norm_kurtosis(g,:) = nan;
                        
                        % Derivative
                        % Gyro mean
                        dSC_Gyr_x_mean(g,:) = nan;
                        dSC_Gyr_y_mean(g,:) = nan;
                        dSC_Gyr_z_mean(g,:) = nan;
                        dSC_Gyr_norm_mean(g,:) = nan;
                        
                        % Range
                        dSC_Gyr_x_range(g,:) = nan;
                        dSC_Gyr_y_range(g,:) = nan;
                        dSC_Gyr_z_range(g,:) = nan;
                        dSC_Gyr_norm_range(g,:) = nan;
                        
                        % RMS
                        dSC_Gyr_x_rms(g,:) = nan;
                        dSC_Gyr_y_rms(g,:) = nan;
                        dSC_Gyr_z_rms(g,:) = nan;
                        dSC_Gyr_norm_rms(g,:) = nan;
                        
                        % Standard Deviation
                        dSC_Gyr_x_std(g,:) = nan;
                        dSC_Gyr_y_std(g,:) = nan;
                        dSC_Gyr_z_std(g,:) = nan;
                        dSC_Gyr_norm_std(g,:) = nan;
                        
                        % Skew
                        dSC_Gyr_x_skew(g,:) = nan;
                        dSC_Gyr_y_skew(g,:) = nan;
                        dSC_Gyr_z_skew(g,:) = nan;
                        dSC_Gyr_norm_skew(g,:) = nan;
                        
                        % Kurtosis
                        dSC_Gyr_x_kurtosis(g,:) = nan;
                        dSC_Gyr_y_kurtosis(g,:) = nan;
                        dSC_Gyr_z_kurtosis(g,:) = nan;
                        dSC_Gyr_norm_kurtosis(g,:) = nan;
                        
                        % Pearson correlation coefficient
                        SC_Gyr_corr_xy(g,:) = nan;
                        SC_Gyr_corr_xz(g,:) = nan;
                        SC_Gyr_corr_yz(g,:) = nan;
                        
                        % Sample Entropy
                        SC_Gyr_x_SamEn(g,:) = nan;
                        SC_Gyr_y_SamEn(g,:) = nan;
                        SC_Gyr_z_SamEn(g,:) = nan;
                        SC_Gyr_norm_SamEn(g,:) = nan;
                        
                        % Frequency Domain
                        SC_Gyr_x_DAmp(g,:) = nan;
                        SC_Gyr_x_DFreq(g,:) = nan;
                        SC_Gyr_x_PSD_mean(g,:) = nan;
                        SC_Gyr_x_PSD_std(g,:) = nan;
                        SC_Gyr_x_PSD_skew(g,:) = nan;
                        SC_Gyr_x_PSD_kurtosis(g,:) = nan;
                        
                        SC_Gyr_y_DAmp(g,:) = nan;
                        SC_Gyr_y_DFreq(g,:) = nan;
                        SC_Gyr_y_PSD_mean(g,:) = nan;
                        SC_Gyr_y_PSD_std(g,:) = nan;
                        SC_Gyr_y_PSD_skew(g,:) = nan;
                        SC_Gyr_y_PSD_kurtosis(g,:) = nan;
                        
                        SC_Gyr_z_DAmp(g,:) = nan;
                        SC_Gyr_z_DFreq(g,:) = nan;
                        SC_Gyr_z_PSD_mean(g,:) = nan;
                        SC_Gyr_z_PSD_std(g,:) = nan;
                        SC_Gyr_z_PSD_skew(g,:) = nan;
                        SC_Gyr_z_PSD_kurtosis(g,:) = nan;
                        
                        SC_Gyr_norm_DAmp(g,:) = nan;
                        SC_Gyr_norm_DFreq(g,:) = nan;
                        SC_Gyr_norm_PSD_mean(g,:) = nan;
                        SC_Gyr_norm_PSD_std(g,:) = nan;
                        SC_Gyr_norm_PSD_skew(g,:) = nan;
                        SC_Gyr_norm_PSD_kurtosis(g,:) = nan;
                        
                        
                    else
                        %Trial Information
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        activity(g,:) = Activity;
                        session(g,:) = s;
                        trial_No(g,:) = TrialNumbers(j);
                        cutoff(g,:) = TrialTimes(g);
                        %                         BBS_TotalScore(g,:) = sum(BBS(1:14));
                        %                         BBS_Subscore(g,:) = BBS(j);
                        
                        
                        % Gyroscope data
                        SC_Gyr = data.Session{SN{s}}.Motion.SC.Gyr{N_index};
                        
                        % Acceleration data
                        SC_Acc = data.Session{SN{s}}.Motion.SC.Acc{N_index};
                        
                        % Demean Acc data
                        SC_Acc = SC_Acc - ones(length(SC_Acc),1)*mean(SC_Acc);
                        
                        % To match final time
                        
                        %final = min([length(SC_Gyr) length(SC_Acc)]);
                        final = round(Hz*TrialTimes(g));
                        if final >= length(data.Session{SN{s}}.Motion.Time{N_index}(1:end,:))
                            Time = data.Session{SN{s}}.Motion.Time{N_index}(1:end,:);
                            SC_Gyr = SC_Gyr(1:end,:);
                            SC_Acc = SC_Acc(1:end,:);
                            trialtime(g,:) = data.Session{SN{s}}.Motion.Time{N_index}(end,:)-data.Session{SN{s}}.Motion.Time{N_index}(1,:);
                        else
                            Time = data.Session{SN{s}}.Motion.Time{N_index}(1:final,:);
                            SC_Gyr = SC_Gyr(1:final,:);
                            SC_Acc = SC_Acc(1:final,:);
                            trialtime(g,:) = data.Session{SN{s}}.Motion.Time{N_index}(final,:)-data.Session{SN{s}}.Motion.Time{N_index}(1,:);
                        end
                        
                        
                        % Norm of Gyro and Acc Data
                        
                        for i = 1:1:length(Time)
                            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
                            SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
                            
                        end
                        
                        %General Features
                        %Gyro Mean
                        SC_Gyr_x_mean(g,:) = mean(SC_Gyr(:,1));
                        SC_Gyr_y_mean(g,:) = mean(SC_Gyr(:,2));
                        SC_Gyr_z_mean(g,:) = mean(SC_Gyr(:,3));
                        SC_Gyr_norm_mean(g,:) = mean(SC_Gyr_norm);
                        
                        % Range
                        SC_Gyr_x_range(g,:) = range(SC_Gyr(:,1));
                        SC_Gyr_y_range(g,:) = range(SC_Gyr(:,2));
                        SC_Gyr_z_range(g,:) = range(SC_Gyr(:,3));
                        SC_Gyr_norm_range(g,:) = range(SC_Gyr_norm);
                        
                        % RMS
                        SC_Gyr_x_rms(g,:) = rms(SC_Gyr(:,1));
                        SC_Gyr_y_rms(g,:) = rms(SC_Gyr(:,2));
                        SC_Gyr_z_rms(g,:) = rms(SC_Gyr(:,3));
                        SC_Gyr_norm_rms(g,:) = rms(SC_Gyr_norm);
                        
                        % Standard Deviation
                        SC_Gyr_x_std(g,:) = std(SC_Gyr(:,1));
                        SC_Gyr_y_std(g,:) = std(SC_Gyr(:,2));
                        SC_Gyr_z_std(g,:) = std(SC_Gyr(:,3));
                        SC_Gyr_norm_std(g,:) = std(SC_Gyr_norm);
                        
                        % Skew
                        SC_Gyr_x_skew(g,:) = skewness(SC_Gyr(:,1));
                        SC_Gyr_y_skew(g,:) = skewness(SC_Gyr(:,2));
                        SC_Gyr_z_skew(g,:) = skewness(SC_Gyr(:,3));
                        SC_Gyr_norm_skew(g,:) = skewness(SC_Gyr_norm);
                        
                        % Kurtosis
                        SC_Gyr_x_kurtosis(g,:) = kurtosis(SC_Gyr(:,1));
                        SC_Gyr_y_kurtosis(g,:) = kurtosis(SC_Gyr(:,2));
                        SC_Gyr_z_kurtosis(g,:) = kurtosis(SC_Gyr(:,3));
                        SC_Gyr_norm_kurtosis(g,:) = kurtosis(SC_Gyr_norm);
                        
                        % Derivative
                        % Gyro mean
                        dSC_Gyr_x_mean(g,:) = mean(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_mean(g,:) = mean(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_mean(g,:) = mean(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_mean(g,:) = mean(diff(SC_Gyr_norm)/dt);
                        
                        % Range
                        dSC_Gyr_x_range(g,:) = range(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_range(g,:) = range(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_range(g,:) = range(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_range(g,:) = range(diff(SC_Gyr_norm)/dt);
                        
                        % RMS
                        dSC_Gyr_x_rms(g,:) = rms(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_rms(g,:) = rms(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_rms(g,:) = rms(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_rms(g,:) = rms(diff(SC_Gyr_norm)/dt);
                        
                        % Standard Deviation
                        dSC_Gyr_x_std(g,:) = std(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_std(g,:) = std(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_std(g,:) = std(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_std(g,:) = std(diff(SC_Gyr_norm)/dt);
                        
                        % Skew
                        dSC_Gyr_x_skew(g,:) = skewness(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_skew(g,:) = skewness(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_skew(g,:) = skewness(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_skew(g,:) = skewness(diff(SC_Gyr_norm)/dt);
                        
                        % Kurtosis
                        dSC_Gyr_x_kurtosis(g,:) = kurtosis(diff(SC_Gyr(:,1))/dt);
                        dSC_Gyr_y_kurtosis(g,:) = kurtosis(diff(SC_Gyr(:,2))/dt);
                        dSC_Gyr_z_kurtosis(g,:) = kurtosis(diff(SC_Gyr(:,3))/dt);
                        dSC_Gyr_norm_kurtosis(g,:) = kurtosis(diff(SC_Gyr_norm)/dt);
                        
                        % Pearson correlation coefficient
                        SC_Gyr_corr_xy(g,:) = corr(SC_Gyr(:,1),SC_Gyr(:,2));
                        SC_Gyr_corr_xz(g,:) = corr(SC_Gyr(:,1),SC_Gyr(:,3));
                        SC_Gyr_corr_yz(g,:) = corr(SC_Gyr(:,2),SC_Gyr(:,3));
                        
                        % Sample Entropy
                        r = 0.2;
                        SC_Gyr_x_SamEn(g,:) = sampen(SC_Gyr(:,1),1,r);
                        SC_Gyr_y_SamEn(g,:) = sampen(SC_Gyr(:,2),1,r);
                        SC_Gyr_z_SamEn(g,:) = sampen(SC_Gyr(:,3),1,r);
                        SC_Gyr_norm_SamEn(g,:) = sampen(SC_Gyr_norm,1,r);
                        
                        % Frequency Domain
                        ff = FFeatures(SC_Gyr(:,1), Hz);
                        SC_Gyr_x_DAmp(g,:) = ff(1);
                        SC_Gyr_x_DFreq(g,:) = ff(2);
                        SC_Gyr_x_PSD_mean(g,:) = ff(3);
                        SC_Gyr_x_PSD_std(g,:) = ff(4);
                        SC_Gyr_x_PSD_skew(g,:) = ff(5);
                        SC_Gyr_x_PSD_kurtosis(g,:) = ff(6);
                        
                        ff = FFeatures(SC_Gyr(:,2), Hz);
                        SC_Gyr_y_DAmp(g,:) = ff(1);
                        SC_Gyr_y_DFreq(g,:) = ff(2);
                        SC_Gyr_y_PSD_mean(g,:) = ff(3);
                        SC_Gyr_y_PSD_std(g,:) = ff(4);
                        SC_Gyr_y_PSD_skew(g,:) = ff(5);
                        SC_Gyr_y_PSD_kurtosis(g,:) = ff(6);
                        
                        ff = FFeatures(SC_Gyr(:,3), Hz);
                        SC_Gyr_z_DAmp(g,:) = ff(1);
                        SC_Gyr_z_DFreq(g,:) = ff(2);
                        SC_Gyr_z_PSD_mean(g,:) = ff(3);
                        SC_Gyr_z_PSD_std(g,:) = ff(4);
                        SC_Gyr_z_PSD_skew(g,:) = ff(5);
                        SC_Gyr_z_PSD_kurtosis(g,:) = ff(6);
                        
                        ff = FFeatures(SC_Gyr_norm, Hz);
                        SC_Gyr_norm_DAmp(g,:) = ff(1);
                        SC_Gyr_norm_DFreq(g,:) = ff(2);
                        SC_Gyr_norm_PSD_mean(g,:) = ff(3);
                        SC_Gyr_norm_PSD_std(g,:) = ff(4);
                        SC_Gyr_norm_PSD_skew(g,:) = ff(5);
                        SC_Gyr_norm_PSD_kurtosis(g,:) = ff(6);
                        
                        % To compute balance features
                        %Code from main_computeFeatures_Sway
                        t = Time;
                        AccData_sacrum = data.Session{SN{s}}.Motion.SC.Acc{N_index};
                        Acc.x = SC_Acc(:,1);
                        Acc.y = SC_Acc(:,2);
                        Acc.z = SC_Acc(:,3);
                        
                        GyrData_sacrum = data.Session{SN{s}}.Motion.SC.Gyr{N_index};
                        Gyr.x = SC_Gyr(:,1);
                        Gyr.y = SC_Gyr(:,2);
                        Gyr.z = SC_Gyr(:,3);
                        
                        temp.Acc = Acc;
                        temp.Gyr = Gyr;
                        
                        transformData = acctransformation(temp,ID(n),SN{s});
                        name = char(TrialNames(g));
                        AllFeatures.(name){1,SN{s}} = staticBalance(transformData,t);
                        
                        %Sway Features
                        f50_ML(g,:) = AllFeatures.(name){1,1}.data.f50_ML;
                        f50_AP(g,:) = AllFeatures.(name){1,1}.data.f50_AP;
                        f95_ML(g,:) = AllFeatures.(name){1,1}.data.f95_ML;
                        f95_AP(g,:) = AllFeatures.(name){1,1}.data.f95_AP;
                        spectral_centroid_AP(g,:) = AllFeatures.(name){1,1}.data.spectral_centroid_AP;
                        spectral_centroid_ML(g,:) = AllFeatures.(name){1,1}.data.spectral_centroid_ML;
                        max_accAP(g,:) = AllFeatures.(name){1,1}.data.max_accAP;
                        max_accML(g,:) = AllFeatures.(name){1,1}.data.max_accML;
                        mean_accAP(g,:) = AllFeatures.(name){1,1}.data.mean_accAP;
                        mean_accML(g,:) = AllFeatures.(name){1,1}.data.mean_accML;
                        rms_AP(g,:) = AllFeatures.(name){1,1}.data.rms_AP;
                        rms_ML(g,:) = AllFeatures.(name){1,1}.data.rms_ML;
                        jerk_AP(g,:) = AllFeatures.(name){1,1}.data.jerk_AP;
                        jerk_ML(g,:) = AllFeatures.(name){1,1}.data.jerk_ML;
                        mean_velAP(g,:) = AllFeatures.(name){1,1}.data.mean_velAP;
                        mean_velML(g,:) = AllFeatures.(name){1,1}.data.mean_velML;
                        length_swayAPAcc(g,:) = AllFeatures.(name){1,1}.data.length_swayAPAcc;
                        length_swayMLAcc(g,:) = AllFeatures.(name){1,1}.data.length_swayMLAcc;
                        
                    end
                end
                
                
                %BBS feature matrix in long format
                BBS_table_temp = table(subject, group, activity, session, trial_No, cutoff, trialtime, ...
                    SC_Gyr_x_mean, SC_Gyr_y_mean, SC_Gyr_z_mean, SC_Gyr_norm_mean, ...
                    SC_Gyr_x_range, SC_Gyr_y_range, SC_Gyr_z_range, SC_Gyr_norm_range, ...
                    SC_Gyr_x_rms, SC_Gyr_y_rms, SC_Gyr_z_rms, SC_Gyr_norm_rms, ...
                    SC_Gyr_x_std, SC_Gyr_y_std, SC_Gyr_z_std, SC_Gyr_norm_std, ...
                    SC_Gyr_x_skew, SC_Gyr_y_skew, SC_Gyr_z_skew, SC_Gyr_norm_skew, ...
                    SC_Gyr_x_kurtosis, SC_Gyr_y_kurtosis, SC_Gyr_z_kurtosis, SC_Gyr_norm_kurtosis, ...
                    dSC_Gyr_x_mean, dSC_Gyr_y_mean, dSC_Gyr_z_mean, dSC_Gyr_norm_mean, ...
                    dSC_Gyr_x_range, dSC_Gyr_y_range, dSC_Gyr_z_range, dSC_Gyr_norm_range, ...
                    dSC_Gyr_x_rms, dSC_Gyr_y_rms, dSC_Gyr_z_rms, dSC_Gyr_norm_rms, ...
                    dSC_Gyr_x_std, dSC_Gyr_y_std, dSC_Gyr_z_std, dSC_Gyr_norm_std, ...
                    dSC_Gyr_x_skew, dSC_Gyr_y_skew, dSC_Gyr_z_skew, dSC_Gyr_norm_skew, ...
                    dSC_Gyr_x_kurtosis, dSC_Gyr_y_kurtosis, dSC_Gyr_z_kurtosis, dSC_Gyr_norm_kurtosis, ...
                    SC_Gyr_corr_xy, SC_Gyr_corr_xz, SC_Gyr_corr_yz, ...
                    SC_Gyr_x_SamEn, SC_Gyr_y_SamEn, SC_Gyr_z_SamEn, SC_Gyr_norm_SamEn, ...
                    SC_Gyr_x_DAmp, SC_Gyr_x_DFreq, SC_Gyr_x_PSD_mean, SC_Gyr_x_PSD_std, SC_Gyr_x_PSD_skew, SC_Gyr_x_PSD_kurtosis, ...
                    SC_Gyr_y_DAmp, SC_Gyr_y_DFreq, SC_Gyr_y_PSD_mean, SC_Gyr_y_PSD_std, SC_Gyr_y_PSD_skew, SC_Gyr_y_PSD_kurtosis, ...
                    SC_Gyr_z_DAmp, SC_Gyr_z_DFreq, SC_Gyr_z_PSD_mean, SC_Gyr_z_PSD_std, SC_Gyr_z_PSD_skew, SC_Gyr_z_PSD_kurtosis, ...
                    SC_Gyr_norm_DAmp, SC_Gyr_norm_DFreq, SC_Gyr_norm_PSD_mean, SC_Gyr_norm_PSD_std, SC_Gyr_norm_PSD_skew, SC_Gyr_norm_PSD_kurtosis, ...
                    f50_ML, f50_AP, f95_ML, f95_AP, spectral_centroid_AP, spectral_centroid_ML, max_accAP,...
                    max_accML, mean_accAP, mean_accML, rms_AP, rms_ML, jerk_AP, jerk_ML,...
                    mean_velAP, mean_velML, length_swayAPAcc, length_swayMLAcc);
                BBS_table = [BBS_table; BBS_table_temp]
                
                % Check if any date orders needed to be changed
                if order(1) == 0
                    DateOrderID = {SubjectID};
                    WrongSession = s;
                    List_WrongDate_temp = table(DateOrderID, WrongSession);
                    List_WrongDate = [List_WrongDate; List_WrongDate_temp];
                end
            end
        end
    end
end

%% Export .csv of feature matrix

writetable(BBS_table,'General_Feature_Matrix_AllSessions_BBS.csv','Delimiter',',','QuoteStrings',true)

