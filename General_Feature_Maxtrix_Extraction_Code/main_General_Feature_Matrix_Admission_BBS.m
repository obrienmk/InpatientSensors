%------------------------%
%
% BBS Feature Matrix
% Adam P. Horin
% 06/08/2020
%
%-------------------------%

clc
clear all
close all

ID = [1:55]; %enter ID or range of ID's to be processed
Type_of_Subject = 'CVA'; % enter CONTROLS or CVA

%Activity Selection
Activity = 'BBS';

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

%directory
path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Data analysis\2_Clean_Data_Extracted\';

SN=1; %admission session

for n = 1:1:length(ID)
    
    file_input = [path Activity '\' Type_of_Subject '_' Activity '_ID' sprintf('%.2d',ID(n)) '.mat'];
    load(file_input);
    data
    
    for g = 1:1:length(TrialNames)
        %Trial N1
        N_index = find(strcmp(data.Session_trials{SN}, TrialNames(g)));
        if isempty(N_index) == 1
            subject(g,:) = ID(n);
            group(g,:) = Type_of_Subject;
            activity(g,:) = Activity;
            trial_No(g,:) = TrialNumbers(g);
            
            SC_Gyr_x_mean(g,:) = nan;
            SC_Gyr_y_mean(g,:) = nan;
            SC_Gyr_z_mean(g,:) = nan;
            trialtime(g,:) = nan;
        else
            subject(g,:) = ID(n);
            group(g,:) = Type_of_Subject;
            activity(g,:) = Activity;
            trial_No(g,:) = TrialNumbers(g);
            
            % Gyroscope data
            SC_Gyr = data.Session{SN}.Motion.SC.Gyr{N_index};
            
            % Acceleration data
            SC_Acc = data.Session{SN}.Motion.SC.Acc{N_index};
            
            % Demean Acc data
            SC_Acc = SC_Acc - ones(length(SC_Acc),1)*mean(SC_Acc);
            
            % To match final time
            final = min([length(SC_Gyr) length(SC_Acc)]);
            %         final = Hz*20;  % Initial 20 sec
            Time = data.Session{SN}.Motion.Time{N_index}(1:final,:);
            SC_Gyr = SC_Gyr(1:final,:);
            SC_Acc = SC_Acc(1:final,:);
            
            
            % Norm of Gyro and Acc Data
            for i = 1:1:length(Time)
                SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
                SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
            end
            
            SC_Gyr_x_mean(g,:) = mean(SC_Gyr(:,1));
            SC_Gyr_y_mean(g,:) = mean(SC_Gyr(:,2));
            SC_Gyr_z_mean(g,:) = mean(SC_Gyr(:,3));
            
            trialtime(g,:) = data.Session{SN}.Motion.Time{N_index}(final,:)-data.Session{SN}.Motion.Time{N_index}(1,:);
        end
        
        
    end
    %BBS feature matrix in long format
    BBS_table_temp = table(subject, group, activity, trial_No, SC_Gyr_x_mean, SC_Gyr_y_mean, SC_Gyr_z_mean, trialtime)
    BBS_table = [BBS_table; BBS_table_temp]
end

%tranpose feature matrix to wide format by trial number
variableList = {'SC_Gyr_x_mean','SC_Gyr_y_mean','SC_Gyr_z_mean', 'trialtime'};
BBS_table_wide = unstack(BBS_table,variableList,'trial_No')

writetable(BBS_table_wide,'General_Feature_Matrix_Admission_BBS.csv','Delimiter',',','QuoteStrings',true)



