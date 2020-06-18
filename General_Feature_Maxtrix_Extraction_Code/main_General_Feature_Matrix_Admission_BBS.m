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


ID = [1:2]; %enter ID or range of ID's to be processed CONTROLS [1:51]; CVA [1:55]
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
    
    %Get Subject ID
    if isnumeric(ID(n))
        str = num2str(ID(n));
        if length(str) == 1
            if strcmp(Type_of_Subject,'CVA')
                Subject=['CVA0' str];
            else strcmp(Type_of_Subject,'CONTROLS')
                Subject=['HC0' str];
            end
        else
            if strcmp(Type_of_Subject,'CVA')
                Subject=['CVA' str];
            else strcmp(Type_of_Subject,'CONTROLS')
                Subject=['HC' str];
            end
        end
    end
    
    for g = 1:1:length(TrialNames)
        
        %Trial N
        if isfield(data, 'Session_trials') == 0
            N_index = [];
        else
            N_index = find(strcmp(data.Session_trials{SN}, TrialNames(1)));
        end
        
        if isempty(N_index) == 1
            
            %Trial Information
            subject(g,:) = Subject;
            group(g,:) = Type_of_Subject;
            activity(g,:) = Activity;
            trial_No(g,:) = TrialNumbers(g);
            trialtime(g,:) = nan;
            
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
            subject(g,:) = Subject;
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
            
            
            trialtime(g,:) = data.Session{SN}.Motion.Time{N_index}(final,:)-data.Session{SN}.Motion.Time{N_index}(1,:);
        end
        
        
    end
    %BBS feature matrix in long format
    BBS_table_temp = table(subject, group, activity, trial_No, trialtime,...
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
        SC_Gyr_norm_DAmp, SC_Gyr_norm_DFreq, SC_Gyr_norm_PSD_mean, SC_Gyr_norm_PSD_std, SC_Gyr_norm_PSD_skew, SC_Gyr_norm_PSD_kurtosis)
    BBS_table = [BBS_table; BBS_table_temp]
end



%tranpose feature matrix to wide format by trial number
variableList = {'trialtime',...
    'SC_Gyr_x_mean', 'SC_Gyr_y_mean', 'SC_Gyr_z_mean', 'SC_Gyr_norm_mean', ...
    'SC_Gyr_x_range', 'SC_Gyr_y_range', 'SC_Gyr_z_range', 'SC_Gyr_norm_range', ...
    'SC_Gyr_x_rms', 'SC_Gyr_y_rms', 'SC_Gyr_z_rms', 'SC_Gyr_norm_rms', ...
    'SC_Gyr_x_std', 'SC_Gyr_y_std', 'SC_Gyr_z_std', 'SC_Gyr_norm_std', ...
    'SC_Gyr_x_skew', 'SC_Gyr_y_skew', 'SC_Gyr_z_skew', 'SC_Gyr_norm_skew', ...
    'SC_Gyr_x_kurtosis', 'SC_Gyr_y_kurtosis', 'SC_Gyr_z_kurtosis', 'SC_Gyr_norm_kurtosis', ...
    'dSC_Gyr_x_mean', 'dSC_Gyr_y_mean', 'dSC_Gyr_z_mean', 'dSC_Gyr_norm_mean', ...
    'dSC_Gyr_x_range', 'dSC_Gyr_y_range', 'dSC_Gyr_z_range', 'dSC_Gyr_norm_range', ...
    'dSC_Gyr_x_rms', 'dSC_Gyr_y_rms', 'dSC_Gyr_z_rms', 'dSC_Gyr_norm_rms', ...
    'dSC_Gyr_x_std', 'dSC_Gyr_y_std', 'dSC_Gyr_z_std', 'dSC_Gyr_norm_std', ...
    'dSC_Gyr_x_skew', 'dSC_Gyr_y_skew', 'dSC_Gyr_z_skew', 'dSC_Gyr_norm_skew', ...
    'dSC_Gyr_x_kurtosis', 'dSC_Gyr_y_kurtosis', 'dSC_Gyr_z_kurtosis', 'dSC_Gyr_norm_kurtosis', ...
    'SC_Gyr_corr_xy', 'SC_Gyr_corr_xz', 'SC_Gyr_corr_yz', ...
    'SC_Gyr_x_SamEn', 'SC_Gyr_y_SamEn', 'SC_Gyr_z_SamEn', 'SC_Gyr_norm_SamEn', ...
    'SC_Gyr_x_DAmp', 'SC_Gyr_x_DFreq', 'SC_Gyr_x_PSD_mean', 'SC_Gyr_x_PSD_std', 'SC_Gyr_x_PSD_skew', 'SC_Gyr_x_PSD_kurtosis', ...
    'SC_Gyr_y_DAmp', 'SC_Gyr_y_DFreq', 'SC_Gyr_y_PSD_mean', 'SC_Gyr_y_PSD_std', 'SC_Gyr_y_PSD_skew', 'SC_Gyr_y_PSD_kurtosis', ...
    'SC_Gyr_z_DAmp', 'SC_Gyr_z_DFreq', 'SC_Gyr_z_PSD_mean', 'SC_Gyr_z_PSD_std', 'SC_Gyr_z_PSD_skew', 'SC_Gyr_z_PSD_kurtosis', ...
    'SC_Gyr_norm_DAmp', 'SC_Gyr_norm_DFreq', 'SC_Gyr_norm_PSD_mean', 'SC_Gyr_norm_PSD_std', 'SC_Gyr_norm_PSD_skew', 'SC_Gyr_norm_PSD_kurtosis'};
BBS_table_wide = unstack(BBS_table,variableList,'trial_No')

% writetable(BBS_table,'General_Feature_Matrix_Admission_BBS_CVA_longformat_061720.csv','Delimiter',',','QuoteStrings',true)
% writetable(BBS_table_wide,'General_Feature_Matrix_Admission_CVA_061720.csv','Delimiter',',','QuoteStrings',true)

