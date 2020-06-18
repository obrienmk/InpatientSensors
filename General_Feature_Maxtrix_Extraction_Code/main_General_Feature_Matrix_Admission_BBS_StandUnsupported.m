%------------------------%
%
% BBS Feature Matrix for N02 Stand Unsupported with different time cutoffs
% Adam P. Horin
% 06/08/2020
%
%-------------------------%


clc
clear all
close all


ID = [1:55]; %enter ID or range of ID's to be processed CONTROLS [1:25,27:51]; CVA [1:55]
Type_of_Subject = 'CVA'; % enter CONTROLS or CVA


%Activity Selection
Activity = 'BBS';

Hz = 31.25;
dt = 1/Hz;
k = 1;
BBS_table = [];
TrialNames = {'N2_STAND_UNSUPPORTED_S1'};
TrialNumbers = {'N02'};
TrialTimes = [15, 30, 45, 60, 75, 90, 105, 120];

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
                SubjectID=['CVA0' str];
            else strcmp(Type_of_Subject,'CONTROLS')
                SubjectID=['HC0' str];
            end
        else
            if strcmp(Type_of_Subject,'CVA')
                SubjectID=['CVA' str];
            else strcmp(Type_of_Subject,'CONTROLS')
                SubjectID=['HC' str];
            end
        end
    end
    
    
    for g = 1:1:length(TrialTimes)
        
        %Trial N
        N_index = find(strcmp(data.Session_trials{SN}, TrialNames(1)));
        if isempty(N_index) == 1
            
            %Trial Information
            subject(g,:) = SubjectID;
            group(g,:) = Type_of_Subject;
            activity(g,:) = Activity;
            trial_No(g,:) = TrialNumbers(1);
            cutoff(g,:) = TrialTimes(g);
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
            subject(g,:) = SubjectID;
            group(g,:) = Type_of_Subject;
            activity(g,:) = Activity;
            trial_No(g,:) = TrialNumbers(1);
            cutoff(g,:) = TrialTimes(g);
            
            
            % Gyroscope data
            SC_Gyr = data.Session{SN}.Motion.SC.Gyr{N_index};
            
            % Acceleration data
            SC_Acc = data.Session{SN}.Motion.SC.Acc{N_index};
            
            % Demean Acc data
            SC_Acc = SC_Acc - ones(length(SC_Acc),1)*mean(SC_Acc);
            
            % To match final time
            
            %final = min([length(SC_Gyr) length(SC_Acc)]);
            final = round(Hz*TrialTimes(g));
            if final >= length(data.Session{SN}.Motion.Time{N_index}(1:end,:))
                Time = data.Session{SN}.Motion.Time{N_index}(1:end,:);
                SC_Gyr = SC_Gyr(1:end,:);
                SC_Acc = SC_Acc(1:end,:);
                trialtime(g,:) = data.Session{SN}.Motion.Time{N_index}(end,:)-data.Session{SN}.Motion.Time{N_index}(1,:);
            else
                Time = data.Session{SN}.Motion.Time{N_index}(1:final,:);
                SC_Gyr = SC_Gyr(1:end,:);
                SC_Acc = SC_Acc(1:end,:);
                trialtime(g,:) = data.Session{SN}.Motion.Time{N_index}(final,:)-data.Session{SN}.Motion.Time{N_index}(1,:);
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
            
        end
    end
    
    
    
    %BBS feature matrix in long format
    BBS_table_temp = table(subject, group, activity, trial_No, cutoff, trialtime,...
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
        SC_Gyr_norm_DAmp, SC_Gyr_norm_DFreq, SC_Gyr_norm_PSD_mean, SC_Gyr_norm_PSD_std, SC_Gyr_norm_PSD_skew, SC_Gyr_norm_PSD_kurtosis);
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
BBS_table_wide = unstack(BBS_table,variableList,'cutoff')

writetable(BBS_table,'General_Feature_Matrix_Admission_BBS_CVA_Times_longformat_061820.csv','Delimiter',',','QuoteStrings',true)
writetable(BBS_table_wide,'General_Feature_Matrix_Admission_BBS_CVA_Times_061820.csv','Delimiter',',','QuoteStrings',true)

