clc
clear all
close all


Sub_ID = [];
Sub_Type = [];
Cut_Off_Time = [];
FIM = [];
FIM_M = [];
BBS = [];
MWT10_SSV = [];
MWT10_FV = [];
MWT6 = [];
TUG = [];

AoM_Pel_tilt = [];
AoM_Pel_ro = [];
AoM_Pel_oblq = [];
AoM_Pel_norm = [];
AoM_Ankle_AS_x = [];  
AoM_Ankle_US_x = [];
AoM_Ankle_AS_y = []; 
AoM_Ankle_US_y = []; 
AoM_Ankle_AS_z = [];  
AoM_Ankle_US_z = [];  
AoM_Ankle_AS_norm = [];  
AoM_Ankle_US_norm = [];  

% General Features
% Gyro mean
SC_Gyr_x_mean = [];
SC_Gyr_y_mean = [];
SC_Gyr_z_mean = [];
SC_Gyr_norm_mean = [];
DLS_R_Gyr_x_mean = [];
DLS_R_Gyr_y_mean = [];
DLS_R_Gyr_z_mean = [];
DLS_R_Gyr_norm_mean = [];
DLS_L_Gyr_x_mean = [];
DLS_L_Gyr_y_mean = [];
DLS_L_Gyr_z_mean = [];
DLS_L_Gyr_norm_mean = [];

% Range
SC_Gyr_x_range = [];
SC_Gyr_y_range = [];
SC_Gyr_z_range = [];
SC_Gyr_norm_range = [];
DLS_R_Gyr_x_range = [];
DLS_R_Gyr_y_range = [];
DLS_R_Gyr_z_range = [];
DLS_R_Gyr_norm_range = [];
DLS_L_Gyr_x_range = [];
DLS_L_Gyr_y_range = [];
DLS_L_Gyr_z_range = [];
DLS_L_Gyr_norm_range = [];

% RMS
SC_Gyr_x_rms = [];
SC_Gyr_y_rms = [];
SC_Gyr_z_rms = [];
SC_Gyr_norm_rms = [];
DLS_R_Gyr_x_rms = [];
DLS_R_Gyr_y_rms = [];
DLS_R_Gyr_z_rms = [];
DLS_R_Gyr_norm_rms = [];
DLS_L_Gyr_x_rms = [];
DLS_L_Gyr_y_rms = [];
DLS_L_Gyr_z_rms = [];
DLS_L_Gyr_norm_rms = [];

% Standard Deviation
SC_Gyr_x_std = [];
SC_Gyr_y_std = [];
SC_Gyr_z_std = [];
SC_Gyr_norm_std = [];
DLS_R_Gyr_x_std = [];
DLS_R_Gyr_y_std = [];
DLS_R_Gyr_z_std = [];
DLS_R_Gyr_norm_std = [];
DLS_L_Gyr_x_std = [];
DLS_L_Gyr_y_std = [];
DLS_L_Gyr_z_std = [];
DLS_L_Gyr_norm_std = [];

% Skew
SC_Gyr_x_skew = [];
SC_Gyr_y_skew = [];
SC_Gyr_z_skew = [];
SC_Gyr_norm_skew = [];
DLS_R_Gyr_x_skew = [];
DLS_R_Gyr_y_skew = [];
DLS_R_Gyr_z_skew = [];
DLS_R_Gyr_norm_skew = [];
DLS_L_Gyr_x_skew = [];
DLS_L_Gyr_y_skew = [];
DLS_L_Gyr_z_skew = [];
DLS_L_Gyr_norm_skew = [];

% Kurtosis
SC_Gyr_x_kurtosis = [];
SC_Gyr_y_kurtosis = [];
SC_Gyr_z_kurtosis = [];
SC_Gyr_norm_kurtosis = [];
DLS_R_Gyr_x_kurtosis = [];
DLS_R_Gyr_y_kurtosis = [];
DLS_R_Gyr_z_kurtosis = [];
DLS_R_Gyr_norm_kurtosis = [];
DLS_L_Gyr_x_kurtosis = [];
DLS_L_Gyr_y_kurtosis = [];
DLS_L_Gyr_z_kurtosis = [];
DLS_L_Gyr_norm_kurtosis = [];

% Pearson correlation coefficient
SC_Gyr_corr_xy = [];
SC_Gyr_corr_xz = [];
SC_Gyr_corr_yz = [];
DLS_R_Gyr_corr_xy = [];
DLS_R_Gyr_corr_xz = [];
DLS_R_Gyr_corr_yz = [];
DLS_L_Gyr_corr_xy = [];
DLS_L_Gyr_corr_xz = [];
DLS_L_Gyr_corr_yz = [];

% Sample Entropy
SC_Gyr_x_SamEn = [];
SC_Gyr_y_SamEn = [];
SC_Gyr_z_SamEn = [];
SC_Gyr_norm_SamEn = [];
DLS_R_Gyr_x_SamEn = [];
DLS_R_Gyr_y_SamEn = [];
DLS_R_Gyr_z_SamEn = [];
DLS_R_Gyr_norm_SamEn = [];
DLS_L_Gyr_x_SamEn = [];
DLS_L_Gyr_y_SamEn = [];
DLS_L_Gyr_z_SamEn = [];
DLS_L_Gyr_norm_SamEn = [];

% Frequency domain
SC_Gyr_x_DAmp = [];
SC_Gyr_x_DFreq = [];
SC_Gyr_x_PSD_mean = [];
SC_Gyr_x_PSD_std = [];
SC_Gyr_x_PSD_skew = [];
SC_Gyr_x_PSD_kurtosis = [];

SC_Gyr_y_DAmp = [];
SC_Gyr_y_DFreq = [];
SC_Gyr_y_PSD_mean = [];
SC_Gyr_y_PSD_std = [];
SC_Gyr_y_PSD_skew = [];
SC_Gyr_y_PSD_kurtosis = [];

SC_Gyr_z_DAmp = [];
SC_Gyr_z_DFreq = [];
SC_Gyr_z_PSD_mean = [];
SC_Gyr_z_PSD_std = [];
SC_Gyr_z_PSD_skew = [];
SC_Gyr_z_PSD_kurtosis = [];

SC_Gyr_norm_DAmp = [];
SC_Gyr_norm_DFreq = [];
SC_Gyr_norm_PSD_mean = [];
SC_Gyr_norm_PSD_std = [];
SC_Gyr_norm_PSD_skew = [];
SC_Gyr_norm_PSD_kurtosis = [];

DLS_R_Gyr_x_DAmp = [];
DLS_R_Gyr_x_DFreq = [];
DLS_R_Gyr_x_PSD_mean = [];
DLS_R_Gyr_x_PSD_std = [];
DLS_R_Gyr_x_PSD_skew = [];
DLS_R_Gyr_x_PSD_kurtosis = [];

DLS_R_Gyr_y_DAmp = [];
DLS_R_Gyr_y_DFreq = [];
DLS_R_Gyr_y_PSD_mean = [];
DLS_R_Gyr_y_PSD_std = [];
DLS_R_Gyr_y_PSD_skew = [];
DLS_R_Gyr_y_PSD_kurtosis = [];

DLS_R_Gyr_z_DAmp = [];
DLS_R_Gyr_z_DFreq = [];
DLS_R_Gyr_z_PSD_mean = [];
DLS_R_Gyr_z_PSD_std = [];
DLS_R_Gyr_z_PSD_skew = [];
DLS_R_Gyr_z_PSD_kurtosis = [];

DLS_R_Gyr_norm_DAmp = [];
DLS_R_Gyr_norm_DFreq = [];
DLS_R_Gyr_norm_PSD_mean = [];
DLS_R_Gyr_norm_PSD_std = [];
DLS_R_Gyr_norm_PSD_skew = [];
DLS_R_Gyr_norm_PSD_kurtosis = [];

DLS_L_Gyr_x_DAmp = [];
DLS_L_Gyr_x_DFreq = [];
DLS_L_Gyr_x_PSD_mean = [];
DLS_L_Gyr_x_PSD_std = [];
DLS_L_Gyr_x_PSD_skew = [];
DLS_L_Gyr_x_PSD_kurtosis = [];

DLS_L_Gyr_y_DAmp = [];
DLS_L_Gyr_y_DFreq = [];
DLS_L_Gyr_y_PSD_mean = [];
DLS_L_Gyr_y_PSD_std = [];
DLS_L_Gyr_y_PSD_skew = [];
DLS_L_Gyr_y_PSD_kurtosis = [];

DLS_L_Gyr_z_DAmp = [];
DLS_L_Gyr_z_DFreq = [];
DLS_L_Gyr_z_PSD_mean = [];
DLS_L_Gyr_z_PSD_std = [];
DLS_L_Gyr_z_PSD_skew = [];
DLS_L_Gyr_z_PSD_kurtosis = [];

DLS_L_Gyr_norm_DAmp = [];
DLS_L_Gyr_norm_DFreq = [];
DLS_L_Gyr_norm_PSD_mean = [];
DLS_L_Gyr_norm_PSD_std = [];
DLS_L_Gyr_norm_PSD_skew = [];
DLS_L_Gyr_norm_PSD_kurtosis = [];


% Acc mean
% General Features
% Gyro mean
SC_Acc_x_mean = [];
SC_Acc_y_mean = [];
SC_Acc_z_mean = [];
SC_Acc_norm_mean = [];
DLS_R_Acc_x_mean = [];
DLS_R_Acc_y_mean = [];
DLS_R_Acc_z_mean = [];
DLS_R_Acc_norm_mean = [];
DLS_L_Acc_x_mean = [];
DLS_L_Acc_y_mean = [];
DLS_L_Acc_z_mean = [];
DLS_L_Acc_norm_mean = [];

% Range
SC_Acc_x_range = [];
SC_Acc_y_range = [];
SC_Acc_z_range = [];
SC_Acc_norm_range = [];
DLS_R_Acc_x_range = [];
DLS_R_Acc_y_range = [];
DLS_R_Acc_z_range = [];
DLS_R_Acc_norm_range = [];
DLS_L_Acc_x_range = [];
DLS_L_Acc_y_range = [];
DLS_L_Acc_z_range = [];
DLS_L_Acc_norm_range = [];

% RMS
SC_Acc_x_rms = [];
SC_Acc_y_rms = [];
SC_Acc_z_rms = [];
SC_Acc_norm_rms = [];
DLS_R_Acc_x_rms = [];
DLS_R_Acc_y_rms = [];
DLS_R_Acc_z_rms = [];
DLS_R_Acc_norm_rms = [];
DLS_L_Acc_x_rms = [];
DLS_L_Acc_y_rms = [];
DLS_L_Acc_z_rms = [];
DLS_L_Acc_norm_rms = [];

% Standard Deviation
SC_Acc_x_std = [];
SC_Acc_y_std = [];
SC_Acc_z_std = [];
SC_Acc_norm_std = [];
DLS_R_Acc_x_std = [];
DLS_R_Acc_y_std = [];
DLS_R_Acc_z_std = [];
DLS_R_Acc_norm_std = [];
DLS_L_Acc_x_std = [];
DLS_L_Acc_y_std = [];
DLS_L_Acc_z_std = [];
DLS_L_Acc_norm_std = [];

% Skew
SC_Acc_x_skew = [];
SC_Acc_y_skew = [];
SC_Acc_z_skew = [];
SC_Acc_norm_skew = [];
DLS_R_Acc_x_skew = [];
DLS_R_Acc_y_skew = [];
DLS_R_Acc_z_skew = [];
DLS_R_Acc_norm_skew = [];
DLS_L_Acc_x_skew = [];
DLS_L_Acc_y_skew = [];
DLS_L_Acc_z_skew = [];
DLS_L_Acc_norm_skew = [];

% Kurtosis
SC_Acc_x_kurtosis = [];
SC_Acc_y_kurtosis = [];
SC_Acc_z_kurtosis = [];
SC_Acc_norm_kurtosis = [];
DLS_R_Acc_x_kurtosis = [];
DLS_R_Acc_y_kurtosis = [];
DLS_R_Acc_z_kurtosis = [];
DLS_R_Acc_norm_kurtosis = [];
DLS_L_Acc_x_kurtosis = [];
DLS_L_Acc_y_kurtosis = [];
DLS_L_Acc_z_kurtosis = [];
DLS_L_Acc_norm_kurtosis = [];

% Pearson correlation coefficient
SC_Acc_corr_xy = [];
SC_Acc_corr_xz = [];
SC_Acc_corr_yz = [];
DLS_R_Acc_corr_xy = [];
DLS_R_Acc_corr_xz = [];
DLS_R_Acc_corr_yz = [];
DLS_L_Acc_corr_xy = [];
DLS_L_Acc_corr_xz = [];
DLS_L_Acc_corr_yz = [];

% Sample Entropy
SC_Acc_x_SamEn = [];
SC_Acc_y_SamEn = [];
SC_Acc_z_SamEn = [];
SC_Acc_norm_SamEn = [];
DLS_R_Acc_x_SamEn = [];
DLS_R_Acc_y_SamEn = [];
DLS_R_Acc_z_SamEn = [];
DLS_R_Acc_norm_SamEn = [];
DLS_L_Acc_x_SamEn = [];
DLS_L_Acc_y_SamEn = [];
DLS_L_Acc_z_SamEn = [];
DLS_L_Acc_norm_SamEn = [];

% Frequency domain
SC_Acc_x_DAmp = [];
SC_Acc_x_DFreq = [];
SC_Acc_x_PSD_mean = [];
SC_Acc_x_PSD_std = [];
SC_Acc_x_PSD_skew = [];
SC_Acc_x_PSD_kurtosis = [];

SC_Acc_y_DAmp = [];
SC_Acc_y_DFreq = [];
SC_Acc_y_PSD_mean = [];
SC_Acc_y_PSD_std = [];
SC_Acc_y_PSD_skew = [];
SC_Acc_y_PSD_kurtosis = [];

SC_Acc_z_DAmp = [];
SC_Acc_z_DFreq = [];
SC_Acc_z_PSD_mean = [];
SC_Acc_z_PSD_std = [];
SC_Acc_z_PSD_skew = [];
SC_Acc_z_PSD_kurtosis = [];

SC_Acc_norm_DAmp = [];
SC_Acc_norm_DFreq = [];
SC_Acc_norm_PSD_mean = [];
SC_Acc_norm_PSD_std = [];
SC_Acc_norm_PSD_skew = [];
SC_Acc_norm_PSD_kurtosis = [];

DLS_R_Acc_x_DAmp = [];
DLS_R_Acc_x_DFreq = [];
DLS_R_Acc_x_PSD_mean = [];
DLS_R_Acc_x_PSD_std = [];
DLS_R_Acc_x_PSD_skew = [];
DLS_R_Acc_x_PSD_kurtosis = [];

DLS_R_Acc_y_DAmp = [];
DLS_R_Acc_y_DFreq = [];
DLS_R_Acc_y_PSD_mean = [];
DLS_R_Acc_y_PSD_std = [];
DLS_R_Acc_y_PSD_skew = [];
DLS_R_Acc_y_PSD_kurtosis = [];

DLS_R_Acc_z_DAmp = [];
DLS_R_Acc_z_DFreq = [];
DLS_R_Acc_z_PSD_mean = [];
DLS_R_Acc_z_PSD_std = [];
DLS_R_Acc_z_PSD_skew = [];
DLS_R_Acc_z_PSD_kurtosis = [];

DLS_R_Acc_norm_DAmp = [];
DLS_R_Acc_norm_DFreq = [];
DLS_R_Acc_norm_PSD_mean = [];
DLS_R_Acc_norm_PSD_std = [];
DLS_R_Acc_norm_PSD_skew = [];
DLS_R_Acc_norm_PSD_kurtosis = [];

DLS_L_Acc_x_DAmp = [];
DLS_L_Acc_x_DFreq = [];
DLS_L_Acc_x_PSD_mean = [];
DLS_L_Acc_x_PSD_std = [];
DLS_L_Acc_x_PSD_skew = [];
DLS_L_Acc_x_PSD_kurtosis = [];

DLS_L_Acc_y_DAmp = [];
DLS_L_Acc_y_DFreq = [];
DLS_L_Acc_y_PSD_mean = [];
DLS_L_Acc_y_PSD_std = [];
DLS_L_Acc_y_PSD_skew = [];
DLS_L_Acc_y_PSD_kurtosis = [];

DLS_L_Acc_z_DAmp = [];
DLS_L_Acc_z_DFreq = [];
DLS_L_Acc_z_PSD_mean = [];
DLS_L_Acc_z_PSD_std = [];
DLS_L_Acc_z_PSD_skew = [];
DLS_L_Acc_z_PSD_kurtosis = [];

DLS_L_Acc_norm_DAmp = [];
DLS_L_Acc_norm_DFreq = [];
DLS_L_Acc_norm_PSD_mean = [];
DLS_L_Acc_norm_PSD_std = [];
DLS_L_Acc_norm_PSD_skew = [];
DLS_L_Acc_norm_PSD_kurtosis = [];


Hz = 31.25;
dt = 1/Hz;
TS = [5 10 20 30 60 90 120 180 240 300 360]

for k = 1:1:length(TS)
    
    % CVA Admission
    % ID: 1, 25  --> no data
    % ID: 7, 40 --> sensor at DLS_L missing
    % ID: 8 --> sensor at DLS_R missing
    % ID: 32 --> long data 3 times*6MWT? 
    % ID: 3, 22, 30 42 44 49 --> No admission sensor data 
    % ID: 35, 43 --> Short data 5595 (179sec, 43), 10167 (325sec, 35)  frames 
    % ID: 26, 54 --> Outlier. Sensor data broken --> fixed code
    
    ID_AD = [2 4:6 9:21 23:24 26:29 31 33:34 36:39 41 45:48 50:55]
    Type_of_Subject = 'CVA'

    for n = 1:1:length(ID_AD)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID_AD(n)) '.mat']
        load(file_input);

        Side = get_side('./CVA_Paretic_Side.csv');
        AS = Side.Side(ID_AD(n));

        SN = 1;     % Session #
        TN = 1;     % Trial #

        if ID_AD(n) == 26 || ID_AD(n) == 54
            SC_Gyr = [data.Session{SN}.Motion.SC.Gyr{1}; data.Session{SN}.Motion.SC.Gyr{2}];
            DLS_R_Gyr = [data.Session{SN}.Motion.DLS_R.Gyr{1}; data.Session{SN}.Motion.DLS_R.Gyr{2}];
            DLS_L_Gyr = [data.Session{SN}.Motion.DLS_L.Gyr{1}; data.Session{SN}.Motion.DLS_L.Gyr{2}];

            % Acceleration data
            SC_Acc = [data.Session{SN}.Motion.SC.Acc{1}; data.Session{SN}.Motion.SC.Acc{2}];
            DLS_R_Acc = [data.Session{SN}.Motion.DLS_R.Acc{1}; data.Session{SN}.Motion.DLS_R.Acc{2}];
            DLS_L_Acc = [data.Session{SN}.Motion.DLS_L.Acc{1}; data.Session{SN}.Motion.DLS_L.Acc{2}];

             % To match final time
%              final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
            final = Hz*TS(k);  
            Time_ = [data.Session{SN}.Motion.Time{1}; data.Session{SN}.Motion.Time{2}];
            Time = Time_(1:final,:);
            SC_Gyr = SC_Gyr(1:final,:);
            DLS_R_Gyr = DLS_R_Gyr(1:final,:);
            DLS_L_Gyr = DLS_L_Gyr(1:final,:);
        else

            % Gyroscope data
            SC_Gyr = data.Session{SN}.Motion.SC.Gyr{TN};
            DLS_R_Gyr = data.Session{SN}.Motion.DLS_R.Gyr{TN};
            DLS_L_Gyr = data.Session{SN}.Motion.DLS_L.Gyr{TN};

            % Acceleration data
            SC_Acc = data.Session{SN}.Motion.SC.Acc{TN};
            DLS_R_Acc = data.Session{SN}.Motion.DLS_R.Acc{TN};
            DLS_L_Acc = data.Session{SN}.Motion.DLS_L.Acc{TN};

             % To match final time
%              final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
            final = Hz*TS(k);  

            Time = data.Session{SN}.Motion.Time{TN}(1:final,:);
            SC_Gyr = SC_Gyr(1:final,:);
            DLS_R_Gyr = DLS_R_Gyr(1:final,:);
            DLS_L_Gyr = DLS_L_Gyr(1:final,:);
        end
        
        % Norm of Gyro & Acc
        for i = 1:1:length(Time)
            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
            DLS_R_Gyr_norm(i,:) = norm(DLS_R_Gyr(i,:));
            DLS_L_Gyr_norm(i,:) = norm(DLS_L_Gyr(i,:));

            SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
            DLS_R_Acc_norm(i,:) = norm(DLS_R_Acc(i,:));
            DLS_L_Acc_norm(i,:) = norm(DLS_L_Acc(i,:));
        end 
        
        % Step Count
        f = Hz;
        f_cutoff = 2;
        f_norm = f_cutoff/(f/2);
        ftype = 'low';
        N = 4;  % Order of low pass Butterworth filter
        [B, A] = butter(N,f_norm,ftype);
        DLS_R_Gyr_filt = filtfilt(B, A, DLS_R_Gyr_norm);
        DLS_L_Gyr_filt = filtfilt(B, A, DLS_L_Gyr_norm);
        
        for i = 1:1:length(DLS_R_Gyr_filt)
            if DLS_R_Gyr_filt(i) <  20
                DLS_R_Gyr_filt(i) = 0;
            else
                DLS_R_Gyr_filt(i) = DLS_R_Gyr_filt(i);
            end

            if DLS_L_Gyr_filt(i) <  20
                DLS_L_Gyr_filt(i) = 0;
            else
                DLS_L_Gyr_filt(i) = DLS_L_Gyr_filt(i);
            end
        end

        [pks_G_R, locs_G_R] = findpeaks(DLS_R_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_R_Gyr_filt)+std(DLS_R_Gyr_filt));
        [pks_G_L, locs_G_L] = findpeaks(DLS_L_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_L_Gyr_filt)+std(DLS_L_Gyr_filt));

        Step_R = length(locs_G_R);
        Step_L = length(locs_G_L);
        
        AD.Steps(n,:) = Step_R + Step_L;

        % Amount of motion
        % Pelvic Amount of motion 
        AD.AoM_Pel_tilt(n,:) = sum(abs(SC_Gyr(:,1))) * dt;
        AD.AoM_Pel_ro(n,:) = sum(abs(SC_Gyr(:,2))) * dt;
        AD.AoM_Pel_oblq(n,:) = sum(abs(SC_Gyr(:,3))) * dt;
        AD.AoM_Pel_norm(n,:) = sum(abs(SC_Gyr_norm)) * dt; 
        
        % Ankle Amount of motion
        if strcmp(AS,'L') == 1
            AD.AoM_Ankle_AS_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;  
            AD.AoM_Ankle_US_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;
            AD.AoM_Ankle_AS_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;  
            AD.AoM_Ankle_US_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;
            AD.AoM_Ankle_AS_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;  
            AD.AoM_Ankle_US_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;
            AD.AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt; 
            AD.AoM_Ankle_US_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt;

        elseif strcmp(AS,'R') == 1
            AD.AoM_Ankle_AS_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;  
            AD.AoM_Ankle_US_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;
            AD.AoM_Ankle_AS_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;  
            AD.AoM_Ankle_US_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;
            AD.AoM_Ankle_AS_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;  
            AD.AoM_Ankle_US_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;
            AD.AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt; 
            AD.AoM_Ankle_US_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt;
        end  

        % General Features
        % Gyro mean
        AD.SC_Gyr_x_mean(n,:) = mean(SC_Gyr(:,1));
        AD.SC_Gyr_y_mean(n,:) = mean(SC_Gyr(:,2));
        AD.SC_Gyr_z_mean(n,:) = mean(SC_Gyr(:,3));
        AD.SC_Gyr_norm_mean(n,:) = mean(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_mean(n,:) = mean(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_mean(n,:) = mean(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_mean(n,:) = mean(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_mean(n,:) = mean(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_mean(n,:) = mean(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_mean(n,:) = mean(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_mean(n,:) = mean(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_mean(n,:) = mean(DLS_L_Gyr_norm);

        % Range
        AD.SC_Gyr_x_range(n,:) = range(SC_Gyr(:,1));
        AD.SC_Gyr_y_range(n,:) = range(SC_Gyr(:,2));
        AD.SC_Gyr_z_range(n,:) = range(SC_Gyr(:,3));
        AD.SC_Gyr_norm_range(n,:) = range(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_range(n,:) = range(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_range(n,:) = range(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_range(n,:) = range(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_range(n,:) = range(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_range(n,:) = range(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_range(n,:) = range(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_range(n,:) = range(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_range(n,:) = range(DLS_L_Gyr_norm);

        % RMS
        AD.SC_Gyr_x_rms(n,:) = rms(SC_Gyr(:,1));
        AD.SC_Gyr_y_rms(n,:) = rms(SC_Gyr(:,2));
        AD.SC_Gyr_z_rms(n,:) = rms(SC_Gyr(:,3));
        AD.SC_Gyr_norm_rms(n,:) = rms(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_rms(n,:) = rms(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_rms(n,:) = rms(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_rms(n,:) = rms(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_rms(n,:) = rms(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_rms(n,:) = rms(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_rms(n,:) = rms(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_rms(n,:) = rms(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_rms(n,:) = rms(DLS_L_Gyr_norm);

        % Standard Deviation
        AD.SC_Gyr_x_std(n,:) = std(SC_Gyr(:,1));
        AD.SC_Gyr_y_std(n,:) = std(SC_Gyr(:,2));
        AD.SC_Gyr_z_std(n,:) = std(SC_Gyr(:,3));
        AD.SC_Gyr_norm_std(n,:) = std(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_std(n,:) = std(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_std(n,:) = std(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_std(n,:) = std(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_std(n,:) = std(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_std(n,:) = std(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_std(n,:) = std(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_std(n,:) = std(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_std(n,:) = std(DLS_L_Gyr_norm);

        % Skew
        AD.SC_Gyr_x_skew(n,:) = skewness(SC_Gyr(:,1));
        AD.SC_Gyr_y_skew(n,:) = skewness(SC_Gyr(:,2));
        AD.SC_Gyr_z_skew(n,:) = skewness(SC_Gyr(:,3));
        AD.SC_Gyr_norm_skew(n,:) = skewness(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_skew(n,:) = skewness(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_skew(n,:) = skewness(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_skew(n,:) = skewness(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_skew(n,:) = skewness(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_skew(n,:) = skewness(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_skew(n,:) = skewness(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_skew(n,:) = skewness(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_skew(n,:) = skewness(DLS_L_Gyr_norm);

        % Kurtosis
        AD.SC_Gyr_x_kurtosis(n,:) = kurtosis(SC_Gyr(:,1));
        AD.SC_Gyr_y_kurtosis(n,:) = kurtosis(SC_Gyr(:,2));
        AD.SC_Gyr_z_kurtosis(n,:) = kurtosis(SC_Gyr(:,3));
        AD.SC_Gyr_norm_kurtosis(n,:) = kurtosis(SC_Gyr_norm);
        AD.DLS_R_Gyr_x_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,1));
        AD.DLS_R_Gyr_y_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_z_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_R_Gyr_norm);
        AD.DLS_L_Gyr_x_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,1));
        AD.DLS_L_Gyr_y_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_z_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_L_Gyr_norm);

        % Pearson correlation coefficient
        AD.SC_Gyr_corr_xy(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,2));
        AD.SC_Gyr_corr_xz(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,3));
        AD.SC_Gyr_corr_yz(n,:) = corr(SC_Gyr(:,2),SC_Gyr(:,3));
        AD.DLS_R_Gyr_corr_xy(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,2));
        AD.DLS_R_Gyr_corr_xz(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,3));
        AD.DLS_R_Gyr_corr_yz(n,:) = corr(DLS_R_Gyr(:,2),DLS_R_Gyr(:,3));
        AD.DLS_L_Gyr_corr_xy(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,2));
        AD.DLS_L_Gyr_corr_xz(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,3));
        AD.DLS_L_Gyr_corr_yz(n,:) = corr(DLS_L_Gyr(:,2),DLS_L_Gyr(:,3));

        % Sample Entropy
        r = 0.2;
        AD.SC_Gyr_x_SamEn(n,:) = sampen(SC_Gyr(:,1),1,r);
        AD.SC_Gyr_y_SamEn(n,:) = sampen(SC_Gyr(:,2),1,r);
        AD.SC_Gyr_z_SamEn(n,:) = sampen(SC_Gyr(:,3),1,r);
        AD.SC_Gyr_norm_SamEn(n,:) = sampen(SC_Gyr_norm,1,r);
        AD.DLS_R_Gyr_x_SamEn(n,:) = sampen(DLS_R_Gyr(:,1),1,r);
        AD.DLS_R_Gyr_y_SamEn(n,:) = sampen(DLS_R_Gyr(:,2),1,r);
        AD.DLS_R_Gyr_z_SamEn(n,:) = sampen(DLS_R_Gyr(:,3),1,r);
        AD.DLS_R_Gyr_norm_SamEn(n,:) = sampen(DLS_R_Gyr_norm,1,r);
        AD.DLS_L_Gyr_x_SamEn(n,:) = sampen(DLS_L_Gyr(:,1),1,r);
        AD.DLS_L_Gyr_y_SamEn(n,:) = sampen(DLS_L_Gyr(:,2),1,r);
        AD.DLS_L_Gyr_z_SamEn(n,:) = sampen(DLS_L_Gyr(:,3),1,r);
        AD.DLS_L_Gyr_norm_SamEn(n,:) = sampen(DLS_L_Gyr_norm,1,r);

        % Frequency Domain
        ff = FFeatures(SC_Gyr(:,1), Hz);
        AD.SC_Gyr_x_DAmp(n,:) = ff(1);
        AD.SC_Gyr_x_DFreq(n,:) = ff(2);
        AD.SC_Gyr_x_PSD_mean(n,:) = ff(3);
        AD.SC_Gyr_x_PSD_std(n,:) = ff(4);
        AD.SC_Gyr_x_PSD_skew(n,:) = ff(5);
        AD.SC_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,2), Hz);
        AD.SC_Gyr_y_DAmp(n,:) = ff(1);
        AD.SC_Gyr_y_DFreq(n,:) = ff(2);
        AD.SC_Gyr_y_PSD_mean(n,:) = ff(3);
        AD.SC_Gyr_y_PSD_std(n,:) = ff(4);
        AD.SC_Gyr_y_PSD_skew(n,:) = ff(5);
        AD.SC_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,3), Hz);
        AD.SC_Gyr_z_DAmp(n,:) = ff(1);
        AD.SC_Gyr_z_DFreq(n,:) = ff(2);
        AD.SC_Gyr_z_PSD_mean(n,:) = ff(3);
        AD.SC_Gyr_z_PSD_std(n,:) = ff(4);
        AD.SC_Gyr_z_PSD_skew(n,:) = ff(5);
        AD.SC_Gyr_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Gyr_norm, Hz);
        AD.SC_Gyr_norm_DAmp(n,:) = ff(1);
        AD.SC_Gyr_norm_DFreq(n,:) = ff(2);
        AD.SC_Gyr_norm_PSD_mean(n,:) = ff(3);
        AD.SC_Gyr_norm_PSD_std(n,:) = ff(4);
        AD.SC_Gyr_norm_PSD_skew(n,:) = ff(5);
        AD.SC_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,1), Hz);
        AD.DLS_R_Gyr_x_DAmp(n,:) = ff(1);
        AD.DLS_R_Gyr_x_DFreq(n,:) = ff(2);
        AD.DLS_R_Gyr_x_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Gyr_x_PSD_std(n,:) = ff(4);
        AD.DLS_R_Gyr_x_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,2), Hz);
        AD.DLS_R_Gyr_y_DAmp(n,:) = ff(1);
        AD.DLS_R_Gyr_y_DFreq(n,:) = ff(2);
        AD.DLS_R_Gyr_y_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Gyr_y_PSD_std(n,:) = ff(4);
        AD.DLS_R_Gyr_y_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,3), Hz);
        AD.DLS_R_Gyr_z_DAmp(n,:) = ff(1);
        AD.DLS_R_Gyr_z_DFreq(n,:) = ff(2);
        AD.DLS_R_Gyr_z_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Gyr_z_PSD_std(n,:) = ff(4);
        AD.DLS_R_Gyr_z_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr_norm, Hz);
        AD.DLS_R_Gyr_norm_DAmp(n,:) = ff(1);
        AD.DLS_R_Gyr_norm_DFreq(n,:) = ff(2);
        AD.DLS_R_Gyr_norm_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Gyr_norm_PSD_std(n,:) = ff(4);
        AD.DLS_R_Gyr_norm_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,1), Hz);
        AD.DLS_L_Gyr_x_DAmp(n,:) = ff(1);
        AD.DLS_L_Gyr_x_DFreq(n,:) = ff(2);
        AD.DLS_L_Gyr_x_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Gyr_x_PSD_std(n,:) = ff(4);
        AD.DLS_L_Gyr_x_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,2), Hz);
        AD.DLS_L_Gyr_y_DAmp(n,:) = ff(1);
        AD.DLS_L_Gyr_y_DFreq(n,:) = ff(2);
        AD.DLS_L_Gyr_y_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Gyr_y_PSD_std(n,:) = ff(4);
        AD.DLS_L_Gyr_y_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,3), Hz);
        AD.DLS_L_Gyr_z_DAmp(n,:) = ff(1);
        AD.DLS_L_Gyr_z_DFreq(n,:) = ff(2);
        AD.DLS_L_Gyr_z_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Gyr_z_PSD_std(n,:) = ff(4);
        AD.DLS_L_Gyr_z_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr_norm, Hz);
        AD.DLS_L_Gyr_norm_DAmp(n,:) = ff(1);
        AD.DLS_L_Gyr_norm_DFreq(n,:) = ff(2);
        AD.DLS_L_Gyr_norm_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Gyr_norm_PSD_std(n,:) = ff(4);
        AD.DLS_L_Gyr_norm_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Gyr_norm_PSD_kurtosis(n,:) = ff(6);


        % Acc mean
        AD.SC_Acc_x_mean(n,:) = mean(SC_Acc(:,1));
        AD.SC_Acc_y_mean(n,:) = mean(SC_Acc(:,2));
        AD.SC_Acc_z_mean(n,:) = mean(SC_Acc(:,3));
        AD.SC_Acc_norm_mean(n,:) = mean(SC_Acc_norm);
        AD.DLS_R_Acc_x_mean(n,:) = mean(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_mean(n,:) = mean(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_mean(n,:) = mean(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_mean(n,:) = mean(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_mean(n,:) = mean(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_mean(n,:) = mean(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_mean(n,:) = mean(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_mean(n,:) = mean(DLS_L_Acc_norm);

        % Range
        AD.SC_Acc_x_range(n,:) = range(SC_Acc(:,1));
        AD.SC_Acc_y_range(n,:) = range(SC_Acc(:,2));
        AD.SC_Acc_z_range(n,:) = range(SC_Acc(:,3));
        AD.SC_Acc_norm_range(n,:) = range(SC_Acc_norm);
        AD.DLS_R_Acc_x_range(n,:) = range(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_range(n,:) = range(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_range(n,:) = range(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_range(n,:) = range(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_range(n,:) = range(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_range(n,:) = range(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_range(n,:) = range(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_range(n,:) = range(DLS_L_Acc_norm);

        % RMS
        AD.SC_Acc_x_rms(n,:) = rms(SC_Acc(:,1));
        AD.SC_Acc_y_rms(n,:) = rms(SC_Acc(:,2));
        AD.SC_Acc_z_rms(n,:) = rms(SC_Acc(:,3));
        AD.SC_Acc_norm_rms(n,:) = rms(SC_Acc_norm);
        AD.DLS_R_Acc_x_rms(n,:) = rms(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_rms(n,:) = rms(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_rms(n,:) = rms(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_rms(n,:) = rms(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_rms(n,:) = rms(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_rms(n,:) = rms(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_rms(n,:) = rms(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_rms(n,:) = rms(DLS_L_Acc_norm);

        % Standard Deviation
        AD.SC_Acc_x_std(n,:) = std(SC_Acc(:,1));
        AD.SC_Acc_y_std(n,:) = std(SC_Acc(:,2));
        AD.SC_Acc_z_std(n,:) = std(SC_Acc(:,3));
        AD.SC_Acc_norm_std(n,:) = std(SC_Acc_norm);
        AD.DLS_R_Acc_x_std(n,:) = std(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_std(n,:) = std(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_std(n,:) = std(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_std(n,:) = std(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_std(n,:) = std(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_std(n,:) = std(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_std(n,:) = std(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_std(n,:) = std(DLS_L_Acc_norm);

        % Skew
        AD.SC_Acc_x_skew(n,:) = skewness(SC_Acc(:,1));
        AD.SC_Acc_y_skew(n,:) = skewness(SC_Acc(:,2));
        AD.SC_Acc_z_skew(n,:) = skewness(SC_Acc(:,3));
        AD.SC_Acc_norm_skew(n,:) = skewness(SC_Acc_norm);
        AD.DLS_R_Acc_x_skew(n,:) = skewness(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_skew(n,:) = skewness(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_skew(n,:) = skewness(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_skew(n,:) = skewness(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_skew(n,:) = skewness(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_skew(n,:) = skewness(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_skew(n,:) = skewness(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_skew(n,:) = skewness(DLS_L_Acc_norm);

        % Kurtosis
        AD.SC_Acc_x_kurtosis(n,:) = kurtosis(SC_Acc(:,1));
        AD.SC_Acc_y_kurtosis(n,:) = kurtosis(SC_Acc(:,2));
        AD.SC_Acc_z_kurtosis(n,:) = kurtosis(SC_Acc(:,3));
        AD.SC_Acc_norm_kurtosis(n,:) = kurtosis(SC_Acc_norm);
        AD.DLS_R_Acc_x_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,1));
        AD.DLS_R_Acc_y_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,2));
        AD.DLS_R_Acc_z_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,3));
        AD.DLS_R_Acc_norm_kurtosis(n,:) = kurtosis(DLS_R_Acc_norm);
        AD.DLS_L_Acc_x_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,1));
        AD.DLS_L_Acc_y_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,2));
        AD.DLS_L_Acc_z_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,3));
        AD.DLS_L_Acc_norm_kurtosis(n,:) = kurtosis(DLS_L_Acc_norm);

        % Pearson correlation coefficient
        AD.SC_Acc_corr_xy(n,:) = corr(SC_Acc(:,1),SC_Acc(:,2));
        AD.SC_Acc_corr_xz(n,:) = corr(SC_Acc(:,1),SC_Acc(:,3));
        AD.SC_Acc_corr_yz(n,:) = corr(SC_Acc(:,2),SC_Acc(:,3));
        AD.DLS_R_Acc_corr_xy(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,2));
        AD.DLS_R_Acc_corr_xz(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,3));
        AD.DLS_R_Acc_corr_yz(n,:) = corr(DLS_R_Acc(:,2),DLS_R_Acc(:,3));
        AD.DLS_L_Acc_corr_xy(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,2));
        AD.DLS_L_Acc_corr_xz(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,3));
        AD.DLS_L_Acc_corr_yz(n,:) = corr(DLS_L_Acc(:,2),DLS_L_Acc(:,3));

        % Sample Entropy
        r = 0.2;
        AD.SC_Acc_x_SamEn(n,:) = sampen(SC_Acc(:,1),1,r);
        AD.SC_Acc_y_SamEn(n,:) = sampen(SC_Acc(:,2),1,r);
        AD.SC_Acc_z_SamEn(n,:) = sampen(SC_Acc(:,3),1,r);
        AD.SC_Acc_norm_SamEn(n,:) = sampen(SC_Acc_norm,1,r);
        AD.DLS_R_Acc_x_SamEn(n,:) = sampen(DLS_R_Acc(:,1),1,r);
        AD.DLS_R_Acc_y_SamEn(n,:) = sampen(DLS_R_Acc(:,2),1,r);
        AD.DLS_R_Acc_z_SamEn(n,:) = sampen(DLS_R_Acc(:,3),1,r);
        AD.DLS_R_Acc_norm_SamEn(n,:) = sampen(DLS_R_Acc_norm,1,r);
        AD.DLS_L_Acc_x_SamEn(n,:) = sampen(DLS_L_Acc(:,1),1,r);
        AD.DLS_L_Acc_y_SamEn(n,:) = sampen(DLS_L_Acc(:,2),1,r);
        AD.DLS_L_Acc_z_SamEn(n,:) = sampen(DLS_L_Acc(:,3),1,r);
        AD.DLS_L_Acc_norm_SamEn(n,:) = sampen(DLS_L_Acc_norm,1,r);


        % Frequency Domain
        ff = FFeatures(SC_Acc(:,1), Hz);
        AD.SC_Acc_x_DAmp(n,:) = ff(1);
        AD.SC_Acc_x_DFreq(n,:) = ff(2);
        AD.SC_Acc_x_PSD_mean(n,:) = ff(3);
        AD.SC_Acc_x_PSD_std(n,:) = ff(4);
        AD.SC_Acc_x_PSD_skew(n,:) = ff(5);
        AD.SC_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,2), Hz);
        AD.SC_Acc_y_DAmp(n,:) = ff(1);
        AD.SC_Acc_y_DFreq(n,:) = ff(2);
        AD.SC_Acc_y_PSD_mean(n,:) = ff(3);
        AD.SC_Acc_y_PSD_std(n,:) = ff(4);
        AD.SC_Acc_y_PSD_skew(n,:) = ff(5);
        AD.SC_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,3), Hz);
        AD.SC_Acc_z_DAmp(n,:) = ff(1);
        AD.SC_Acc_z_DFreq(n,:) = ff(2);
        AD.SC_Acc_z_PSD_mean(n,:) = ff(3);
        AD.SC_Acc_z_PSD_std(n,:) = ff(4);
        AD.SC_Acc_z_PSD_skew(n,:) = ff(5);
        AD.SC_Acc_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Acc_norm, Hz);
        AD.SC_Acc_norm_DAmp(n,:) = ff(1);
        AD.SC_Acc_norm_DFreq(n,:) = ff(2);
        AD.SC_Acc_norm_PSD_mean(n,:) = ff(3);
        AD.SC_Acc_norm_PSD_std(n,:) = ff(4);
        AD.SC_Acc_norm_PSD_skew(n,:) = ff(5);
        AD.SC_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,1), Hz);
        AD.DLS_R_Acc_x_DAmp(n,:) = ff(1);
        AD.DLS_R_Acc_x_DFreq(n,:) = ff(2);
        AD.DLS_R_Acc_x_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Acc_x_PSD_std(n,:) = ff(4);
        AD.DLS_R_Acc_x_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,2), Hz);
        AD.DLS_R_Acc_y_DAmp(n,:) = ff(1);
        AD.DLS_R_Acc_y_DFreq(n,:) = ff(2);
        AD.DLS_R_Acc_y_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Acc_y_PSD_std(n,:) = ff(4);
        AD.DLS_R_Acc_y_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,3), Hz);
        AD.DLS_R_Acc_z_DAmp(n,:) = ff(1);
        AD.DLS_R_Acc_z_DFreq(n,:) = ff(2);
        AD.DLS_R_Acc_z_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Acc_z_PSD_std(n,:) = ff(4);
        AD.DLS_R_Acc_z_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc_norm, Hz);
        AD.DLS_R_Acc_norm_DAmp(n,:) = ff(1);
        AD.DLS_R_Acc_norm_DFreq(n,:) = ff(2);
        AD.DLS_R_Acc_norm_PSD_mean(n,:) = ff(3);
        AD.DLS_R_Acc_norm_PSD_std(n,:) = ff(4);
        AD.DLS_R_Acc_norm_PSD_skew(n,:) = ff(5);
        AD.DLS_R_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,1), Hz);
        AD.DLS_L_Acc_x_DAmp(n,:) = ff(1);
        AD.DLS_L_Acc_x_DFreq(n,:) = ff(2);
        AD.DLS_L_Acc_x_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Acc_x_PSD_std(n,:) = ff(4);
        AD.DLS_L_Acc_x_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,2), Hz);
        AD.DLS_L_Acc_y_DAmp(n,:) = ff(1);
        AD.DLS_L_Acc_y_DFreq(n,:) = ff(2);
        AD.DLS_L_Acc_y_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Acc_y_PSD_std(n,:) = ff(4);
        AD.DLS_L_Acc_y_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,3), Hz);
        AD.DLS_L_Acc_z_DAmp(n,:) = ff(1);
        AD.DLS_L_Acc_z_DFreq(n,:) = ff(2);
        AD.DLS_L_Acc_z_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Acc_z_PSD_std(n,:) = ff(4);
        AD.DLS_L_Acc_z_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc_norm, Hz);
        AD.DLS_L_Acc_norm_DAmp(n,:) = ff(1);
        AD.DLS_L_Acc_norm_DFreq(n,:) = ff(2);
        AD.DLS_L_Acc_norm_PSD_mean(n,:) = ff(3);
        AD.DLS_L_Acc_norm_PSD_std(n,:) = ff(4);
        AD.DLS_L_Acc_norm_PSD_skew(n,:) = ff(5);
        AD.DLS_L_Acc_norm_PSD_kurtosis(n,:) = ff(6);


        clear data Time NoSteps ...
         SC_Gyr DLS_L_Gyr DLS_R_Gyr SC_Gyr_norm DLS_L_Gyr_norm DLS_R_Gyr_norm ...
         SC_Acc DLS_L_Acc DLS_R_Acc SC_Acc_norm DLS_L_Acc_norm DLS_R_Acc_norm ...
         step_R step_L DLS_L_Gyr_filt DLS_R_Gyr_filt pks_G_L locs_G_L pks_G_R locs_G_R 
      
        
    end
    
    % CVA
    % ID: 1, 25  --> no data
    % ID: 2 --> partially missing DLS_L_Gyr data (only 2561 (81.92 sec))
    % ID: 7 --> sensor at DLS_L missing
    % ID: 32 --> long data 3 times*6MWT? 
    % Outlier --> 19, 33 (not anymore) - corrected
    
    ID_DC = [3:6 8:24 26:31 33:55]
    Type_of_Subject = 'CVA'

    for n = 1:1:length(ID_DC)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID_DC(n)) '.mat']
        load(file_input);

        Side = get_side('./CVA_Paretic_Side.csv');
        AS = Side.Side(ID_DC(n));

        SN = length(data.Session);     % Session #
        TN = 1;     % Trial #
     
        % Gyroscope data
        SC_Gyr = data.Session{SN}.Motion.SC.Gyr{TN};
        DLS_R_Gyr = data.Session{SN}.Motion.DLS_R.Gyr{TN};
        DLS_L_Gyr = data.Session{SN}.Motion.DLS_L.Gyr{TN};

        % Acceleration data
        SC_Acc = data.Session{SN}.Motion.SC.Acc{TN};
        DLS_R_Acc = data.Session{SN}.Motion.DLS_R.Acc{TN};
        DLS_L_Acc = data.Session{SN}.Motion.DLS_L.Acc{TN};

         % To match final time
%              final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
        final = Hz*TS(k);  

        Time = data.Session{SN}.Motion.Time{TN}(1:final,:);
        SC_Gyr = SC_Gyr(1:final,:);
        DLS_R_Gyr = DLS_R_Gyr(1:final,:);
        DLS_L_Gyr = DLS_L_Gyr(1:final,:);

        
        % Norm of Gyro & Acc
        for i = 1:1:length(Time)
            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
            DLS_R_Gyr_norm(i,:) = norm(DLS_R_Gyr(i,:));
            DLS_L_Gyr_norm(i,:) = norm(DLS_L_Gyr(i,:));

            SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
            DLS_R_Acc_norm(i,:) = norm(DLS_R_Acc(i,:));
            DLS_L_Acc_norm(i,:) = norm(DLS_L_Acc(i,:));
        end 

        % Step Count
        f = Hz;
        f_cutoff = 2;
        f_norm = f_cutoff/(f/2);
        ftype = 'low';
        N = 4;  % Order of low pass Butterworth filter
        [B, A] = butter(N,f_norm,ftype);
        DLS_R_Gyr_filt = filtfilt(B, A, DLS_R_Gyr_norm);
        DLS_L_Gyr_filt = filtfilt(B, A, DLS_L_Gyr_norm);


        for i = 1:1:length(DLS_R_Gyr_filt)
            if DLS_R_Gyr_filt(i) <  20
                DLS_R_Gyr_filt(i) = 0;
            else
                DLS_R_Gyr_filt(i) = DLS_R_Gyr_filt(i);
            end

            if DLS_L_Gyr_filt(i) <  20
                DLS_L_Gyr_filt(i) = 0;
            else
                DLS_L_Gyr_filt(i) = DLS_L_Gyr_filt(i);
            end
        end

        [pks_G_R, locs_G_R] = findpeaks(DLS_R_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_R_Gyr_filt)+std(DLS_R_Gyr_filt));
        [pks_G_L, locs_G_L] = findpeaks(DLS_L_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_L_Gyr_filt)+std(DLS_L_Gyr_filt));

        Step_R = length(locs_G_R);
        Step_L = length(locs_G_L);
        
        DC.Steps(n,:) = Step_R + Step_L;

        % Amount of motion
        % Pelvic Amount of motion 
        DC.AoM_Pel_tilt(n,:) = sum(abs(SC_Gyr(:,1))) * dt;
        DC.AoM_Pel_ro(n,:) = sum(abs(SC_Gyr(:,2))) * dt;
        DC.AoM_Pel_oblq(n,:) = sum(abs(SC_Gyr(:,3))) * dt;
        DC.AoM_Pel_norm(n,:) = sum(abs(SC_Gyr_norm)) * dt; 
        
        % Ankle Amount of motion
        if strcmp(AS,'L') == 1
            DC.AoM_Ankle_AS_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;  
            DC.AoM_Ankle_US_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;
            DC.AoM_Ankle_AS_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;  
            DC.AoM_Ankle_US_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;
            DC.AoM_Ankle_AS_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;  
            DC.AoM_Ankle_US_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;
            DC.AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt; 
            DC.AoM_Ankle_US_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt;

        elseif strcmp(AS,'R') == 1
            DC.AoM_Ankle_AS_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;  
            DC.AoM_Ankle_US_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;
            DC.AoM_Ankle_AS_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;  
            DC.AoM_Ankle_US_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;
            DC.AoM_Ankle_AS_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;  
            DC.AoM_Ankle_US_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;
            DC.AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt; 
            DC.AoM_Ankle_US_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt;
        end  

        % General Features
        % Gyro mean
        DC.SC_Gyr_x_mean(n,:) = mean(SC_Gyr(:,1));
        DC.SC_Gyr_y_mean(n,:) = mean(SC_Gyr(:,2));
        DC.SC_Gyr_z_mean(n,:) = mean(SC_Gyr(:,3));
        DC.SC_Gyr_norm_mean(n,:) = mean(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_mean(n,:) = mean(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_mean(n,:) = mean(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_mean(n,:) = mean(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_mean(n,:) = mean(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_mean(n,:) = mean(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_mean(n,:) = mean(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_mean(n,:) = mean(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_mean(n,:) = mean(DLS_L_Gyr_norm);

        % Range
        DC.SC_Gyr_x_range(n,:) = range(SC_Gyr(:,1));
        DC.SC_Gyr_y_range(n,:) = range(SC_Gyr(:,2));
        DC.SC_Gyr_z_range(n,:) = range(SC_Gyr(:,3));
        DC.SC_Gyr_norm_range(n,:) = range(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_range(n,:) = range(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_range(n,:) = range(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_range(n,:) = range(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_range(n,:) = range(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_range(n,:) = range(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_range(n,:) = range(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_range(n,:) = range(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_range(n,:) = range(DLS_L_Gyr_norm);

        % RMS
        DC.SC_Gyr_x_rms(n,:) = rms(SC_Gyr(:,1));
        DC.SC_Gyr_y_rms(n,:) = rms(SC_Gyr(:,2));
        DC.SC_Gyr_z_rms(n,:) = rms(SC_Gyr(:,3));
        DC.SC_Gyr_norm_rms(n,:) = rms(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_rms(n,:) = rms(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_rms(n,:) = rms(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_rms(n,:) = rms(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_rms(n,:) = rms(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_rms(n,:) = rms(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_rms(n,:) = rms(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_rms(n,:) = rms(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_rms(n,:) = rms(DLS_L_Gyr_norm);

        % Standard Deviation
        DC.SC_Gyr_x_std(n,:) = std(SC_Gyr(:,1));
        DC.SC_Gyr_y_std(n,:) = std(SC_Gyr(:,2));
        DC.SC_Gyr_z_std(n,:) = std(SC_Gyr(:,3));
        DC.SC_Gyr_norm_std(n,:) = std(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_std(n,:) = std(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_std(n,:) = std(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_std(n,:) = std(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_std(n,:) = std(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_std(n,:) = std(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_std(n,:) = std(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_std(n,:) = std(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_std(n,:) = std(DLS_L_Gyr_norm);

        % Skew
        DC.SC_Gyr_x_skew(n,:) = skewness(SC_Gyr(:,1));
        DC.SC_Gyr_y_skew(n,:) = skewness(SC_Gyr(:,2));
        DC.SC_Gyr_z_skew(n,:) = skewness(SC_Gyr(:,3));
        DC.SC_Gyr_norm_skew(n,:) = skewness(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_skew(n,:) = skewness(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_skew(n,:) = skewness(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_skew(n,:) = skewness(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_skew(n,:) = skewness(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_skew(n,:) = skewness(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_skew(n,:) = skewness(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_skew(n,:) = skewness(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_skew(n,:) = skewness(DLS_L_Gyr_norm);

        % Kurtosis
        DC.SC_Gyr_x_kurtosis(n,:) = kurtosis(SC_Gyr(:,1));
        DC.SC_Gyr_y_kurtosis(n,:) = kurtosis(SC_Gyr(:,2));
        DC.SC_Gyr_z_kurtosis(n,:) = kurtosis(SC_Gyr(:,3));
        DC.SC_Gyr_norm_kurtosis(n,:) = kurtosis(SC_Gyr_norm);
        DC.DLS_R_Gyr_x_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,1));
        DC.DLS_R_Gyr_y_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_z_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_R_Gyr_norm);
        DC.DLS_L_Gyr_x_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,1));
        DC.DLS_L_Gyr_y_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_z_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_L_Gyr_norm);

        % Pearson correlation coefficient
        DC.SC_Gyr_corr_xy(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,2));
        DC.SC_Gyr_corr_xz(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,3));
        DC.SC_Gyr_corr_yz(n,:) = corr(SC_Gyr(:,2),SC_Gyr(:,3));
        DC.DLS_R_Gyr_corr_xy(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,2));
        DC.DLS_R_Gyr_corr_xz(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,3));
        DC.DLS_R_Gyr_corr_yz(n,:) = corr(DLS_R_Gyr(:,2),DLS_R_Gyr(:,3));
        DC.DLS_L_Gyr_corr_xy(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,2));
        DC.DLS_L_Gyr_corr_xz(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,3));
        DC.DLS_L_Gyr_corr_yz(n,:) = corr(DLS_L_Gyr(:,2),DLS_L_Gyr(:,3));

        % Sample Entropy
        r = 0.2;
        DC.SC_Gyr_x_SamEn(n,:) = sampen(SC_Gyr(:,1),1,r);
        DC.SC_Gyr_y_SamEn(n,:) = sampen(SC_Gyr(:,2),1,r);
        DC.SC_Gyr_z_SamEn(n,:) = sampen(SC_Gyr(:,3),1,r);
        DC.SC_Gyr_norm_SamEn(n,:) = sampen(SC_Gyr_norm,1,r);
        DC.DLS_R_Gyr_x_SamEn(n,:) = sampen(DLS_R_Gyr(:,1),1,r);
        DC.DLS_R_Gyr_y_SamEn(n,:) = sampen(DLS_R_Gyr(:,2),1,r);
        DC.DLS_R_Gyr_z_SamEn(n,:) = sampen(DLS_R_Gyr(:,3),1,r);
        DC.DLS_R_Gyr_norm_SamEn(n,:) = sampen(DLS_R_Gyr_norm,1,r);
        DC.DLS_L_Gyr_x_SamEn(n,:) = sampen(DLS_L_Gyr(:,1),1,r);
        DC.DLS_L_Gyr_y_SamEn(n,:) = sampen(DLS_L_Gyr(:,2),1,r);
        DC.DLS_L_Gyr_z_SamEn(n,:) = sampen(DLS_L_Gyr(:,3),1,r);
        DC.DLS_L_Gyr_norm_SamEn(n,:) = sampen(DLS_L_Gyr_norm,1,r);

        % Frequency Domain
        ff = FFeatures(SC_Gyr(:,1), Hz);
        DC.SC_Gyr_x_DAmp(n,:) = ff(1);
        DC.SC_Gyr_x_DFreq(n,:) = ff(2);
        DC.SC_Gyr_x_PSD_mean(n,:) = ff(3);
        DC.SC_Gyr_x_PSD_std(n,:) = ff(4);
        DC.SC_Gyr_x_PSD_skew(n,:) = ff(5);
        DC.SC_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,2), Hz);
        DC.SC_Gyr_y_DAmp(n,:) = ff(1);
        DC.SC_Gyr_y_DFreq(n,:) = ff(2);
        DC.SC_Gyr_y_PSD_mean(n,:) = ff(3);
        DC.SC_Gyr_y_PSD_std(n,:) = ff(4);
        DC.SC_Gyr_y_PSD_skew(n,:) = ff(5);
        DC.SC_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,3), Hz);
        DC.SC_Gyr_z_DAmp(n,:) = ff(1);
        DC.SC_Gyr_z_DFreq(n,:) = ff(2);
        DC.SC_Gyr_z_PSD_mean(n,:) = ff(3);
        DC.SC_Gyr_z_PSD_std(n,:) = ff(4);
        DC.SC_Gyr_z_PSD_skew(n,:) = ff(5);
        DC.SC_Gyr_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Gyr_norm, Hz);
        DC.SC_Gyr_norm_DAmp(n,:) = ff(1);
        DC.SC_Gyr_norm_DFreq(n,:) = ff(2);
        DC.SC_Gyr_norm_PSD_mean(n,:) = ff(3);
        DC.SC_Gyr_norm_PSD_std(n,:) = ff(4);
        DC.SC_Gyr_norm_PSD_skew(n,:) = ff(5);
        DC.SC_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,1), Hz);
        DC.DLS_R_Gyr_x_DAmp(n,:) = ff(1);
        DC.DLS_R_Gyr_x_DFreq(n,:) = ff(2);
        DC.DLS_R_Gyr_x_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Gyr_x_PSD_std(n,:) = ff(4);
        DC.DLS_R_Gyr_x_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,2), Hz);
        DC.DLS_R_Gyr_y_DAmp(n,:) = ff(1);
        DC.DLS_R_Gyr_y_DFreq(n,:) = ff(2);
        DC.DLS_R_Gyr_y_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Gyr_y_PSD_std(n,:) = ff(4);
        DC.DLS_R_Gyr_y_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,3), Hz);
        DC.DLS_R_Gyr_z_DAmp(n,:) = ff(1);
        DC.DLS_R_Gyr_z_DFreq(n,:) = ff(2);
        DC.DLS_R_Gyr_z_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Gyr_z_PSD_std(n,:) = ff(4);
        DC.DLS_R_Gyr_z_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr_norm, Hz);
        DC.DLS_R_Gyr_norm_DAmp(n,:) = ff(1);
        DC.DLS_R_Gyr_norm_DFreq(n,:) = ff(2);
        DC.DLS_R_Gyr_norm_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Gyr_norm_PSD_std(n,:) = ff(4);
        DC.DLS_R_Gyr_norm_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,1), Hz);
        DC.DLS_L_Gyr_x_DAmp(n,:) = ff(1);
        DC.DLS_L_Gyr_x_DFreq(n,:) = ff(2);
        DC.DLS_L_Gyr_x_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Gyr_x_PSD_std(n,:) = ff(4);
        DC.DLS_L_Gyr_x_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,2), Hz);
        DC.DLS_L_Gyr_y_DAmp(n,:) = ff(1);
        DC.DLS_L_Gyr_y_DFreq(n,:) = ff(2);
        DC.DLS_L_Gyr_y_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Gyr_y_PSD_std(n,:) = ff(4);
        DC.DLS_L_Gyr_y_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,3), Hz);
        DC.DLS_L_Gyr_z_DAmp(n,:) = ff(1);
        DC.DLS_L_Gyr_z_DFreq(n,:) = ff(2);
        DC.DLS_L_Gyr_z_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Gyr_z_PSD_std(n,:) = ff(4);
        DC.DLS_L_Gyr_z_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr_norm, Hz);
        DC.DLS_L_Gyr_norm_DAmp(n,:) = ff(1);
        DC.DLS_L_Gyr_norm_DFreq(n,:) = ff(2);
        DC.DLS_L_Gyr_norm_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Gyr_norm_PSD_std(n,:) = ff(4);
        DC.DLS_L_Gyr_norm_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Gyr_norm_PSD_kurtosis(n,:) = ff(6);


        % Acc mean
        DC.SC_Acc_x_mean(n,:) = mean(SC_Acc(:,1));
        DC.SC_Acc_y_mean(n,:) = mean(SC_Acc(:,2));
        DC.SC_Acc_z_mean(n,:) = mean(SC_Acc(:,3));
        DC.SC_Acc_norm_mean(n,:) = mean(SC_Acc_norm);
        DC.DLS_R_Acc_x_mean(n,:) = mean(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_mean(n,:) = mean(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_mean(n,:) = mean(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_mean(n,:) = mean(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_mean(n,:) = mean(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_mean(n,:) = mean(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_mean(n,:) = mean(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_mean(n,:) = mean(DLS_L_Acc_norm);

        % Range
        DC.SC_Acc_x_range(n,:) = range(SC_Acc(:,1));
        DC.SC_Acc_y_range(n,:) = range(SC_Acc(:,2));
        DC.SC_Acc_z_range(n,:) = range(SC_Acc(:,3));
        DC.SC_Acc_norm_range(n,:) = range(SC_Acc_norm);
        DC.DLS_R_Acc_x_range(n,:) = range(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_range(n,:) = range(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_range(n,:) = range(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_range(n,:) = range(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_range(n,:) = range(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_range(n,:) = range(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_range(n,:) = range(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_range(n,:) = range(DLS_L_Acc_norm);

        % RMS
        DC.SC_Acc_x_rms(n,:) = rms(SC_Acc(:,1));
        DC.SC_Acc_y_rms(n,:) = rms(SC_Acc(:,2));
        DC.SC_Acc_z_rms(n,:) = rms(SC_Acc(:,3));
        DC.SC_Acc_norm_rms(n,:) = rms(SC_Acc_norm);
        DC.DLS_R_Acc_x_rms(n,:) = rms(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_rms(n,:) = rms(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_rms(n,:) = rms(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_rms(n,:) = rms(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_rms(n,:) = rms(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_rms(n,:) = rms(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_rms(n,:) = rms(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_rms(n,:) = rms(DLS_L_Acc_norm);

        % Standard Deviation
        DC.SC_Acc_x_std(n,:) = std(SC_Acc(:,1));
        DC.SC_Acc_y_std(n,:) = std(SC_Acc(:,2));
        DC.SC_Acc_z_std(n,:) = std(SC_Acc(:,3));
        DC.SC_Acc_norm_std(n,:) = std(SC_Acc_norm);
        DC.DLS_R_Acc_x_std(n,:) = std(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_std(n,:) = std(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_std(n,:) = std(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_std(n,:) = std(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_std(n,:) = std(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_std(n,:) = std(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_std(n,:) = std(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_std(n,:) = std(DLS_L_Acc_norm);

        % Skew
        DC.SC_Acc_x_skew(n,:) = skewness(SC_Acc(:,1));
        DC.SC_Acc_y_skew(n,:) = skewness(SC_Acc(:,2));
        DC.SC_Acc_z_skew(n,:) = skewness(SC_Acc(:,3));
        DC.SC_Acc_norm_skew(n,:) = skewness(SC_Acc_norm);
        DC.DLS_R_Acc_x_skew(n,:) = skewness(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_skew(n,:) = skewness(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_skew(n,:) = skewness(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_skew(n,:) = skewness(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_skew(n,:) = skewness(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_skew(n,:) = skewness(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_skew(n,:) = skewness(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_skew(n,:) = skewness(DLS_L_Acc_norm);

        % Kurtosis
        DC.SC_Acc_x_kurtosis(n,:) = kurtosis(SC_Acc(:,1));
        DC.SC_Acc_y_kurtosis(n,:) = kurtosis(SC_Acc(:,2));
        DC.SC_Acc_z_kurtosis(n,:) = kurtosis(SC_Acc(:,3));
        DC.SC_Acc_norm_kurtosis(n,:) = kurtosis(SC_Acc_norm);
        DC.DLS_R_Acc_x_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,1));
        DC.DLS_R_Acc_y_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,2));
        DC.DLS_R_Acc_z_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,3));
        DC.DLS_R_Acc_norm_kurtosis(n,:) = kurtosis(DLS_R_Acc_norm);
        DC.DLS_L_Acc_x_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,1));
        DC.DLS_L_Acc_y_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,2));
        DC.DLS_L_Acc_z_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,3));
        DC.DLS_L_Acc_norm_kurtosis(n,:) = kurtosis(DLS_L_Acc_norm);

        % Pearson correlation coefficient
        DC.SC_Acc_corr_xy(n,:) = corr(SC_Acc(:,1),SC_Acc(:,2));
        DC.SC_Acc_corr_xz(n,:) = corr(SC_Acc(:,1),SC_Acc(:,3));
        DC.SC_Acc_corr_yz(n,:) = corr(SC_Acc(:,2),SC_Acc(:,3));
        DC.DLS_R_Acc_corr_xy(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,2));
        DC.DLS_R_Acc_corr_xz(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,3));
        DC.DLS_R_Acc_corr_yz(n,:) = corr(DLS_R_Acc(:,2),DLS_R_Acc(:,3));
        DC.DLS_L_Acc_corr_xy(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,2));
        DC.DLS_L_Acc_corr_xz(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,3));
        DC.DLS_L_Acc_corr_yz(n,:) = corr(DLS_L_Acc(:,2),DLS_L_Acc(:,3));

        % Sample Entropy
        r = 0.2;
        DC.SC_Acc_x_SamEn(n,:) = sampen(SC_Acc(:,1),1,r);
        DC.SC_Acc_y_SamEn(n,:) = sampen(SC_Acc(:,2),1,r);
        DC.SC_Acc_z_SamEn(n,:) = sampen(SC_Acc(:,3),1,r);
        DC.SC_Acc_norm_SamEn(n,:) = sampen(SC_Acc_norm,1,r);
        DC.DLS_R_Acc_x_SamEn(n,:) = sampen(DLS_R_Acc(:,1),1,r);
        DC.DLS_R_Acc_y_SamEn(n,:) = sampen(DLS_R_Acc(:,2),1,r);
        DC.DLS_R_Acc_z_SamEn(n,:) = sampen(DLS_R_Acc(:,3),1,r);
        DC.DLS_R_Acc_norm_SamEn(n,:) = sampen(DLS_R_Acc_norm,1,r);
        DC.DLS_L_Acc_x_SamEn(n,:) = sampen(DLS_L_Acc(:,1),1,r);
        DC.DLS_L_Acc_y_SamEn(n,:) = sampen(DLS_L_Acc(:,2),1,r);
        DC.DLS_L_Acc_z_SamEn(n,:) = sampen(DLS_L_Acc(:,3),1,r);
        DC.DLS_L_Acc_norm_SamEn(n,:) = sampen(DLS_L_Acc_norm,1,r);

        % Frequency Domain
        ff = FFeatures(SC_Acc(:,1), Hz);
        DC.SC_Acc_x_DAmp(n,:) = ff(1);
        DC.SC_Acc_x_DFreq(n,:) = ff(2);
        DC.SC_Acc_x_PSD_mean(n,:) = ff(3);
        DC.SC_Acc_x_PSD_std(n,:) = ff(4);
        DC.SC_Acc_x_PSD_skew(n,:) = ff(5);
        DC.SC_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,2), Hz);
        DC.SC_Acc_y_DAmp(n,:) = ff(1);
        DC.SC_Acc_y_DFreq(n,:) = ff(2);
        DC.SC_Acc_y_PSD_mean(n,:) = ff(3);
        DC.SC_Acc_y_PSD_std(n,:) = ff(4);
        DC.SC_Acc_y_PSD_skew(n,:) = ff(5);
        DC.SC_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,3), Hz);
        DC.SC_Acc_z_DAmp(n,:) = ff(1);
        DC.SC_Acc_z_DFreq(n,:) = ff(2);
        DC.SC_Acc_z_PSD_mean(n,:) = ff(3);
        DC.SC_Acc_z_PSD_std(n,:) = ff(4);
        DC.SC_Acc_z_PSD_skew(n,:) = ff(5);
        DC.SC_Acc_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Acc_norm, Hz);
        DC.SC_Acc_norm_DAmp(n,:) = ff(1);
        DC.SC_Acc_norm_DFreq(n,:) = ff(2);
        DC.SC_Acc_norm_PSD_mean(n,:) = ff(3);
        DC.SC_Acc_norm_PSD_std(n,:) = ff(4);
        DC.SC_Acc_norm_PSD_skew(n,:) = ff(5);
        DC.SC_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,1), Hz);
        DC.DLS_R_Acc_x_DAmp(n,:) = ff(1);
        DC.DLS_R_Acc_x_DFreq(n,:) = ff(2);
        DC.DLS_R_Acc_x_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Acc_x_PSD_std(n,:) = ff(4);
        DC.DLS_R_Acc_x_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,2), Hz);
        DC.DLS_R_Acc_y_DAmp(n,:) = ff(1);
        DC.DLS_R_Acc_y_DFreq(n,:) = ff(2);
        DC.DLS_R_Acc_y_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Acc_y_PSD_std(n,:) = ff(4);
        DC.DLS_R_Acc_y_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,3), Hz);
        DC.DLS_R_Acc_z_DAmp(n,:) = ff(1);
        DC.DLS_R_Acc_z_DFreq(n,:) = ff(2);
        DC.DLS_R_Acc_z_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Acc_z_PSD_std(n,:) = ff(4);
        DC.DLS_R_Acc_z_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc_norm, Hz);
        DC.DLS_R_Acc_norm_DAmp(n,:) = ff(1);
        DC.DLS_R_Acc_norm_DFreq(n,:) = ff(2);
        DC.DLS_R_Acc_norm_PSD_mean(n,:) = ff(3);
        DC.DLS_R_Acc_norm_PSD_std(n,:) = ff(4);
        DC.DLS_R_Acc_norm_PSD_skew(n,:) = ff(5);
        DC.DLS_R_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,1), Hz);
        DC.DLS_L_Acc_x_DAmp(n,:) = ff(1);
        DC.DLS_L_Acc_x_DFreq(n,:) = ff(2);
        DC.DLS_L_Acc_x_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Acc_x_PSD_std(n,:) = ff(4);
        DC.DLS_L_Acc_x_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,2), Hz);
        DC.DLS_L_Acc_y_DAmp(n,:) = ff(1);
        DC.DLS_L_Acc_y_DFreq(n,:) = ff(2);
        DC.DLS_L_Acc_y_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Acc_y_PSD_std(n,:) = ff(4);
        DC.DLS_L_Acc_y_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,3), Hz);
        DC.DLS_L_Acc_z_DAmp(n,:) = ff(1);
        DC.DLS_L_Acc_z_DFreq(n,:) = ff(2);
        DC.DLS_L_Acc_z_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Acc_z_PSD_std(n,:) = ff(4);
        DC.DLS_L_Acc_z_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc_norm, Hz);
        DC.DLS_L_Acc_norm_DAmp(n,:) = ff(1);
        DC.DLS_L_Acc_norm_DFreq(n,:) = ff(2);
        DC.DLS_L_Acc_norm_PSD_mean(n,:) = ff(3);
        DC.DLS_L_Acc_norm_PSD_std(n,:) = ff(4);
        DC.DLS_L_Acc_norm_PSD_skew(n,:) = ff(5);
        DC.DLS_L_Acc_norm_PSD_kurtosis(n,:) = ff(6);


         clear data Time NoSteps ...
         SC_Gyr DLS_L_Gyr DLS_R_Gyr SC_Gyr_norm DLS_L_Gyr_norm DLS_R_Gyr_norm ...
         SC_Acc DLS_L_Acc DLS_R_Acc SC_Acc_norm DLS_L_Acc_norm DLS_R_Acc_norm ...
         step_R step_L DLS_L_Gyr_filt DLS_R_Gyr_filt pks_G_L locs_G_L pks_G_R locs_G_R 

    end
  
    
    % Healthy Controls (HC)
    % ID: 26, 38 --> no data
    % ID: 4 --> no sensor at DLS_L
    
    ID_HC = [1:3 5:25 27:37 39:51] 
    Type_of_Subject = 'CONTROLS'
    
    for n = 1:1:length(ID_HC)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID_HC(n)) '.mat']
        load(file_input);

        SN = length(data.Session);     % Session #
        TN = 1;     % Trial #

        % Gyroscope data
        SC_Gyr = data.Session{SN}.Motion.SC.Gyr{TN};
        DLS_R_Gyr = data.Session{SN}.Motion.DLS_R.Gyr{TN};
        DLS_L_Gyr = data.Session{SN}.Motion.DLS_L.Gyr{TN};

        % Acceleration data
        SC_Acc = data.Session{SN}.Motion.SC.Acc{TN};
        DLS_R_Acc = data.Session{SN}.Motion.DLS_R.Acc{TN};
        DLS_L_Acc = data.Session{SN}.Motion.DLS_L.Acc{TN};

         % To match final time
%              final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
        final = Hz*TS(k);  

        Time = data.Session{SN}.Motion.Time{TN}(1:final,:);
        SC_Gyr = SC_Gyr(1:final,:);
        DLS_R_Gyr = DLS_R_Gyr(1:final,:);
        DLS_L_Gyr = DLS_L_Gyr(1:final,:);

        
        % Norm of Gyro & Acc
        for i = 1:1:length(Time)
            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
            DLS_R_Gyr_norm(i,:) = norm(DLS_R_Gyr(i,:));
            DLS_L_Gyr_norm(i,:) = norm(DLS_L_Gyr(i,:));

            SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
            DLS_R_Acc_norm(i,:) = norm(DLS_R_Acc(i,:));
            DLS_L_Acc_norm(i,:) = norm(DLS_L_Acc(i,:));
        end 

        % Step Count
        f = Hz;
        f_cutoff = 2;
        f_norm = f_cutoff/(f/2);
        ftype = 'low';
        N = 4;  % Order of low pass Butterworth filter
        [B, A] = butter(N,f_norm,ftype);
        DLS_R_Gyr_filt = filtfilt(B, A, DLS_R_Gyr_norm);
        DLS_L_Gyr_filt = filtfilt(B, A, DLS_L_Gyr_norm);


        for i = 1:1:length(DLS_R_Gyr_filt)
            if DLS_R_Gyr_filt(i) <  20
                DLS_R_Gyr_filt(i) = 0;
            else
                DLS_R_Gyr_filt(i) = DLS_R_Gyr_filt(i);
            end

            if DLS_L_Gyr_filt(i) <  20
                DLS_L_Gyr_filt(i) = 0;
            else
                DLS_L_Gyr_filt(i) = DLS_L_Gyr_filt(i);
            end
        end

        [pks_G_R, locs_G_R] = findpeaks(DLS_R_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_R_Gyr_filt)+std(DLS_R_Gyr_filt));
        [pks_G_L, locs_G_L] = findpeaks(DLS_L_Gyr_filt, 'MinPeakDistance',20,'MinPeakHeight',mean(DLS_L_Gyr_filt)+std(DLS_L_Gyr_filt));

        Step_R = length(locs_G_R);
        Step_L = length(locs_G_L);
        
        HC.Steps(n,:) = Step_R + Step_L;

        % Amount of motion
        % Pelvic Amount of motion 
        HC.AoM_Pel_tilt(n,:) = sum(abs(SC_Gyr(:,1))) * dt;
        HC.AoM_Pel_ro(n,:) = sum(abs(SC_Gyr(:,2))) * dt;
        HC.AoM_Pel_oblq(n,:) = sum(abs(SC_Gyr(:,3))) * dt;
        HC.AoM_Pel_norm(n,:) = sum(abs(SC_Gyr_norm)) * dt; 
        
        % Ankle sensor
        HC.AoM_Ankle_L_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;  
        HC.AoM_Ankle_R_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;
        HC.AoM_Ankle_L_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;  
        HC.AoM_Ankle_R_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;
        HC.AoM_Ankle_L_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;  
        HC.AoM_Ankle_R_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;
        HC.AoM_Ankle_L_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt; 
        HC.AoM_Ankle_R_norm (n,:) = sum(abs(DLS_R_Gyr_norm)) * dt;
        
        % General Features
        % Gyro mean
        HC.SC_Gyr_x_mean(n,:) = mean(SC_Gyr(:,1));
        HC.SC_Gyr_y_mean(n,:) = mean(SC_Gyr(:,2));
        HC.SC_Gyr_z_mean(n,:) = mean(SC_Gyr(:,3));
        HC.SC_Gyr_norm_mean(n,:) = mean(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_mean(n,:) = mean(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_mean(n,:) = mean(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_mean(n,:) = mean(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_mean(n,:) = mean(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_mean(n,:) = mean(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_mean(n,:) = mean(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_mean(n,:) = mean(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_mean(n,:) = mean(DLS_L_Gyr_norm);

        % Range
        HC.SC_Gyr_x_range(n,:) = range(SC_Gyr(:,1));
        HC.SC_Gyr_y_range(n,:) = range(SC_Gyr(:,2));
        HC.SC_Gyr_z_range(n,:) = range(SC_Gyr(:,3));
        HC.SC_Gyr_norm_range(n,:) = range(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_range(n,:) = range(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_range(n,:) = range(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_range(n,:) = range(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_range(n,:) = range(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_range(n,:) = range(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_range(n,:) = range(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_range(n,:) = range(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_range(n,:) = range(DLS_L_Gyr_norm);

        % RMS
        HC.SC_Gyr_x_rms(n,:) = rms(SC_Gyr(:,1));
        HC.SC_Gyr_y_rms(n,:) = rms(SC_Gyr(:,2));
        HC.SC_Gyr_z_rms(n,:) = rms(SC_Gyr(:,3));
        HC.SC_Gyr_norm_rms(n,:) = rms(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_rms(n,:) = rms(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_rms(n,:) = rms(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_rms(n,:) = rms(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_rms(n,:) = rms(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_rms(n,:) = rms(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_rms(n,:) = rms(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_rms(n,:) = rms(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_rms(n,:) = rms(DLS_L_Gyr_norm);

        % Standard Deviation
        HC.SC_Gyr_x_std(n,:) = std(SC_Gyr(:,1));
        HC.SC_Gyr_y_std(n,:) = std(SC_Gyr(:,2));
        HC.SC_Gyr_z_std(n,:) = std(SC_Gyr(:,3));
        HC.SC_Gyr_norm_std(n,:) = std(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_std(n,:) = std(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_std(n,:) = std(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_std(n,:) = std(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_std(n,:) = std(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_std(n,:) = std(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_std(n,:) = std(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_std(n,:) = std(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_std(n,:) = std(DLS_L_Gyr_norm);

        % Skew
        HC.SC_Gyr_x_skew(n,:) = skewness(SC_Gyr(:,1));
        HC.SC_Gyr_y_skew(n,:) = skewness(SC_Gyr(:,2));
        HC.SC_Gyr_z_skew(n,:) = skewness(SC_Gyr(:,3));
        HC.SC_Gyr_norm_skew(n,:) = skewness(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_skew(n,:) = skewness(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_skew(n,:) = skewness(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_skew(n,:) = skewness(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_skew(n,:) = skewness(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_skew(n,:) = skewness(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_skew(n,:) = skewness(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_skew(n,:) = skewness(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_skew(n,:) = skewness(DLS_L_Gyr_norm);

        % Kurtosis
        HC.SC_Gyr_x_kurtosis(n,:) = kurtosis(SC_Gyr(:,1));
        HC.SC_Gyr_y_kurtosis(n,:) = kurtosis(SC_Gyr(:,2));
        HC.SC_Gyr_z_kurtosis(n,:) = kurtosis(SC_Gyr(:,3));
        HC.SC_Gyr_norm_kurtosis(n,:) = kurtosis(SC_Gyr_norm);
        HC.DLS_R_Gyr_x_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,1));
        HC.DLS_R_Gyr_y_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_z_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_R_Gyr_norm);
        HC.DLS_L_Gyr_x_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,1));
        HC.DLS_L_Gyr_y_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_z_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_L_Gyr_norm);

        % Pearson correlation coefficient
        HC.SC_Gyr_corr_xy(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,2));
        HC.SC_Gyr_corr_xz(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,3));
        HC.SC_Gyr_corr_yz(n,:) = corr(SC_Gyr(:,2),SC_Gyr(:,3));
        HC.DLS_R_Gyr_corr_xy(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,2));
        HC.DLS_R_Gyr_corr_xz(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,3));
        HC.DLS_R_Gyr_corr_yz(n,:) = corr(DLS_R_Gyr(:,2),DLS_R_Gyr(:,3));
        HC.DLS_L_Gyr_corr_xy(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,2));
        HC.DLS_L_Gyr_corr_xz(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,3));
        HC.DLS_L_Gyr_corr_yz(n,:) = corr(DLS_L_Gyr(:,2),DLS_L_Gyr(:,3));

        % Sample Entropy
        r = 0.2;
        HC.SC_Gyr_x_SamEn(n,:) = sampen(SC_Gyr(:,1),1,r);
        HC.SC_Gyr_y_SamEn(n,:) = sampen(SC_Gyr(:,2),1,r);
        HC.SC_Gyr_z_SamEn(n,:) = sampen(SC_Gyr(:,3),1,r);
        HC.SC_Gyr_norm_SamEn(n,:) = sampen(SC_Gyr_norm,1,r);
        HC.DLS_R_Gyr_x_SamEn(n,:) = sampen(DLS_R_Gyr(:,1),1,r);
        HC.DLS_R_Gyr_y_SamEn(n,:) = sampen(DLS_R_Gyr(:,2),1,r);
        HC.DLS_R_Gyr_z_SamEn(n,:) = sampen(DLS_R_Gyr(:,3),1,r);
        HC.DLS_R_Gyr_norm_SamEn(n,:) = sampen(DLS_R_Gyr_norm,1,r);
        HC.DLS_L_Gyr_x_SamEn(n,:) = sampen(DLS_L_Gyr(:,1),1,r);
        HC.DLS_L_Gyr_y_SamEn(n,:) = sampen(DLS_L_Gyr(:,2),1,r);
        HC.DLS_L_Gyr_z_SamEn(n,:) = sampen(DLS_L_Gyr(:,3),1,r);
        HC.DLS_L_Gyr_norm_SamEn(n,:) = sampen(DLS_L_Gyr_norm,1,r);

        % Frequency Domain
        ff = FFeatures(SC_Gyr(:,1), Hz);
        HC.SC_Gyr_x_DAmp(n,:) = ff(1);
        HC.SC_Gyr_x_DFreq(n,:) = ff(2);
        HC.SC_Gyr_x_PSD_mean(n,:) = ff(3);
        HC.SC_Gyr_x_PSD_std(n,:) = ff(4);
        HC.SC_Gyr_x_PSD_skew(n,:) = ff(5);
        HC.SC_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,2), Hz);
        HC.SC_Gyr_y_DAmp(n,:) = ff(1);
        HC.SC_Gyr_y_DFreq(n,:) = ff(2);
        HC.SC_Gyr_y_PSD_mean(n,:) = ff(3);
        HC.SC_Gyr_y_PSD_std(n,:) = ff(4);
        HC.SC_Gyr_y_PSD_skew(n,:) = ff(5);
        HC.SC_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Gyr(:,3), Hz);
        HC.SC_Gyr_z_DAmp(n,:) = ff(1);
        HC.SC_Gyr_z_DFreq(n,:) = ff(2);
        HC.SC_Gyr_z_PSD_mean(n,:) = ff(3);
        HC.SC_Gyr_z_PSD_std(n,:) = ff(4);
        HC.SC_Gyr_z_PSD_skew(n,:) = ff(5);
        HC.SC_Gyr_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Gyr_norm, Hz);
        HC.SC_Gyr_norm_DAmp(n,:) = ff(1);
        HC.SC_Gyr_norm_DFreq(n,:) = ff(2);
        HC.SC_Gyr_norm_PSD_mean(n,:) = ff(3);
        HC.SC_Gyr_norm_PSD_std(n,:) = ff(4);
        HC.SC_Gyr_norm_PSD_skew(n,:) = ff(5);
        HC.SC_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,1), Hz);
        HC.DLS_R_Gyr_x_DAmp(n,:) = ff(1);
        HC.DLS_R_Gyr_x_DFreq(n,:) = ff(2);
        HC.DLS_R_Gyr_x_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Gyr_x_PSD_std(n,:) = ff(4);
        HC.DLS_R_Gyr_x_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,2), Hz);
        HC.DLS_R_Gyr_y_DAmp(n,:) = ff(1);
        HC.DLS_R_Gyr_y_DFreq(n,:) = ff(2);
        HC.DLS_R_Gyr_y_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Gyr_y_PSD_std(n,:) = ff(4);
        HC.DLS_R_Gyr_y_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr(:,3), Hz);
        HC.DLS_R_Gyr_z_DAmp(n,:) = ff(1);
        HC.DLS_R_Gyr_z_DFreq(n,:) = ff(2);
        HC.DLS_R_Gyr_z_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Gyr_z_PSD_std(n,:) = ff(4);
        HC.DLS_R_Gyr_z_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Gyr_norm, Hz);
        HC.DLS_R_Gyr_norm_DAmp(n,:) = ff(1);
        HC.DLS_R_Gyr_norm_DFreq(n,:) = ff(2);
        HC.DLS_R_Gyr_norm_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Gyr_norm_PSD_std(n,:) = ff(4);
        HC.DLS_R_Gyr_norm_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,1), Hz);
        HC.DLS_L_Gyr_x_DAmp(n,:) = ff(1);
        HC.DLS_L_Gyr_x_DFreq(n,:) = ff(2);
        HC.DLS_L_Gyr_x_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Gyr_x_PSD_std(n,:) = ff(4);
        HC.DLS_L_Gyr_x_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Gyr_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,2), Hz);
        HC.DLS_L_Gyr_y_DAmp(n,:) = ff(1);
        HC.DLS_L_Gyr_y_DFreq(n,:) = ff(2);
        HC.DLS_L_Gyr_y_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Gyr_y_PSD_std(n,:) = ff(4);
        HC.DLS_L_Gyr_y_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Gyr_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr(:,3), Hz);
        HC.DLS_L_Gyr_z_DAmp(n,:) = ff(1);
        HC.DLS_L_Gyr_z_DFreq(n,:) = ff(2);
        HC.DLS_L_Gyr_z_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Gyr_z_PSD_std(n,:) = ff(4);
        HC.DLS_L_Gyr_z_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Gyr_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Gyr_norm, Hz);
        HC.DLS_L_Gyr_norm_DAmp(n,:) = ff(1);
        HC.DLS_L_Gyr_norm_DFreq(n,:) = ff(2);
        HC.DLS_L_Gyr_norm_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Gyr_norm_PSD_std(n,:) = ff(4);
        HC.DLS_L_Gyr_norm_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Gyr_norm_PSD_kurtosis(n,:) = ff(6);


        % Acc mean
        HC.SC_Acc_x_mean(n,:) = mean(SC_Acc(:,1));
        HC.SC_Acc_y_mean(n,:) = mean(SC_Acc(:,2));
        HC.SC_Acc_z_mean(n,:) = mean(SC_Acc(:,3));
        HC.SC_Acc_norm_mean(n,:) = mean(SC_Acc_norm);
        HC.DLS_R_Acc_x_mean(n,:) = mean(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_mean(n,:) = mean(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_mean(n,:) = mean(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_mean(n,:) = mean(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_mean(n,:) = mean(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_mean(n,:) = mean(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_mean(n,:) = mean(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_mean(n,:) = mean(DLS_L_Acc_norm);

        % Range
        HC.SC_Acc_x_range(n,:) = range(SC_Acc(:,1));
        HC.SC_Acc_y_range(n,:) = range(SC_Acc(:,2));
        HC.SC_Acc_z_range(n,:) = range(SC_Acc(:,3));
        HC.SC_Acc_norm_range(n,:) = range(SC_Acc_norm);
        HC.DLS_R_Acc_x_range(n,:) = range(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_range(n,:) = range(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_range(n,:) = range(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_range(n,:) = range(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_range(n,:) = range(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_range(n,:) = range(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_range(n,:) = range(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_range(n,:) = range(DLS_L_Acc_norm);

        % RMS
        HC.SC_Acc_x_rms(n,:) = rms(SC_Acc(:,1));
        HC.SC_Acc_y_rms(n,:) = rms(SC_Acc(:,2));
        HC.SC_Acc_z_rms(n,:) = rms(SC_Acc(:,3));
        HC.SC_Acc_norm_rms(n,:) = rms(SC_Acc_norm);
        HC.DLS_R_Acc_x_rms(n,:) = rms(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_rms(n,:) = rms(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_rms(n,:) = rms(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_rms(n,:) = rms(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_rms(n,:) = rms(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_rms(n,:) = rms(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_rms(n,:) = rms(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_rms(n,:) = rms(DLS_L_Acc_norm);

        % Standard Deviation
        HC.SC_Acc_x_std(n,:) = std(SC_Acc(:,1));
        HC.SC_Acc_y_std(n,:) = std(SC_Acc(:,2));
        HC.SC_Acc_z_std(n,:) = std(SC_Acc(:,3));
        HC.SC_Acc_norm_std(n,:) = std(SC_Acc_norm);
        HC.DLS_R_Acc_x_std(n,:) = std(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_std(n,:) = std(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_std(n,:) = std(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_std(n,:) = std(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_std(n,:) = std(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_std(n,:) = std(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_std(n,:) = std(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_std(n,:) = std(DLS_L_Acc_norm);

        % Skew
        HC.SC_Acc_x_skew(n,:) = skewness(SC_Acc(:,1));
        HC.SC_Acc_y_skew(n,:) = skewness(SC_Acc(:,2));
        HC.SC_Acc_z_skew(n,:) = skewness(SC_Acc(:,3));
        HC.SC_Acc_norm_skew(n,:) = skewness(SC_Acc_norm);
        HC.DLS_R_Acc_x_skew(n,:) = skewness(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_skew(n,:) = skewness(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_skew(n,:) = skewness(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_skew(n,:) = skewness(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_skew(n,:) = skewness(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_skew(n,:) = skewness(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_skew(n,:) = skewness(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_skew(n,:) = skewness(DLS_L_Acc_norm);

        % Kurtosis
        HC.SC_Acc_x_kurtosis(n,:) = kurtosis(SC_Acc(:,1));
        HC.SC_Acc_y_kurtosis(n,:) = kurtosis(SC_Acc(:,2));
        HC.SC_Acc_z_kurtosis(n,:) = kurtosis(SC_Acc(:,3));
        HC.SC_Acc_norm_kurtosis(n,:) = kurtosis(SC_Acc_norm);
        HC.DLS_R_Acc_x_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,1));
        HC.DLS_R_Acc_y_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,2));
        HC.DLS_R_Acc_z_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,3));
        HC.DLS_R_Acc_norm_kurtosis(n,:) = kurtosis(DLS_R_Acc_norm);
        HC.DLS_L_Acc_x_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,1));
        HC.DLS_L_Acc_y_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,2));
        HC.DLS_L_Acc_z_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,3));
        HC.DLS_L_Acc_norm_kurtosis(n,:) = kurtosis(DLS_L_Acc_norm);

        % Pearson correlation coefficient
        HC.SC_Acc_corr_xy(n,:) = corr(SC_Acc(:,1),SC_Acc(:,2));
        HC.SC_Acc_corr_xz(n,:) = corr(SC_Acc(:,1),SC_Acc(:,3));
        HC.SC_Acc_corr_yz(n,:) = corr(SC_Acc(:,2),SC_Acc(:,3));
        HC.DLS_R_Acc_corr_xy(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,2));
        HC.DLS_R_Acc_corr_xz(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,3));
        HC.DLS_R_Acc_corr_yz(n,:) = corr(DLS_R_Acc(:,2),DLS_R_Acc(:,3));
        HC.DLS_L_Acc_corr_xy(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,2));
        HC.DLS_L_Acc_corr_xz(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,3));
        HC.DLS_L_Acc_corr_yz(n,:) = corr(DLS_L_Acc(:,2),DLS_L_Acc(:,3));

        % Sample Entropy
        r = 0.2;
        HC.SC_Acc_x_SamEn(n,:) = sampen(SC_Acc(:,1),1,r);
        HC.SC_Acc_y_SamEn(n,:) = sampen(SC_Acc(:,2),1,r);
        HC.SC_Acc_z_SamEn(n,:) = sampen(SC_Acc(:,3),1,r);
        HC.SC_Acc_norm_SamEn(n,:) = sampen(SC_Acc_norm,1,r);
        HC.DLS_R_Acc_x_SamEn(n,:) = sampen(DLS_R_Acc(:,1),1,r);
        HC.DLS_R_Acc_y_SamEn(n,:) = sampen(DLS_R_Acc(:,2),1,r);
        HC.DLS_R_Acc_z_SamEn(n,:) = sampen(DLS_R_Acc(:,3),1,r);
        HC.DLS_R_Acc_norm_SamEn(n,:) = sampen(DLS_R_Acc_norm,1,r);
        HC.DLS_L_Acc_x_SamEn(n,:) = sampen(DLS_L_Acc(:,1),1,r);
        HC.DLS_L_Acc_y_SamEn(n,:) = sampen(DLS_L_Acc(:,2),1,r);
        HC.DLS_L_Acc_z_SamEn(n,:) = sampen(DLS_L_Acc(:,3),1,r);
        HC.DLS_L_Acc_norm_SamEn(n,:) = sampen(DLS_L_Acc_norm,1,r);

        % Frequency Domain
        ff = FFeatures(SC_Acc(:,1), Hz);
        HC.SC_Acc_x_DAmp(n,:) = ff(1);
        HC.SC_Acc_x_DFreq(n,:) = ff(2);
        HC.SC_Acc_x_PSD_mean(n,:) = ff(3);
        HC.SC_Acc_x_PSD_std(n,:) = ff(4);
        HC.SC_Acc_x_PSD_skew(n,:) = ff(5);
        HC.SC_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,2), Hz);
        HC.SC_Acc_y_DAmp(n,:) = ff(1);
        HC.SC_Acc_y_DFreq(n,:) = ff(2);
        HC.SC_Acc_y_PSD_mean(n,:) = ff(3);
        HC.SC_Acc_y_PSD_std(n,:) = ff(4);
        HC.SC_Acc_y_PSD_skew(n,:) = ff(5);
        HC.SC_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(SC_Acc(:,3), Hz);
        HC.SC_Acc_z_DAmp(n,:) = ff(1);
        HC.SC_Acc_z_DFreq(n,:) = ff(2);
        HC.SC_Acc_z_PSD_mean(n,:) = ff(3);
        HC.SC_Acc_z_PSD_std(n,:) = ff(4);
        HC.SC_Acc_z_PSD_skew(n,:) = ff(5);
        HC.SC_Acc_z_PSD_kurtosis(n,:) =ff(6);

        ff = FFeatures(SC_Acc_norm, Hz);
        HC.SC_Acc_norm_DAmp(n,:) = ff(1);
        HC.SC_Acc_norm_DFreq(n,:) = ff(2);
        HC.SC_Acc_norm_PSD_mean(n,:) = ff(3);
        HC.SC_Acc_norm_PSD_std(n,:) = ff(4);
        HC.SC_Acc_norm_PSD_skew(n,:) = ff(5);
        HC.SC_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,1), Hz);
        HC.DLS_R_Acc_x_DAmp(n,:) = ff(1);
        HC.DLS_R_Acc_x_DFreq(n,:) = ff(2);
        HC.DLS_R_Acc_x_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Acc_x_PSD_std(n,:) = ff(4);
        HC.DLS_R_Acc_x_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,2), Hz);
        HC.DLS_R_Acc_y_DAmp(n,:) = ff(1);
        HC.DLS_R_Acc_y_DFreq(n,:) = ff(2);
        HC.DLS_R_Acc_y_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Acc_y_PSD_std(n,:) = ff(4);
        HC.DLS_R_Acc_y_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc(:,3), Hz);
        HC.DLS_R_Acc_z_DAmp(n,:) = ff(1);
        HC.DLS_R_Acc_z_DFreq(n,:) = ff(2);
        HC.DLS_R_Acc_z_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Acc_z_PSD_std(n,:) = ff(4);
        HC.DLS_R_Acc_z_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_R_Acc_norm, Hz);
        HC.DLS_R_Acc_norm_DAmp(n,:) = ff(1);
        HC.DLS_R_Acc_norm_DFreq(n,:) = ff(2);
        HC.DLS_R_Acc_norm_PSD_mean(n,:) = ff(3);
        HC.DLS_R_Acc_norm_PSD_std(n,:) = ff(4);
        HC.DLS_R_Acc_norm_PSD_skew(n,:) = ff(5);
        HC.DLS_R_Acc_norm_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,1), Hz);
        HC.DLS_L_Acc_x_DAmp(n,:) = ff(1);
        HC.DLS_L_Acc_x_DFreq(n,:) = ff(2);
        HC.DLS_L_Acc_x_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Acc_x_PSD_std(n,:) = ff(4);
        HC.DLS_L_Acc_x_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Acc_x_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,2), Hz);
        HC.DLS_L_Acc_y_DAmp(n,:) = ff(1);
        HC.DLS_L_Acc_y_DFreq(n,:) = ff(2);
        HC.DLS_L_Acc_y_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Acc_y_PSD_std(n,:) = ff(4);
        HC.DLS_L_Acc_y_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Acc_y_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc(:,3), Hz);
        HC.DLS_L_Acc_z_DAmp(n,:) = ff(1);
        HC.DLS_L_Acc_z_DFreq(n,:) = ff(2);
        HC.DLS_L_Acc_z_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Acc_z_PSD_std(n,:) = ff(4);
        HC.DLS_L_Acc_z_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Acc_z_PSD_kurtosis(n,:) = ff(6);

        ff = FFeatures(DLS_L_Acc_norm, Hz);
        HC.DLS_L_Acc_norm_DAmp(n,:) = ff(1);
        HC.DLS_L_Acc_norm_DFreq(n,:) = ff(2);
        HC.DLS_L_Acc_norm_PSD_mean(n,:) = ff(3);
        HC.DLS_L_Acc_norm_PSD_std(n,:) = ff(4);
        HC.DLS_L_Acc_norm_PSD_skew(n,:) = ff(5);
        HC.DLS_L_Acc_norm_PSD_kurtosis(n,:) = ff(6);

         clear data Time NoSteps ...
         SC_Gyr DLS_L_Gyr DLS_R_Gyr SC_Gyr_norm DLS_L_Gyr_norm DLS_R_Gyr_norm ...
         SC_Acc DLS_L_Acc DLS_R_Acc SC_Acc_norm DLS_L_Acc_norm DLS_R_Acc_norm ...
         step_R step_L DLS_L_Gyr_filt DLS_R_Gyr_filt pks_G_L locs_G_L pks_G_R locs_G_R 


    end
    
    Steps_ = [AD.Steps; DC.Steps; HC.Steps];

    % Amount of motion
    % Pelvic Amount of motion 
    AoM_Pel_tilt_ = [AD.AoM_Pel_tilt; DC.AoM_Pel_tilt; HC.AoM_Pel_tilt]
    AoM_Pel_ro_ = [AD.AoM_Pel_ro; DC.AoM_Pel_ro; HC.AoM_Pel_ro]
    AoM_Pel_oblq_ = [AD.AoM_Pel_oblq; DC.AoM_Pel_oblq; HC.AoM_Pel_oblq]
    AoM_Pel_norm_ = [AD.AoM_Pel_oblq; DC.AoM_Pel_oblq; HC.AoM_Pel_oblq]
    AoM_Ankle_AS_x_ = [AD.AoM_Ankle_AS_x; DC.AoM_Ankle_AS_x; HC.AoM_Ankle_L_x]  
    AoM_Ankle_US_x_ = [AD.AoM_Ankle_US_x; DC.AoM_Ankle_US_x; HC.AoM_Ankle_R_x]  
    AoM_Ankle_AS_y_ = [AD.AoM_Ankle_AS_y; DC.AoM_Ankle_AS_y; HC.AoM_Ankle_L_y]  
    AoM_Ankle_US_y_ = [AD.AoM_Ankle_US_y; DC.AoM_Ankle_US_y; HC.AoM_Ankle_R_y]  
    AoM_Ankle_AS_z_ = [AD.AoM_Ankle_AS_z; DC.AoM_Ankle_AS_z; HC.AoM_Ankle_L_z]  
    AoM_Ankle_US_z_ = [AD.AoM_Ankle_US_z; DC.AoM_Ankle_US_z; HC.AoM_Ankle_R_z]  
    AoM_Ankle_AS_norm_ = [AD.AoM_Ankle_AS_norm; DC.AoM_Ankle_AS_norm; HC.AoM_Ankle_L_norm]  
    AoM_Ankle_US_norm_ = [AD.AoM_Ankle_US_norm; DC.AoM_Ankle_US_norm; HC.AoM_Ankle_R_norm]  

    % General Features
    % Gyro mean
    SC_Gyr_x_mean_ = [AD.SC_Gyr_x_mean; DC.SC_Gyr_x_mean; HC.SC_Gyr_x_mean];
    SC_Gyr_y_mean_ = [AD.SC_Gyr_y_mean; DC.SC_Gyr_y_mean; HC.SC_Gyr_y_mean];
    SC_Gyr_z_mean_ = [AD.SC_Gyr_z_mean; DC.SC_Gyr_z_mean; HC.SC_Gyr_z_mean];
    SC_Gyr_norm_mean_ = [AD.SC_Gyr_norm_mean; DC.SC_Gyr_norm_mean; HC.SC_Gyr_norm_mean];
    DLS_R_Gyr_x_mean_ = [AD.DLS_R_Gyr_x_mean; DC.DLS_R_Gyr_x_mean; HC.DLS_R_Gyr_x_mean];
    DLS_R_Gyr_y_mean_ = [AD.DLS_R_Gyr_y_mean; DC.DLS_R_Gyr_y_mean; HC.DLS_R_Gyr_y_mean];
    DLS_R_Gyr_z_mean_ = [AD.DLS_R_Gyr_z_mean; DC.DLS_R_Gyr_z_mean; HC.DLS_R_Gyr_z_mean];
    DLS_R_Gyr_norm_mean_ = [AD.DLS_R_Gyr_norm_mean; DC.DLS_R_Gyr_norm_mean; HC.DLS_R_Gyr_norm_mean];
    DLS_L_Gyr_x_mean_ = [AD.DLS_L_Gyr_x_mean; DC.DLS_L_Gyr_x_mean; HC.DLS_L_Gyr_x_mean];
    DLS_L_Gyr_y_mean_ = [AD.DLS_L_Gyr_y_mean; DC.DLS_L_Gyr_y_mean; HC.DLS_L_Gyr_y_mean];
    DLS_L_Gyr_z_mean_ = [AD.DLS_L_Gyr_z_mean; DC.DLS_L_Gyr_z_mean; HC.DLS_L_Gyr_z_mean];
    DLS_L_Gyr_norm_mean_ = [AD.DLS_L_Gyr_norm_mean; DC.DLS_L_Gyr_norm_mean; HC.DLS_L_Gyr_norm_mean];
    
    % Range
    SC_Gyr_x_range_ = [AD.SC_Gyr_x_range; DC.SC_Gyr_x_range; HC.SC_Gyr_x_range];
    SC_Gyr_y_range_ = [AD.SC_Gyr_y_range; DC.SC_Gyr_y_range; HC.SC_Gyr_y_range];
    SC_Gyr_z_range_ = [AD.SC_Gyr_z_range; DC.SC_Gyr_z_range; HC.SC_Gyr_z_range];
    SC_Gyr_norm_range_ = [AD.SC_Gyr_norm_range; DC.SC_Gyr_norm_range; HC.SC_Gyr_norm_range];
    DLS_R_Gyr_x_range_ = [AD.DLS_R_Gyr_x_range; DC.DLS_R_Gyr_x_range; HC.DLS_R_Gyr_x_range];
    DLS_R_Gyr_y_range_ = [AD.DLS_R_Gyr_y_range; DC.DLS_R_Gyr_y_range; HC.DLS_R_Gyr_y_range];
    DLS_R_Gyr_z_range_ = [AD.DLS_R_Gyr_z_range; DC.DLS_R_Gyr_z_range; HC.DLS_R_Gyr_z_range];
    DLS_R_Gyr_norm_range_ = [AD.DLS_R_Gyr_norm_range; DC.DLS_R_Gyr_norm_range; HC.DLS_R_Gyr_norm_range];
    DLS_L_Gyr_x_range_ = [AD.DLS_L_Gyr_x_range; DC.DLS_L_Gyr_x_range; HC.DLS_L_Gyr_x_range];
    DLS_L_Gyr_y_range_ = [AD.DLS_L_Gyr_y_range; DC.DLS_L_Gyr_y_range; HC.DLS_L_Gyr_y_range];
    DLS_L_Gyr_z_range_ = [AD.DLS_L_Gyr_z_range; DC.DLS_L_Gyr_z_range; HC.DLS_L_Gyr_z_range];
    DLS_L_Gyr_norm_range_ = [AD.DLS_L_Gyr_norm_range; DC.DLS_L_Gyr_norm_range; HC.DLS_L_Gyr_norm_range];

    % RMS
    SC_Gyr_x_rms_ = [AD.SC_Gyr_x_rms; DC.SC_Gyr_x_rms; HC.SC_Gyr_x_rms];
    SC_Gyr_y_rms_ = [AD.SC_Gyr_y_rms; DC.SC_Gyr_y_rms; HC.SC_Gyr_y_rms];
    SC_Gyr_z_rms_ = [AD.SC_Gyr_z_rms; DC.SC_Gyr_z_rms; HC.SC_Gyr_z_rms];
    SC_Gyr_norm_rms_ = [AD.SC_Gyr_norm_rms; DC.SC_Gyr_norm_rms; HC.SC_Gyr_norm_rms];
    DLS_R_Gyr_x_rms_ = [AD.DLS_R_Gyr_x_rms; DC.DLS_R_Gyr_x_rms; HC.DLS_R_Gyr_x_rms];
    DLS_R_Gyr_y_rms_ = [AD.DLS_R_Gyr_y_rms; DC.DLS_R_Gyr_y_rms; HC.DLS_R_Gyr_y_rms];
    DLS_R_Gyr_z_rms_ = [AD.DLS_R_Gyr_z_rms; DC.DLS_R_Gyr_z_rms; HC.DLS_R_Gyr_z_rms];
    DLS_R_Gyr_norm_rms_ = [AD.DLS_R_Gyr_norm_rms; DC.DLS_R_Gyr_norm_rms; HC.DLS_R_Gyr_norm_rms];
    DLS_L_Gyr_x_rms_ = [AD.DLS_L_Gyr_x_rms; DC.DLS_L_Gyr_x_rms; HC.DLS_L_Gyr_x_rms];
    DLS_L_Gyr_y_rms_ = [AD.DLS_L_Gyr_y_rms; DC.DLS_L_Gyr_y_rms; HC.DLS_L_Gyr_y_rms];
    DLS_L_Gyr_z_rms_ = [AD.DLS_L_Gyr_z_rms; DC.DLS_L_Gyr_z_rms; HC.DLS_L_Gyr_z_rms];
    DLS_L_Gyr_norm_rms_ = [AD.DLS_L_Gyr_norm_rms; DC.DLS_L_Gyr_norm_rms; HC.DLS_L_Gyr_norm_rms];

    % Standard Deviation
    SC_Gyr_x_std_ = [AD.SC_Gyr_x_std; DC.SC_Gyr_x_std; HC.SC_Gyr_x_std];
    SC_Gyr_y_std_ = [AD.SC_Gyr_y_std; DC.SC_Gyr_y_std; HC.SC_Gyr_y_std];
    SC_Gyr_z_std_ = [AD.SC_Gyr_z_std; DC.SC_Gyr_z_std; HC.SC_Gyr_z_std];
    SC_Gyr_norm_std_ = [AD.SC_Gyr_norm_std; DC.SC_Gyr_norm_std; HC.SC_Gyr_norm_std];
    DLS_R_Gyr_x_std_ = [AD.DLS_R_Gyr_x_std; DC.DLS_R_Gyr_x_std; HC.DLS_R_Gyr_x_std];
    DLS_R_Gyr_y_std_ = [AD.DLS_R_Gyr_y_std; DC.DLS_R_Gyr_y_std; HC.DLS_R_Gyr_y_std];
    DLS_R_Gyr_z_std_ = [AD.DLS_R_Gyr_z_std; DC.DLS_R_Gyr_z_std; HC.DLS_R_Gyr_z_std];
    DLS_R_Gyr_norm_std_ = [AD.DLS_R_Gyr_norm_std; DC.DLS_R_Gyr_norm_std; HC.DLS_R_Gyr_norm_std];
    DLS_L_Gyr_x_std_ = [AD.DLS_L_Gyr_x_std; DC.DLS_L_Gyr_x_std; HC.DLS_L_Gyr_x_std];
    DLS_L_Gyr_y_std_ = [AD.DLS_L_Gyr_y_std; DC.DLS_L_Gyr_y_std; HC.DLS_L_Gyr_y_std];
    DLS_L_Gyr_z_std_ = [AD.DLS_L_Gyr_z_std; DC.DLS_L_Gyr_z_std; HC.DLS_L_Gyr_z_std];
    DLS_L_Gyr_norm_std_ = [AD.DLS_L_Gyr_norm_std; DC.DLS_L_Gyr_norm_std; HC.DLS_L_Gyr_norm_std];

    % Skew
    SC_Gyr_x_skew_ = [AD.SC_Gyr_x_skew; DC.SC_Gyr_x_skew; HC.SC_Gyr_x_skew];
    SC_Gyr_y_skew_ = [AD.SC_Gyr_y_skew; DC.SC_Gyr_y_skew; HC.SC_Gyr_y_skew];
    SC_Gyr_z_skew_ = [AD.SC_Gyr_z_skew; DC.SC_Gyr_z_skew; HC.SC_Gyr_z_skew];
    SC_Gyr_norm_skew_ = [AD.SC_Gyr_norm_skew; DC.SC_Gyr_norm_skew; HC.SC_Gyr_norm_skew];
    DLS_R_Gyr_x_skew_ = [AD.DLS_R_Gyr_x_skew; DC.DLS_R_Gyr_x_skew; HC.DLS_R_Gyr_x_skew];
    DLS_R_Gyr_y_skew_ = [AD.DLS_R_Gyr_y_skew; DC.DLS_R_Gyr_y_skew; HC.DLS_R_Gyr_y_skew];
    DLS_R_Gyr_z_skew_ = [AD.DLS_R_Gyr_z_skew; DC.DLS_R_Gyr_z_skew; HC.DLS_R_Gyr_z_skew];
    DLS_R_Gyr_norm_skew_ = [AD.DLS_R_Gyr_norm_skew; DC.DLS_R_Gyr_norm_skew; HC.DLS_R_Gyr_norm_skew];
    DLS_L_Gyr_x_skew_ = [AD.DLS_L_Gyr_x_skew; DC.DLS_L_Gyr_x_skew; HC.DLS_L_Gyr_x_skew];
    DLS_L_Gyr_y_skew_ = [AD.DLS_L_Gyr_y_skew; DC.DLS_L_Gyr_y_skew; HC.DLS_L_Gyr_y_skew];
    DLS_L_Gyr_z_skew_ = [AD.DLS_L_Gyr_z_skew; DC.DLS_L_Gyr_z_skew; HC.DLS_L_Gyr_z_skew];
    DLS_L_Gyr_norm_skew_ = [AD.DLS_L_Gyr_norm_skew; DC.DLS_L_Gyr_norm_skew; HC.DLS_L_Gyr_norm_skew];

    % Kurtosis
    SC_Gyr_x_kurtosis_ = [AD.SC_Gyr_x_kurtosis; DC.SC_Gyr_x_kurtosis; HC.SC_Gyr_x_kurtosis];
    SC_Gyr_y_kurtosis_ = [AD.SC_Gyr_y_kurtosis; DC.SC_Gyr_y_kurtosis; HC.SC_Gyr_y_kurtosis];
    SC_Gyr_z_kurtosis_ = [AD.SC_Gyr_z_kurtosis; DC.SC_Gyr_z_kurtosis; HC.SC_Gyr_z_kurtosis];
    SC_Gyr_norm_kurtosis_ = [AD.SC_Gyr_norm_kurtosis; DC.SC_Gyr_norm_kurtosis; HC.SC_Gyr_norm_kurtosis];
    DLS_R_Gyr_x_kurtosis_ = [AD.DLS_R_Gyr_x_kurtosis; DC.DLS_R_Gyr_x_kurtosis; HC.DLS_R_Gyr_x_kurtosis];
    DLS_R_Gyr_y_kurtosis_ = [AD.DLS_R_Gyr_y_kurtosis; DC.DLS_R_Gyr_y_kurtosis; HC.DLS_R_Gyr_y_kurtosis];
    DLS_R_Gyr_z_kurtosis_ = [AD.DLS_R_Gyr_z_kurtosis; DC.DLS_R_Gyr_z_kurtosis; HC.DLS_R_Gyr_z_kurtosis];
    DLS_R_Gyr_norm_kurtosis_ = [AD.DLS_R_Gyr_norm_kurtosis; DC.DLS_R_Gyr_norm_kurtosis; HC.DLS_R_Gyr_norm_kurtosis];
    DLS_L_Gyr_x_kurtosis_ = [AD.DLS_L_Gyr_x_kurtosis; DC.DLS_L_Gyr_x_kurtosis; HC.DLS_L_Gyr_x_kurtosis];
    DLS_L_Gyr_y_kurtosis_ = [AD.DLS_L_Gyr_y_kurtosis; DC.DLS_L_Gyr_y_kurtosis; HC.DLS_L_Gyr_y_kurtosis];
    DLS_L_Gyr_z_kurtosis_ = [AD.DLS_L_Gyr_z_kurtosis; DC.DLS_L_Gyr_z_kurtosis; HC.DLS_L_Gyr_z_kurtosis];
    DLS_L_Gyr_norm_kurtosis_ = [AD.DLS_L_Gyr_norm_kurtosis; DC.DLS_L_Gyr_norm_kurtosis; HC.DLS_L_Gyr_norm_kurtosis];

    % Pearson correlation coefficient
    SC_Gyr_corr_xy_ = [AD.SC_Gyr_corr_xy; DC.SC_Gyr_corr_xy; HC.SC_Gyr_corr_xy];
    SC_Gyr_corr_xz_ = [AD.SC_Gyr_corr_xz; DC.SC_Gyr_corr_xz; HC.SC_Gyr_corr_xz];
    SC_Gyr_corr_yz_ = [AD.SC_Gyr_corr_yz; DC.SC_Gyr_corr_yz; HC.SC_Gyr_corr_yz];
    DLS_R_Gyr_corr_xy_ = [AD.DLS_R_Gyr_corr_xy; DC.DLS_R_Gyr_corr_xy; HC.DLS_R_Gyr_corr_xy];
    DLS_R_Gyr_corr_xz_ = [AD.DLS_R_Gyr_corr_xz; DC.DLS_R_Gyr_corr_xz; HC.DLS_R_Gyr_corr_xz];
    DLS_R_Gyr_corr_yz_ = [AD.DLS_R_Gyr_corr_yz; DC.DLS_R_Gyr_corr_yz; HC.DLS_R_Gyr_corr_yz];
    DLS_L_Gyr_corr_xy_ = [AD.DLS_L_Gyr_corr_xy; DC.DLS_L_Gyr_corr_xy; HC.DLS_L_Gyr_corr_xy];
    DLS_L_Gyr_corr_xz_ = [AD.DLS_L_Gyr_corr_xz; DC.DLS_L_Gyr_corr_xz; HC.DLS_L_Gyr_corr_xz];
    DLS_L_Gyr_corr_yz_ = [AD.DLS_L_Gyr_corr_yz; DC.DLS_L_Gyr_corr_yz; HC.DLS_L_Gyr_corr_yz];

    % Sample Entropy
    SC_Gyr_x_SamEn_ = [AD.SC_Gyr_x_SamEn; DC.SC_Gyr_x_SamEn; HC.SC_Gyr_x_SamEn];
    SC_Gyr_y_SamEn_ = [AD.SC_Gyr_y_SamEn; DC.SC_Gyr_y_SamEn; HC.SC_Gyr_y_SamEn];
    SC_Gyr_z_SamEn_ = [AD.SC_Gyr_z_SamEn; DC.SC_Gyr_z_SamEn; HC.SC_Gyr_z_SamEn];
    SC_Gyr_norm_SamEn_ = [AD.SC_Gyr_norm_SamEn; DC.SC_Gyr_norm_SamEn; HC.SC_Gyr_norm_SamEn];
    DLS_R_Gyr_x_SamEn_ = [AD.DLS_R_Gyr_x_SamEn; DC.DLS_R_Gyr_x_SamEn; HC.DLS_R_Gyr_x_SamEn];
    DLS_R_Gyr_y_SamEn_ = [AD.DLS_R_Gyr_y_SamEn; DC.DLS_R_Gyr_y_SamEn; HC.DLS_R_Gyr_y_SamEn];
    DLS_R_Gyr_z_SamEn_ = [AD.DLS_R_Gyr_z_SamEn; DC.DLS_R_Gyr_z_SamEn; HC.DLS_R_Gyr_z_SamEn];
    DLS_R_Gyr_norm_SamEn_ = [AD.DLS_R_Gyr_norm_SamEn; DC.DLS_R_Gyr_norm_SamEn; HC.DLS_R_Gyr_norm_SamEn];
    DLS_L_Gyr_x_SamEn_ = [AD.DLS_L_Gyr_x_SamEn; DC.DLS_L_Gyr_x_SamEn; HC.DLS_L_Gyr_x_SamEn];
    DLS_L_Gyr_y_SamEn_ = [AD.DLS_L_Gyr_y_SamEn; DC.DLS_L_Gyr_y_SamEn; HC.DLS_L_Gyr_y_SamEn];
    DLS_L_Gyr_z_SamEn_ = [AD.DLS_L_Gyr_z_SamEn; DC.DLS_L_Gyr_z_SamEn; HC.DLS_L_Gyr_z_SamEn];
    DLS_L_Gyr_norm_SamEn_ = [AD.DLS_L_Gyr_norm_SamEn; DC.DLS_L_Gyr_norm_SamEn; HC.DLS_L_Gyr_norm_SamEn];

    % Frequency domain
    SC_Gyr_x_DAmp_ = [AD.SC_Gyr_x_DAmp; DC.SC_Gyr_x_DAmp; HC.SC_Gyr_x_DAmp];
    SC_Gyr_x_DFreq_ = [AD.SC_Gyr_x_DFreq; DC.SC_Gyr_x_DFreq; HC.SC_Gyr_x_DFreq];
    SC_Gyr_x_PSD_mean_ = [AD.SC_Gyr_x_PSD_mean; DC.SC_Gyr_x_PSD_mean; HC.SC_Gyr_x_PSD_mean];
    SC_Gyr_x_PSD_std_ = [AD.SC_Gyr_x_PSD_std; DC.SC_Gyr_x_PSD_std; HC.SC_Gyr_x_PSD_std];
    SC_Gyr_x_PSD_skew_ = [AD.SC_Gyr_x_PSD_skew; DC.SC_Gyr_x_PSD_skew; HC.SC_Gyr_x_PSD_skew];
    SC_Gyr_x_PSD_kurtosis_ = [AD.SC_Gyr_x_PSD_kurtosis; DC.SC_Gyr_x_PSD_kurtosis; HC.SC_Gyr_x_PSD_kurtosis];

    SC_Gyr_y_DAmp_ = [AD.SC_Gyr_y_DAmp; DC.SC_Gyr_y_DAmp; HC.SC_Gyr_y_DAmp];
    SC_Gyr_y_DFreq_ = [AD.SC_Gyr_y_DFreq; DC.SC_Gyr_y_DFreq; HC.SC_Gyr_y_DFreq];
    SC_Gyr_y_PSD_mean_ = [AD.SC_Gyr_y_PSD_mean; DC.SC_Gyr_y_PSD_mean; HC.SC_Gyr_y_PSD_mean];
    SC_Gyr_y_PSD_std_ = [AD.SC_Gyr_y_PSD_std; DC.SC_Gyr_y_PSD_std; HC.SC_Gyr_y_PSD_std];
    SC_Gyr_y_PSD_skew_ = [AD.SC_Gyr_y_PSD_skew; DC.SC_Gyr_y_PSD_skew; HC.SC_Gyr_y_PSD_skew];
    SC_Gyr_y_PSD_kurtosis_ = [AD.SC_Gyr_y_PSD_kurtosis; DC.SC_Gyr_y_PSD_kurtosis; HC.SC_Gyr_y_PSD_kurtosis];

    SC_Gyr_z_DAmp_ = [AD.SC_Gyr_z_DAmp; DC.SC_Gyr_z_DAmp; HC.SC_Gyr_z_DAmp];
    SC_Gyr_z_DFreq_ = [AD.SC_Gyr_z_DFreq; DC.SC_Gyr_z_DFreq; HC.SC_Gyr_z_DFreq];
    SC_Gyr_z_PSD_mean_ = [AD.SC_Gyr_z_PSD_mean; DC.SC_Gyr_z_PSD_mean; HC.SC_Gyr_z_PSD_mean];
    SC_Gyr_z_PSD_std_ = [AD.SC_Gyr_z_PSD_std; DC.SC_Gyr_z_PSD_std; HC.SC_Gyr_z_PSD_std];
    SC_Gyr_z_PSD_skew_ = [AD.SC_Gyr_z_PSD_skew; DC.SC_Gyr_z_PSD_skew; HC.SC_Gyr_z_PSD_skew];
    SC_Gyr_z_PSD_kurtosis_ = [AD.SC_Gyr_z_PSD_kurtosis; DC.SC_Gyr_z_PSD_kurtosis; HC.SC_Gyr_z_PSD_kurtosis];
    
    SC_Gyr_norm_DAmp_ = [AD.SC_Gyr_norm_DAmp; DC.SC_Gyr_norm_DAmp; HC.SC_Gyr_norm_DAmp];
    SC_Gyr_norm_DFreq_ = [AD.SC_Gyr_norm_DFreq; DC.SC_Gyr_norm_DFreq; HC.SC_Gyr_norm_DFreq];
    SC_Gyr_norm_PSD_mean_ = [AD.SC_Gyr_norm_PSD_mean; DC.SC_Gyr_norm_PSD_mean; HC.SC_Gyr_norm_PSD_mean];
    SC_Gyr_norm_PSD_std_ = [AD.SC_Gyr_norm_PSD_std; DC.SC_Gyr_norm_PSD_std; HC.SC_Gyr_norm_PSD_std];
    SC_Gyr_norm_PSD_skew_ = [AD.SC_Gyr_norm_PSD_skew; DC.SC_Gyr_norm_PSD_skew; HC.SC_Gyr_norm_PSD_skew];
    SC_Gyr_norm_PSD_kurtosis_ = [AD.SC_Gyr_norm_PSD_kurtosis; DC.SC_Gyr_norm_PSD_kurtosis; HC.SC_Gyr_norm_PSD_kurtosis];

    DLS_R_Gyr_x_DAmp_ = [AD.DLS_R_Gyr_x_DAmp; DC.DLS_R_Gyr_x_DAmp; HC.DLS_R_Gyr_x_DAmp];
    DLS_R_Gyr_x_DFreq_ = [AD.DLS_R_Gyr_x_DFreq; DC.DLS_R_Gyr_x_DFreq; HC.DLS_R_Gyr_x_DFreq];
    DLS_R_Gyr_x_PSD_mean_ = [AD.DLS_R_Gyr_x_PSD_mean; DC.DLS_R_Gyr_x_PSD_mean; HC.DLS_R_Gyr_x_PSD_mean];
    DLS_R_Gyr_x_PSD_std_ = [AD.DLS_R_Gyr_x_PSD_std; DC.DLS_R_Gyr_x_PSD_std; HC.DLS_R_Gyr_x_PSD_std];
    DLS_R_Gyr_x_PSD_skew_ = [AD.DLS_R_Gyr_x_PSD_skew; DC.DLS_R_Gyr_x_PSD_skew; HC.DLS_R_Gyr_x_PSD_skew];
    DLS_R_Gyr_x_PSD_kurtosis_ = [AD.DLS_R_Gyr_x_PSD_kurtosis; DC.DLS_R_Gyr_x_PSD_kurtosis; HC.DLS_R_Gyr_x_PSD_kurtosis];

    DLS_R_Gyr_y_DAmp_ = [AD.DLS_R_Gyr_y_DAmp; DC.DLS_R_Gyr_y_DAmp; HC.DLS_R_Gyr_y_DAmp];
    DLS_R_Gyr_y_DFreq_ = [AD.DLS_R_Gyr_y_DFreq; DC.DLS_R_Gyr_y_DFreq; HC.DLS_R_Gyr_y_DFreq];
    DLS_R_Gyr_y_PSD_mean_ = [AD.DLS_R_Gyr_y_PSD_mean; DC.DLS_R_Gyr_y_PSD_mean; HC.DLS_R_Gyr_y_PSD_mean];
    DLS_R_Gyr_y_PSD_std_ = [AD.DLS_R_Gyr_y_PSD_std; DC.DLS_R_Gyr_y_PSD_std; HC.DLS_R_Gyr_y_PSD_std];
    DLS_R_Gyr_y_PSD_skew_ = [AD.DLS_R_Gyr_y_PSD_skew; DC.DLS_R_Gyr_y_PSD_skew; HC.DLS_R_Gyr_y_PSD_skew];
    DLS_R_Gyr_y_PSD_kurtosis_ = [AD.DLS_R_Gyr_y_PSD_kurtosis; DC.DLS_R_Gyr_y_PSD_kurtosis; HC.DLS_R_Gyr_y_PSD_kurtosis];

    DLS_R_Gyr_z_DAmp_ = [AD.DLS_R_Gyr_z_DAmp; DC.DLS_R_Gyr_z_DAmp; HC.DLS_R_Gyr_z_DAmp];
    DLS_R_Gyr_z_DFreq_ = [AD.DLS_R_Gyr_z_DFreq; DC.DLS_R_Gyr_z_DFreq; HC.DLS_R_Gyr_z_DFreq];
    DLS_R_Gyr_z_PSD_mean_ = [AD.DLS_R_Gyr_z_PSD_mean; DC.DLS_R_Gyr_z_PSD_mean; HC.DLS_R_Gyr_z_PSD_mean];
    DLS_R_Gyr_z_PSD_std_ = [AD.DLS_R_Gyr_z_PSD_std; DC.DLS_R_Gyr_z_PSD_std; HC.DLS_R_Gyr_z_PSD_std];
    DLS_R_Gyr_z_PSD_skew_ = [AD.DLS_R_Gyr_z_PSD_skew; DC.DLS_R_Gyr_z_PSD_skew; HC.DLS_R_Gyr_z_PSD_skew];
    DLS_R_Gyr_z_PSD_kurtosis_ = [AD.DLS_R_Gyr_z_PSD_kurtosis; DC.DLS_R_Gyr_z_PSD_kurtosis; HC.DLS_R_Gyr_z_PSD_kurtosis];

    DLS_R_Gyr_norm_DAmp_ = [AD.DLS_R_Gyr_norm_DAmp; DC.DLS_R_Gyr_norm_DAmp; HC.DLS_R_Gyr_norm_DAmp];
    DLS_R_Gyr_norm_DFreq_ = [AD.DLS_R_Gyr_norm_DFreq; DC.DLS_R_Gyr_norm_DFreq; HC.DLS_R_Gyr_norm_DFreq];
    DLS_R_Gyr_norm_PSD_mean_ = [AD.DLS_R_Gyr_norm_PSD_mean; DC.DLS_R_Gyr_norm_PSD_mean; HC.DLS_R_Gyr_norm_PSD_mean];
    DLS_R_Gyr_norm_PSD_std_ = [AD.DLS_R_Gyr_norm_PSD_std; DC.DLS_R_Gyr_norm_PSD_std; HC.DLS_R_Gyr_norm_PSD_std];
    DLS_R_Gyr_norm_PSD_skew_ = [AD.DLS_R_Gyr_norm_PSD_skew; DC.DLS_R_Gyr_norm_PSD_skew; HC.DLS_R_Gyr_norm_PSD_skew];
    DLS_R_Gyr_norm_PSD_kurtosis_ = [AD.DLS_R_Gyr_norm_PSD_kurtosis; DC.DLS_R_Gyr_norm_PSD_kurtosis; HC.DLS_R_Gyr_norm_PSD_kurtosis];
    
    DLS_L_Gyr_x_DAmp_ = [AD.DLS_L_Gyr_x_DAmp; DC.DLS_L_Gyr_x_DAmp; HC.DLS_L_Gyr_x_DAmp];
    DLS_L_Gyr_x_DFreq_ = [AD.DLS_L_Gyr_x_DFreq; DC.DLS_L_Gyr_x_DFreq; HC.DLS_L_Gyr_x_DFreq];
    DLS_L_Gyr_x_PSD_mean_ = [AD.DLS_L_Gyr_x_PSD_mean; DC.DLS_L_Gyr_x_PSD_mean; HC.DLS_L_Gyr_x_PSD_mean];
    DLS_L_Gyr_x_PSD_std_ = [AD.DLS_L_Gyr_x_PSD_std; DC.DLS_L_Gyr_x_PSD_std; HC.DLS_L_Gyr_x_PSD_std];
    DLS_L_Gyr_x_PSD_skew_ = [AD.DLS_L_Gyr_x_PSD_skew; DC.DLS_L_Gyr_x_PSD_skew; HC.DLS_L_Gyr_x_PSD_skew];
    DLS_L_Gyr_x_PSD_kurtosis_ = [AD.DLS_L_Gyr_x_PSD_kurtosis; DC.DLS_L_Gyr_x_PSD_kurtosis; HC.DLS_L_Gyr_x_PSD_kurtosis];

    DLS_L_Gyr_y_DAmp_ = [AD.DLS_L_Gyr_y_DAmp; DC.DLS_L_Gyr_y_DAmp; HC.DLS_L_Gyr_y_DAmp];
    DLS_L_Gyr_y_DFreq_ = [AD.DLS_L_Gyr_y_DFreq; DC.DLS_L_Gyr_y_DFreq; HC.DLS_L_Gyr_y_DFreq];
    DLS_L_Gyr_y_PSD_mean_ = [AD.DLS_L_Gyr_y_PSD_mean; DC.DLS_L_Gyr_y_PSD_mean; HC.DLS_L_Gyr_y_PSD_mean];
    DLS_L_Gyr_y_PSD_std_ = [AD.DLS_L_Gyr_y_PSD_std; DC.DLS_L_Gyr_y_PSD_std; HC.DLS_L_Gyr_y_PSD_std];
    DLS_L_Gyr_y_PSD_skew_ = [AD.DLS_L_Gyr_y_PSD_skew; DC.DLS_L_Gyr_y_PSD_skew; HC.DLS_L_Gyr_y_PSD_skew];
    DLS_L_Gyr_y_PSD_kurtosis_ = [AD.DLS_L_Gyr_y_PSD_kurtosis; DC.DLS_L_Gyr_y_PSD_kurtosis; HC.DLS_L_Gyr_y_PSD_kurtosis];

    DLS_L_Gyr_z_DAmp_ = [AD.DLS_L_Gyr_z_DAmp; DC.DLS_L_Gyr_z_DAmp; HC.DLS_L_Gyr_z_DAmp];
    DLS_L_Gyr_z_DFreq_ = [AD.DLS_L_Gyr_z_DFreq; DC.DLS_L_Gyr_z_DFreq; HC.DLS_L_Gyr_z_DFreq];
    DLS_L_Gyr_z_PSD_mean_ = [AD.DLS_L_Gyr_z_PSD_mean; DC.DLS_L_Gyr_z_PSD_mean; HC.DLS_L_Gyr_z_PSD_mean];
    DLS_L_Gyr_z_PSD_std_ = [AD.DLS_L_Gyr_z_PSD_std; DC.DLS_L_Gyr_z_PSD_std; HC.DLS_L_Gyr_z_PSD_std];
    DLS_L_Gyr_z_PSD_skew_ = [AD.DLS_L_Gyr_z_PSD_skew; DC.DLS_L_Gyr_z_PSD_skew; HC.DLS_L_Gyr_z_PSD_skew];
    DLS_L_Gyr_z_PSD_kurtosis_ = [AD.DLS_L_Gyr_z_PSD_kurtosis; DC.DLS_L_Gyr_z_PSD_kurtosis; HC.DLS_L_Gyr_z_PSD_kurtosis];

    DLS_L_Gyr_norm_DAmp_ = [AD.DLS_L_Gyr_norm_DAmp; DC.DLS_L_Gyr_norm_DAmp; HC.DLS_L_Gyr_norm_DAmp];
    DLS_L_Gyr_norm_DFreq_ = [AD.DLS_L_Gyr_norm_DFreq; DC.DLS_L_Gyr_norm_DFreq; HC.DLS_L_Gyr_norm_DFreq];
    DLS_L_Gyr_norm_PSD_mean_ = [AD.DLS_L_Gyr_norm_PSD_mean; DC.DLS_L_Gyr_norm_PSD_mean; HC.DLS_L_Gyr_norm_PSD_mean];
    DLS_L_Gyr_norm_PSD_std_ = [AD.DLS_L_Gyr_norm_PSD_std; DC.DLS_L_Gyr_norm_PSD_std; HC.DLS_L_Gyr_norm_PSD_std];
    DLS_L_Gyr_norm_PSD_skew_ = [AD.DLS_L_Gyr_norm_PSD_skew; DC.DLS_L_Gyr_norm_PSD_skew; HC.DLS_L_Gyr_norm_PSD_skew];
    DLS_L_Gyr_norm_PSD_kurtosis_ = [AD.DLS_L_Gyr_norm_PSD_kurtosis; DC.DLS_L_Gyr_norm_PSD_kurtosis; HC.DLS_L_Gyr_norm_PSD_kurtosis];



    % Acc mean
    % General Features
    SC_Acc_x_mean_ = [AD.SC_Acc_x_mean; DC.SC_Acc_x_mean; HC.SC_Acc_x_mean];
    SC_Acc_y_mean_ = [AD.SC_Acc_y_mean; DC.SC_Acc_y_mean; HC.SC_Acc_y_mean];
    SC_Acc_z_mean_ = [AD.SC_Acc_z_mean; DC.SC_Acc_z_mean; HC.SC_Acc_z_mean];
    SC_Acc_norm_mean_ = [AD.SC_Acc_norm_mean; DC.SC_Acc_norm_mean; HC.SC_Acc_norm_mean];
    DLS_R_Acc_x_mean_ = [AD.DLS_R_Acc_x_mean; DC.DLS_R_Acc_x_mean; HC.DLS_R_Acc_x_mean];
    DLS_R_Acc_y_mean_ = [AD.DLS_R_Acc_y_mean; DC.DLS_R_Acc_y_mean; HC.DLS_R_Acc_y_mean];
    DLS_R_Acc_z_mean_ = [AD.DLS_R_Acc_z_mean; DC.DLS_R_Acc_z_mean; HC.DLS_R_Acc_z_mean];
    DLS_R_Acc_norm_mean_ = [AD.DLS_R_Acc_norm_mean; DC.DLS_R_Acc_norm_mean; HC.DLS_R_Acc_norm_mean];
    DLS_L_Acc_x_mean_ = [AD.DLS_L_Acc_x_mean; DC.DLS_L_Acc_x_mean; HC.DLS_L_Acc_x_mean];
    DLS_L_Acc_y_mean_ = [AD.DLS_L_Acc_y_mean; DC.DLS_L_Acc_y_mean; HC.DLS_L_Acc_y_mean];
    DLS_L_Acc_z_mean_ = [AD.DLS_L_Acc_z_mean; DC.DLS_L_Acc_z_mean; HC.DLS_L_Acc_z_mean];
    DLS_L_Acc_norm_mean_ = [AD.DLS_L_Acc_norm_mean; DC.DLS_L_Acc_norm_mean; HC.DLS_L_Acc_norm_mean];
    
    % Range
    SC_Acc_x_range_ = [AD.SC_Acc_x_range; DC.SC_Acc_x_range; HC.SC_Acc_x_range];
    SC_Acc_y_range_ = [AD.SC_Acc_y_range; DC.SC_Acc_y_range; HC.SC_Acc_y_range];
    SC_Acc_z_range_ = [AD.SC_Acc_z_range; DC.SC_Acc_z_range; HC.SC_Acc_z_range];
    SC_Acc_norm_range_ = [AD.SC_Acc_norm_range; DC.SC_Acc_norm_range; HC.SC_Acc_norm_range];
    DLS_R_Acc_x_range_ = [AD.DLS_R_Acc_x_range; DC.DLS_R_Acc_x_range; HC.DLS_R_Acc_x_range];
    DLS_R_Acc_y_range_ = [AD.DLS_R_Acc_y_range; DC.DLS_R_Acc_y_range; HC.DLS_R_Acc_y_range];
    DLS_R_Acc_z_range_ = [AD.DLS_R_Acc_z_range; DC.DLS_R_Acc_z_range; HC.DLS_R_Acc_z_range];
    DLS_R_Acc_norm_range_ = [AD.DLS_R_Acc_norm_range; DC.DLS_R_Acc_norm_range; HC.DLS_R_Acc_norm_range];
    DLS_L_Acc_x_range_ = [AD.DLS_L_Acc_x_range; DC.DLS_L_Acc_x_range; HC.DLS_L_Acc_x_range];
    DLS_L_Acc_y_range_ = [AD.DLS_L_Acc_y_range; DC.DLS_L_Acc_y_range; HC.DLS_L_Acc_y_range];
    DLS_L_Acc_z_range_ = [AD.DLS_L_Acc_z_range; DC.DLS_L_Acc_z_range; HC.DLS_L_Acc_z_range];
    DLS_L_Acc_norm_range_ = [AD.DLS_L_Acc_norm_range; DC.DLS_L_Acc_norm_range; HC.DLS_L_Acc_norm_range];

    % RMS
    SC_Acc_x_rms_ = [AD.SC_Acc_x_rms; DC.SC_Acc_x_rms; HC.SC_Acc_x_rms];
    SC_Acc_y_rms_ = [AD.SC_Acc_y_rms; DC.SC_Acc_y_rms; HC.SC_Acc_y_rms];
    SC_Acc_z_rms_ = [AD.SC_Acc_z_rms; DC.SC_Acc_z_rms; HC.SC_Acc_z_rms];
    SC_Acc_norm_rms_ = [AD.SC_Acc_norm_rms; DC.SC_Acc_norm_rms; HC.SC_Acc_norm_rms];
    DLS_R_Acc_x_rms_ = [AD.DLS_R_Acc_x_rms; DC.DLS_R_Acc_x_rms; HC.DLS_R_Acc_x_rms];
    DLS_R_Acc_y_rms_ = [AD.DLS_R_Acc_y_rms; DC.DLS_R_Acc_y_rms; HC.DLS_R_Acc_y_rms];
    DLS_R_Acc_z_rms_ = [AD.DLS_R_Acc_z_rms; DC.DLS_R_Acc_z_rms; HC.DLS_R_Acc_z_rms];
    DLS_R_Acc_norm_rms_ = [AD.DLS_R_Acc_norm_rms; DC.DLS_R_Acc_norm_rms; HC.DLS_R_Acc_norm_rms];
    DLS_L_Acc_x_rms_ = [AD.DLS_L_Acc_x_rms; DC.DLS_L_Acc_x_rms; HC.DLS_L_Acc_x_rms];
    DLS_L_Acc_y_rms_ = [AD.DLS_L_Acc_y_rms; DC.DLS_L_Acc_y_rms; HC.DLS_L_Acc_y_rms];
    DLS_L_Acc_z_rms_ = [AD.DLS_L_Acc_z_rms; DC.DLS_L_Acc_z_rms; HC.DLS_L_Acc_z_rms];
    DLS_L_Acc_norm_rms_ = [AD.DLS_L_Acc_norm_rms; DC.DLS_L_Acc_norm_rms; HC.DLS_L_Acc_norm_rms];

    % Standard Deviation
    SC_Acc_x_std_ = [AD.SC_Acc_x_std; DC.SC_Acc_x_std; HC.SC_Acc_x_std];
    SC_Acc_y_std_ = [AD.SC_Acc_y_std; DC.SC_Acc_y_std; HC.SC_Acc_y_std];
    SC_Acc_z_std_ = [AD.SC_Acc_z_std; DC.SC_Acc_z_std; HC.SC_Acc_z_std];
    SC_Acc_norm_std_ = [AD.SC_Acc_norm_std; DC.SC_Acc_norm_std; HC.SC_Acc_norm_std];
    DLS_R_Acc_x_std_ = [AD.DLS_R_Acc_x_std; DC.DLS_R_Acc_x_std; HC.DLS_R_Acc_x_std];
    DLS_R_Acc_y_std_ = [AD.DLS_R_Acc_y_std; DC.DLS_R_Acc_y_std; HC.DLS_R_Acc_y_std];
    DLS_R_Acc_z_std_ = [AD.DLS_R_Acc_z_std; DC.DLS_R_Acc_z_std; HC.DLS_R_Acc_z_std];
    DLS_R_Acc_norm_std_ = [AD.DLS_R_Acc_norm_std; DC.DLS_R_Acc_norm_std; HC.DLS_R_Acc_norm_std];
    DLS_L_Acc_x_std_ = [AD.DLS_L_Acc_x_std; DC.DLS_L_Acc_x_std; HC.DLS_L_Acc_x_std];
    DLS_L_Acc_y_std_ = [AD.DLS_L_Acc_y_std; DC.DLS_L_Acc_y_std; HC.DLS_L_Acc_y_std];
    DLS_L_Acc_z_std_ = [AD.DLS_L_Acc_z_std; DC.DLS_L_Acc_z_std; HC.DLS_L_Acc_z_std];
    DLS_L_Acc_norm_std_ = [AD.DLS_L_Acc_norm_std; DC.DLS_L_Acc_norm_std; HC.DLS_L_Acc_norm_std];

    % Skew
    SC_Acc_x_skew_ = [AD.SC_Acc_x_skew; DC.SC_Acc_x_skew; HC.SC_Acc_x_skew];
    SC_Acc_y_skew_ = [AD.SC_Acc_y_skew; DC.SC_Acc_y_skew; HC.SC_Acc_y_skew];
    SC_Acc_z_skew_ = [AD.SC_Acc_z_skew; DC.SC_Acc_z_skew; HC.SC_Acc_z_skew];
    SC_Acc_norm_skew_ = [AD.SC_Acc_norm_skew; DC.SC_Acc_norm_skew; HC.SC_Acc_norm_skew];
    DLS_R_Acc_x_skew_ = [AD.DLS_R_Acc_x_skew; DC.DLS_R_Acc_x_skew; HC.DLS_R_Acc_x_skew];
    DLS_R_Acc_y_skew_ = [AD.DLS_R_Acc_y_skew; DC.DLS_R_Acc_y_skew; HC.DLS_R_Acc_y_skew];
    DLS_R_Acc_z_skew_ = [AD.DLS_R_Acc_z_skew; DC.DLS_R_Acc_z_skew; HC.DLS_R_Acc_z_skew];
    DLS_R_Acc_norm_skew_ = [AD.DLS_R_Acc_norm_skew; DC.DLS_R_Acc_norm_skew; HC.DLS_R_Acc_norm_skew];
    DLS_L_Acc_x_skew_ = [AD.DLS_L_Acc_x_skew; DC.DLS_L_Acc_x_skew; HC.DLS_L_Acc_x_skew];
    DLS_L_Acc_y_skew_ = [AD.DLS_L_Acc_y_skew; DC.DLS_L_Acc_y_skew; HC.DLS_L_Acc_y_skew];
    DLS_L_Acc_z_skew_ = [AD.DLS_L_Acc_z_skew; DC.DLS_L_Acc_z_skew; HC.DLS_L_Acc_z_skew];
    DLS_L_Acc_norm_skew_ = [AD.DLS_L_Acc_norm_skew; DC.DLS_L_Acc_norm_skew; HC.DLS_L_Acc_norm_skew];

    % Kurtosis
    SC_Acc_x_kurtosis_ = [AD.SC_Acc_x_kurtosis; DC.SC_Acc_x_kurtosis; HC.SC_Acc_x_kurtosis];
    SC_Acc_y_kurtosis_ = [AD.SC_Acc_y_kurtosis; DC.SC_Acc_y_kurtosis; HC.SC_Acc_y_kurtosis];
    SC_Acc_z_kurtosis_ = [AD.SC_Acc_z_kurtosis; DC.SC_Acc_z_kurtosis; HC.SC_Acc_z_kurtosis];
    SC_Acc_norm_kurtosis_ = [AD.SC_Acc_norm_kurtosis; DC.SC_Acc_norm_kurtosis; HC.SC_Acc_norm_kurtosis];
    DLS_R_Acc_x_kurtosis_ = [AD.DLS_R_Acc_x_kurtosis; DC.DLS_R_Acc_x_kurtosis; HC.DLS_R_Acc_x_kurtosis];
    DLS_R_Acc_y_kurtosis_ = [AD.DLS_R_Acc_y_kurtosis; DC.DLS_R_Acc_y_kurtosis; HC.DLS_R_Acc_y_kurtosis];
    DLS_R_Acc_z_kurtosis_ = [AD.DLS_R_Acc_z_kurtosis; DC.DLS_R_Acc_z_kurtosis; HC.DLS_R_Acc_z_kurtosis];
    DLS_R_Acc_norm_kurtosis_ = [AD.DLS_R_Acc_norm_kurtosis; DC.DLS_R_Acc_norm_kurtosis; HC.DLS_R_Acc_norm_kurtosis];
    DLS_L_Acc_x_kurtosis_ = [AD.DLS_L_Acc_x_kurtosis; DC.DLS_L_Acc_x_kurtosis; HC.DLS_L_Acc_x_kurtosis];
    DLS_L_Acc_y_kurtosis_ = [AD.DLS_L_Acc_y_kurtosis; DC.DLS_L_Acc_y_kurtosis; HC.DLS_L_Acc_y_kurtosis];
    DLS_L_Acc_z_kurtosis_ = [AD.DLS_L_Acc_z_kurtosis; DC.DLS_L_Acc_z_kurtosis; HC.DLS_L_Acc_z_kurtosis];
    DLS_L_Acc_norm_kurtosis_ = [AD.DLS_L_Acc_norm_kurtosis; DC.DLS_L_Acc_norm_kurtosis; HC.DLS_L_Acc_norm_kurtosis];

    % Pearson correlation coefficient
    SC_Acc_corr_xy_ = [AD.SC_Acc_corr_xy; DC.SC_Acc_corr_xy; HC.SC_Acc_corr_xy];
    SC_Acc_corr_xz_ = [AD.SC_Acc_corr_xz; DC.SC_Acc_corr_xz; HC.SC_Acc_corr_xz];
    SC_Acc_corr_yz_ = [AD.SC_Acc_corr_yz; DC.SC_Acc_corr_yz; HC.SC_Acc_corr_yz];
    DLS_R_Acc_corr_xy_ = [AD.DLS_R_Acc_corr_xy; DC.DLS_R_Acc_corr_xy; HC.DLS_R_Acc_corr_xy];
    DLS_R_Acc_corr_xz_ = [AD.DLS_R_Acc_corr_xz; DC.DLS_R_Acc_corr_xz; HC.DLS_R_Acc_corr_xz];
    DLS_R_Acc_corr_yz_ = [AD.DLS_R_Acc_corr_yz; DC.DLS_R_Acc_corr_yz; HC.DLS_R_Acc_corr_yz];
    DLS_L_Acc_corr_xy_ = [AD.DLS_L_Acc_corr_xy; DC.DLS_L_Acc_corr_xy; HC.DLS_L_Acc_corr_xy];
    DLS_L_Acc_corr_xz_ = [AD.DLS_L_Acc_corr_xz; DC.DLS_L_Acc_corr_xz; HC.DLS_L_Acc_corr_xz];
    DLS_L_Acc_corr_yz_ = [AD.DLS_L_Acc_corr_yz; DC.DLS_L_Acc_corr_yz; HC.DLS_L_Acc_corr_yz];

    % Sample Entropy
    SC_Acc_x_SamEn_ = [AD.SC_Acc_x_SamEn; DC.SC_Acc_x_SamEn; HC.SC_Acc_x_SamEn];
    SC_Acc_y_SamEn_ = [AD.SC_Acc_y_SamEn; DC.SC_Acc_y_SamEn; HC.SC_Acc_y_SamEn];
    SC_Acc_z_SamEn_ = [AD.SC_Acc_z_SamEn; DC.SC_Acc_z_SamEn; HC.SC_Acc_z_SamEn];
    SC_Acc_norm_SamEn_ = [AD.SC_Acc_norm_SamEn; DC.SC_Acc_norm_SamEn; HC.SC_Acc_norm_SamEn];
    DLS_R_Acc_x_SamEn_ = [AD.DLS_R_Acc_x_SamEn; DC.DLS_R_Acc_x_SamEn; HC.DLS_R_Acc_x_SamEn];
    DLS_R_Acc_y_SamEn_ = [AD.DLS_R_Acc_y_SamEn; DC.DLS_R_Acc_y_SamEn; HC.DLS_R_Acc_y_SamEn];
    DLS_R_Acc_z_SamEn_ = [AD.DLS_R_Acc_z_SamEn; DC.DLS_R_Acc_z_SamEn; HC.DLS_R_Acc_z_SamEn];
    DLS_R_Acc_norm_SamEn_ = [AD.DLS_R_Acc_norm_SamEn; DC.DLS_R_Acc_norm_SamEn; HC.DLS_R_Acc_norm_SamEn];
    DLS_L_Acc_x_SamEn_ = [AD.DLS_L_Acc_x_SamEn; DC.DLS_L_Acc_x_SamEn; HC.DLS_L_Acc_x_SamEn];
    DLS_L_Acc_y_SamEn_ = [AD.DLS_L_Acc_y_SamEn; DC.DLS_L_Acc_y_SamEn; HC.DLS_L_Acc_y_SamEn];
    DLS_L_Acc_z_SamEn_ = [AD.DLS_L_Acc_z_SamEn; DC.DLS_L_Acc_z_SamEn; HC.DLS_L_Acc_z_SamEn];
    DLS_L_Acc_norm_SamEn_ = [AD.DLS_L_Acc_norm_SamEn; DC.DLS_L_Acc_norm_SamEn; HC.DLS_L_Acc_norm_SamEn];

    % Frequency domain
    SC_Acc_x_DAmp_ = [AD.SC_Acc_x_DAmp; DC.SC_Acc_x_DAmp; HC.SC_Acc_x_DAmp];
    SC_Acc_x_DFreq_ = [AD.SC_Acc_x_DFreq; DC.SC_Acc_x_DFreq; HC.SC_Acc_x_DFreq];
    SC_Acc_x_PSD_mean_ = [AD.SC_Acc_x_PSD_mean; DC.SC_Acc_x_PSD_mean; HC.SC_Acc_x_PSD_mean];
    SC_Acc_x_PSD_std_ = [AD.SC_Acc_x_PSD_std; DC.SC_Acc_x_PSD_std; HC.SC_Acc_x_PSD_std];
    SC_Acc_x_PSD_skew_ = [AD.SC_Acc_x_PSD_skew; DC.SC_Acc_x_PSD_skew; HC.SC_Acc_x_PSD_skew];
    SC_Acc_x_PSD_kurtosis_ = [AD.SC_Acc_x_PSD_kurtosis; DC.SC_Acc_x_PSD_kurtosis; HC.SC_Acc_x_PSD_kurtosis];

    SC_Acc_y_DAmp_ = [AD.SC_Acc_y_DAmp; DC.SC_Acc_y_DAmp; HC.SC_Acc_y_DAmp];
    SC_Acc_y_DFreq_ = [AD.SC_Acc_y_DFreq; DC.SC_Acc_y_DFreq; HC.SC_Acc_y_DFreq];
    SC_Acc_y_PSD_mean_ = [AD.SC_Acc_y_PSD_mean; DC.SC_Acc_y_PSD_mean; HC.SC_Acc_y_PSD_mean];
    SC_Acc_y_PSD_std_ = [AD.SC_Acc_y_PSD_std; DC.SC_Acc_y_PSD_std; HC.SC_Acc_y_PSD_std];
    SC_Acc_y_PSD_skew_ = [AD.SC_Acc_y_PSD_skew; DC.SC_Acc_y_PSD_skew; HC.SC_Acc_y_PSD_skew];
    SC_Acc_y_PSD_kurtosis_ = [AD.SC_Acc_y_PSD_kurtosis; DC.SC_Acc_y_PSD_kurtosis; HC.SC_Acc_y_PSD_kurtosis];

    SC_Acc_z_DAmp_ = [AD.SC_Acc_z_DAmp; DC.SC_Acc_z_DAmp; HC.SC_Acc_z_DAmp];
    SC_Acc_z_DFreq_ = [AD.SC_Acc_z_DFreq; DC.SC_Acc_z_DFreq; HC.SC_Acc_z_DFreq];
    SC_Acc_z_PSD_mean_ = [AD.SC_Acc_z_PSD_mean; DC.SC_Acc_z_PSD_mean; HC.SC_Acc_z_PSD_mean];
    SC_Acc_z_PSD_std_ = [AD.SC_Acc_z_PSD_std; DC.SC_Acc_z_PSD_std; HC.SC_Acc_z_PSD_std];
    SC_Acc_z_PSD_skew_ = [AD.SC_Acc_z_PSD_skew; DC.SC_Acc_z_PSD_skew; HC.SC_Acc_z_PSD_skew];
    SC_Acc_z_PSD_kurtosis_ = [AD.SC_Acc_z_PSD_kurtosis; DC.SC_Acc_z_PSD_kurtosis; HC.SC_Acc_z_PSD_kurtosis];
    
    SC_Acc_norm_DAmp_ = [AD.SC_Acc_norm_DAmp; DC.SC_Acc_norm_DAmp; HC.SC_Acc_norm_DAmp];
    SC_Acc_norm_DFreq_ = [AD.SC_Acc_norm_DFreq; DC.SC_Acc_norm_DFreq; HC.SC_Acc_norm_DFreq];
    SC_Acc_norm_PSD_mean_ = [AD.SC_Acc_norm_PSD_mean; DC.SC_Acc_norm_PSD_mean; HC.SC_Acc_norm_PSD_mean];
    SC_Acc_norm_PSD_std_ = [AD.SC_Acc_norm_PSD_std; DC.SC_Acc_norm_PSD_std; HC.SC_Acc_norm_PSD_std];
    SC_Acc_norm_PSD_skew_ = [AD.SC_Acc_norm_PSD_skew; DC.SC_Acc_norm_PSD_skew; HC.SC_Acc_norm_PSD_skew];
    SC_Acc_norm_PSD_kurtosis_ = [AD.SC_Acc_norm_PSD_kurtosis; DC.SC_Acc_norm_PSD_kurtosis; HC.SC_Acc_norm_PSD_kurtosis];

    DLS_R_Acc_x_DAmp_ = [AD.DLS_R_Acc_x_DAmp; DC.DLS_R_Acc_x_DAmp; HC.DLS_R_Acc_x_DAmp];
    DLS_R_Acc_x_DFreq_ = [AD.DLS_R_Acc_x_DFreq; DC.DLS_R_Acc_x_DFreq; HC.DLS_R_Acc_x_DFreq];
    DLS_R_Acc_x_PSD_mean_ = [AD.DLS_R_Acc_x_PSD_mean; DC.DLS_R_Acc_x_PSD_mean; HC.DLS_R_Acc_x_PSD_mean];
    DLS_R_Acc_x_PSD_std_ = [AD.DLS_R_Acc_x_PSD_std; DC.DLS_R_Acc_x_PSD_std; HC.DLS_R_Acc_x_PSD_std];
    DLS_R_Acc_x_PSD_skew_ = [AD.DLS_R_Acc_x_PSD_skew; DC.DLS_R_Acc_x_PSD_skew; HC.DLS_R_Acc_x_PSD_skew];
    DLS_R_Acc_x_PSD_kurtosis_ = [AD.DLS_R_Acc_x_PSD_kurtosis; DC.DLS_R_Acc_x_PSD_kurtosis; HC.DLS_R_Acc_x_PSD_kurtosis];

    DLS_R_Acc_y_DAmp_ = [AD.DLS_R_Acc_y_DAmp; DC.DLS_R_Acc_y_DAmp; HC.DLS_R_Acc_y_DAmp];
    DLS_R_Acc_y_DFreq_ = [AD.DLS_R_Acc_y_DFreq; DC.DLS_R_Acc_y_DFreq; HC.DLS_R_Acc_y_DFreq];
    DLS_R_Acc_y_PSD_mean_ = [AD.DLS_R_Acc_y_PSD_mean; DC.DLS_R_Acc_y_PSD_mean; HC.DLS_R_Acc_y_PSD_mean];
    DLS_R_Acc_y_PSD_std_ = [AD.DLS_R_Acc_y_PSD_std; DC.DLS_R_Acc_y_PSD_std; HC.DLS_R_Acc_y_PSD_std];
    DLS_R_Acc_y_PSD_skew_ = [AD.DLS_R_Acc_y_PSD_skew; DC.DLS_R_Acc_y_PSD_skew; HC.DLS_R_Acc_y_PSD_skew];
    DLS_R_Acc_y_PSD_kurtosis_ = [AD.DLS_R_Acc_y_PSD_kurtosis; DC.DLS_R_Acc_y_PSD_kurtosis; HC.DLS_R_Acc_y_PSD_kurtosis];

    DLS_R_Acc_z_DAmp_ = [AD.DLS_R_Acc_z_DAmp; DC.DLS_R_Acc_z_DAmp; HC.DLS_R_Acc_z_DAmp];
    DLS_R_Acc_z_DFreq_ = [AD.DLS_R_Acc_z_DFreq; DC.DLS_R_Acc_z_DFreq; HC.DLS_R_Acc_z_DFreq];
    DLS_R_Acc_z_PSD_mean_ = [AD.DLS_R_Acc_z_PSD_mean; DC.DLS_R_Acc_z_PSD_mean; HC.DLS_R_Acc_z_PSD_mean];
    DLS_R_Acc_z_PSD_std_ = [AD.DLS_R_Acc_z_PSD_std; DC.DLS_R_Acc_z_PSD_std; HC.DLS_R_Acc_z_PSD_std];
    DLS_R_Acc_z_PSD_skew_ = [AD.DLS_R_Acc_z_PSD_skew; DC.DLS_R_Acc_z_PSD_skew; HC.DLS_R_Acc_z_PSD_skew];
    DLS_R_Acc_z_PSD_kurtosis_ = [AD.DLS_R_Acc_z_PSD_kurtosis; DC.DLS_R_Acc_z_PSD_kurtosis; HC.DLS_R_Acc_z_PSD_kurtosis];

    DLS_R_Acc_norm_DAmp_ = [AD.DLS_R_Acc_norm_DAmp; DC.DLS_R_Acc_norm_DAmp; HC.DLS_R_Acc_norm_DAmp];
    DLS_R_Acc_norm_DFreq_ = [AD.DLS_R_Acc_norm_DFreq; DC.DLS_R_Acc_norm_DFreq; HC.DLS_R_Acc_norm_DFreq];
    DLS_R_Acc_norm_PSD_mean_ = [AD.DLS_R_Acc_norm_PSD_mean; DC.DLS_R_Acc_norm_PSD_mean; HC.DLS_R_Acc_norm_PSD_mean];
    DLS_R_Acc_norm_PSD_std_ = [AD.DLS_R_Acc_norm_PSD_std; DC.DLS_R_Acc_norm_PSD_std; HC.DLS_R_Acc_norm_PSD_std];
    DLS_R_Acc_norm_PSD_skew_ = [AD.DLS_R_Acc_norm_PSD_skew; DC.DLS_R_Acc_norm_PSD_skew; HC.DLS_R_Acc_norm_PSD_skew];
    DLS_R_Acc_norm_PSD_kurtosis_ = [AD.DLS_R_Acc_norm_PSD_kurtosis; DC.DLS_R_Acc_norm_PSD_kurtosis; HC.DLS_R_Acc_norm_PSD_kurtosis];
    
    DLS_L_Acc_x_DAmp_ = [AD.DLS_L_Acc_x_DAmp; DC.DLS_L_Acc_x_DAmp; HC.DLS_L_Acc_x_DAmp];
    DLS_L_Acc_x_DFreq_ = [AD.DLS_L_Acc_x_DFreq; DC.DLS_L_Acc_x_DFreq; HC.DLS_L_Acc_x_DFreq];
    DLS_L_Acc_x_PSD_mean_ = [AD.DLS_L_Acc_x_PSD_mean; DC.DLS_L_Acc_x_PSD_mean; HC.DLS_L_Acc_x_PSD_mean];
    DLS_L_Acc_x_PSD_std_ = [AD.DLS_L_Acc_x_PSD_std; DC.DLS_L_Acc_x_PSD_std; HC.DLS_L_Acc_x_PSD_std];
    DLS_L_Acc_x_PSD_skew_ = [AD.DLS_L_Acc_x_PSD_skew; DC.DLS_L_Acc_x_PSD_skew; HC.DLS_L_Acc_x_PSD_skew];
    DLS_L_Acc_x_PSD_kurtosis_ = [AD.DLS_L_Acc_x_PSD_kurtosis; DC.DLS_L_Acc_x_PSD_kurtosis; HC.DLS_L_Acc_x_PSD_kurtosis];

    DLS_L_Acc_y_DAmp_ = [AD.DLS_L_Acc_y_DAmp; DC.DLS_L_Acc_y_DAmp; HC.DLS_L_Acc_y_DAmp];
    DLS_L_Acc_y_DFreq_ = [AD.DLS_L_Acc_y_DFreq; DC.DLS_L_Acc_y_DFreq; HC.DLS_L_Acc_y_DFreq];
    DLS_L_Acc_y_PSD_mean_ = [AD.DLS_L_Acc_y_PSD_mean; DC.DLS_L_Acc_y_PSD_mean; HC.DLS_L_Acc_y_PSD_mean];
    DLS_L_Acc_y_PSD_std_ = [AD.DLS_L_Acc_y_PSD_std; DC.DLS_L_Acc_y_PSD_std; HC.DLS_L_Acc_y_PSD_std];
    DLS_L_Acc_y_PSD_skew_ = [AD.DLS_L_Acc_y_PSD_skew; DC.DLS_L_Acc_y_PSD_skew; HC.DLS_L_Acc_y_PSD_skew];
    DLS_L_Acc_y_PSD_kurtosis_ = [AD.DLS_L_Acc_y_PSD_kurtosis; DC.DLS_L_Acc_y_PSD_kurtosis; HC.DLS_L_Acc_y_PSD_kurtosis];

    DLS_L_Acc_z_DAmp_ = [AD.DLS_L_Acc_z_DAmp; DC.DLS_L_Acc_z_DAmp; HC.DLS_L_Acc_z_DAmp];
    DLS_L_Acc_z_DFreq_ = [AD.DLS_L_Acc_z_DFreq; DC.DLS_L_Acc_z_DFreq; HC.DLS_L_Acc_z_DFreq];
    DLS_L_Acc_z_PSD_mean_ = [AD.DLS_L_Acc_z_PSD_mean; DC.DLS_L_Acc_z_PSD_mean; HC.DLS_L_Acc_z_PSD_mean];
    DLS_L_Acc_z_PSD_std_ = [AD.DLS_L_Acc_z_PSD_std; DC.DLS_L_Acc_z_PSD_std; HC.DLS_L_Acc_z_PSD_std];
    DLS_L_Acc_z_PSD_skew_ = [AD.DLS_L_Acc_z_PSD_skew; DC.DLS_L_Acc_z_PSD_skew; HC.DLS_L_Acc_z_PSD_skew];
    DLS_L_Acc_z_PSD_kurtosis_ = [AD.DLS_L_Acc_z_PSD_kurtosis; DC.DLS_L_Acc_z_PSD_kurtosis; HC.DLS_L_Acc_z_PSD_kurtosis];

    DLS_L_Acc_norm_DAmp_ = [AD.DLS_L_Acc_norm_DAmp; DC.DLS_L_Acc_norm_DAmp; HC.DLS_L_Acc_norm_DAmp];
    DLS_L_Acc_norm_DFreq_ = [AD.DLS_L_Acc_norm_DFreq; DC.DLS_L_Acc_norm_DFreq; HC.DLS_L_Acc_norm_DFreq];
    DLS_L_Acc_norm_PSD_mean_ = [AD.DLS_L_Acc_norm_PSD_mean; DC.DLS_L_Acc_norm_PSD_mean; HC.DLS_L_Acc_norm_PSD_mean];
    DLS_L_Acc_norm_PSD_std_ = [AD.DLS_L_Acc_norm_PSD_std; DC.DLS_L_Acc_norm_PSD_std; HC.DLS_L_Acc_norm_PSD_std];
    DLS_L_Acc_norm_PSD_skew_ = [AD.DLS_L_Acc_norm_PSD_skew; DC.DLS_L_Acc_norm_PSD_skew; HC.DLS_L_Acc_norm_PSD_skew];
    DLS_L_Acc_norm_PSD_kurtosis_ = [AD.DLS_L_Acc_norm_PSD_kurtosis; DC.DLS_L_Acc_norm_PSD_kurtosis; HC.DLS_L_Acc_norm_PSD_kurtosis];


 
  
    CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_Final.csv');
    % Clinical Outcomes
    Out_ID = [32 25 7 2 1];

    for i = 1:1:length(Out_ID)
        CO(Out_ID(i),:) = [];
    end

    FIM_DC = CO(:,2);
    BBS_DC = CO(:,3);
    MWT10_SSV_DC = CO(:,4);
    MWT10_FV_DC = CO(:,5);
    MWT6_DC = CO(:,6);
    TUG_DC = CO(:,7);
    FIM_M_DC = CO(:,8);
    
    CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_First.csv');
    % Baseline Clinical Outcomes
    Out_ID = [49 44 43 42 40 35 32 30 25 22 8 7 3 1];

    for i = 1:1:length(Out_ID)
        CO(Out_ID(i),:) = [];
    end

    FIM_AD = CO(:,2);
    BBS_AD = CO(:,3);
    MWT10_SSV_AD = CO(:,4);
    MWT10_FV_AD = CO(:,5);
    MWT6_AD = CO(:,6);
    TUG_AD = CO(:,7);
    FIM_M_AD = CO(:,8);

    CO = import_Clincal_Outcomes_Discharge('./HC_Clinical_Outcome.csv');
    % Clinical Outcomes
    Out_ID = [38 26 4];

    for i = 1:1:length(Out_ID)
        CO(Out_ID(i),:) = [];
    end

    FIM_HC = CO(:,2);
    BBS_HC = CO(:,3);
    MWT10_SSV_HC = CO(:,4);
    MWT10_FV_HC = CO(:,5);
    MWT6_HC = CO(:,6);
    TUG_HC = CO(:,7);
    FIM_M_HC = CO(:,8);


    FIM_ = [FIM_AD; FIM_DC; FIM_HC];
    FIM_M_ = [FIM_M_AD; FIM_M_DC; FIM_M_HC];
    BBS_ = [BBS_AD; BBS_DC; BBS_HC];
    MWT10_SSV_ = [MWT10_SSV_AD; MWT10_SSV_DC; MWT10_SSV_HC];
    MWT10_FV_ = [MWT10_FV_AD; MWT10_FV_DC; MWT10_FV_HC];
    MWT6_ = [MWT6_AD; MWT6_DC; MWT6_HC];
    TUG_ = [TUG_AD; TUG_DC; TUG_HC]; 
    
    
    ID = [ID_AD ID_DC ID_HC]';
    for i = 1:1:length(FIM_)
        if i >= 1 && i <= 41
            Sub_Type_{i} = 'CVA_AD'
            Sub_ID_{i} = ['CVA' num2str(ID(i))];
        elseif i >= 42 && i <= 91
            Sub_Type_{i} = 'CVA_DC'
            Sub_ID_{i} = ['CVA' num2str(ID(i))];
        elseif i >= 92 && i <= 139
            Sub_Type_{i} = 'HC'
            Sub_ID_{i} = ['HC' num2str(ID(i))];
        end
    end
    
    Sub_ID_ = Sub_ID_';
    Sub_Type_ = Sub_Type_';
    Cut_Off_Time_ = TS(k)*ones(139,1);

  
    Sub_ID = [Sub_ID; Sub_ID_];
    Sub_Type = [Sub_Type; Sub_Type_];
    Cut_Off_Time = [Cut_Off_Time; Cut_Off_Time_];
    FIM = [FIM; FIM_];
    FIM_M = [FIM_M; FIM_M_];
    BBS = [BBS; BBS_];
    MWT10_SSV = [MWT10_SSV; MWT10_SSV_];
    MWT10_FV = [MWT10_FV; MWT10_FV_];
    MWT6 = [MWT6; MWT6_];
    TUG = [TUG; TUG_];
    
    clear Sub_Type_ Sub_ID_ Time_

    AoM_Pel_tilt = [AoM_Pel_tilt; AoM_Pel_tilt_];
    AoM_Pel_ro = [AoM_Pel_ro; AoM_Pel_ro_];
    AoM_Pel_oblq = [AoM_Pel_oblq; AoM_Pel_oblq_];
    AoM_Pel_norm = [AoM_Pel_norm; AoM_Pel_norm_];
    AoM_Ankle_AS_x = [AoM_Ankle_AS_x; AoM_Ankle_AS_x_];  
    AoM_Ankle_US_x = [AoM_Ankle_US_x; AoM_Ankle_US_x_]; 
    AoM_Ankle_AS_y = [AoM_Ankle_AS_y; AoM_Ankle_AS_y_]; 
    AoM_Ankle_US_y = [AoM_Ankle_US_y; AoM_Ankle_US_y_]; 
    AoM_Ankle_AS_z = [AoM_Ankle_AS_z; AoM_Ankle_AS_z_];  
    AoM_Ankle_US_z = [AoM_Ankle_US_z; AoM_Ankle_US_z_]; 
    AoM_Ankle_AS_norm = [AoM_Ankle_AS_norm; AoM_Ankle_AS_norm_];  
    AoM_Ankle_US_norm = [AoM_Ankle_US_norm; AoM_Ankle_US_norm_];  
    
    % General Features
    % Gyro mean
    SC_Gyr_x_mean = [SC_Gyr_x_mean; SC_Gyr_x_mean_]
    
    SC_Gyr_y_mean = [SC_Gyr_y_mean; SC_Gyr_y_mean_]
    SC_Gyr_z_mean = [SC_Gyr_z_mean; SC_Gyr_z_mean_]
    SC_Gyr_norm_mean = [SC_Gyr_norm_mean; SC_Gyr_norm_mean_]
    DLS_R_Gyr_x_mean = [DLS_R_Gyr_x_mean; DLS_R_Gyr_x_mean_];
    DLS_R_Gyr_y_mean = [DLS_R_Gyr_y_mean; DLS_R_Gyr_y_mean_];
    DLS_R_Gyr_z_mean = [DLS_R_Gyr_z_mean; DLS_R_Gyr_z_mean_];
    DLS_R_Gyr_norm_mean = [DLS_R_Gyr_norm_mean; DLS_R_Gyr_norm_mean_];
    DLS_L_Gyr_x_mean = [DLS_L_Gyr_x_mean; DLS_L_Gyr_x_mean_];
    DLS_L_Gyr_y_mean = [DLS_L_Gyr_y_mean; DLS_L_Gyr_y_mean_];
    DLS_L_Gyr_z_mean = [DLS_L_Gyr_z_mean; DLS_L_Gyr_z_mean_];
    DLS_L_Gyr_norm_mean = [DLS_L_Gyr_norm_mean; DLS_L_Gyr_norm_mean_];
    
    % Range
    SC_Gyr_x_range = [SC_Gyr_x_range; SC_Gyr_x_range_]
    SC_Gyr_y_range = [SC_Gyr_y_range; SC_Gyr_y_range_]
    SC_Gyr_z_range = [SC_Gyr_z_range; SC_Gyr_z_range_]
    SC_Gyr_norm_range = [SC_Gyr_norm_range; SC_Gyr_norm_range_]
    DLS_R_Gyr_x_range = [DLS_R_Gyr_x_range; DLS_R_Gyr_x_range_];
    DLS_R_Gyr_y_range = [DLS_R_Gyr_y_range; DLS_R_Gyr_y_range_];
    DLS_R_Gyr_z_range = [DLS_R_Gyr_z_range; DLS_R_Gyr_z_range_];
    DLS_R_Gyr_norm_range = [DLS_R_Gyr_norm_range; DLS_R_Gyr_norm_range_];
    DLS_L_Gyr_x_range = [DLS_L_Gyr_x_range; DLS_L_Gyr_x_range_];
    DLS_L_Gyr_y_range = [DLS_L_Gyr_y_range; DLS_L_Gyr_y_range_];
    DLS_L_Gyr_z_range = [DLS_L_Gyr_z_range; DLS_L_Gyr_z_range_];
    DLS_L_Gyr_norm_range = [DLS_L_Gyr_norm_range; DLS_L_Gyr_norm_range_];

    % RMS
    SC_Gyr_x_rms = [SC_Gyr_x_rms; SC_Gyr_x_rms_]
    SC_Gyr_y_rms = [SC_Gyr_y_rms; SC_Gyr_y_rms_]
    SC_Gyr_z_rms = [SC_Gyr_z_rms; SC_Gyr_z_rms_]
    SC_Gyr_norm_rms = [SC_Gyr_norm_rms; SC_Gyr_norm_rms_]
    DLS_R_Gyr_x_rms = [DLS_R_Gyr_x_rms; DLS_R_Gyr_x_rms_];
    DLS_R_Gyr_y_rms = [DLS_R_Gyr_y_rms; DLS_R_Gyr_y_rms_];
    DLS_R_Gyr_z_rms = [DLS_R_Gyr_z_rms; DLS_R_Gyr_z_rms_];
    DLS_R_Gyr_norm_rms = [DLS_R_Gyr_norm_rms; DLS_R_Gyr_norm_rms_];
    DLS_L_Gyr_x_rms = [DLS_L_Gyr_x_rms; DLS_L_Gyr_x_rms_];
    DLS_L_Gyr_y_rms = [DLS_L_Gyr_y_rms; DLS_L_Gyr_y_rms_];
    DLS_L_Gyr_z_rms = [DLS_L_Gyr_z_rms; DLS_L_Gyr_z_rms_];
    DLS_L_Gyr_norm_rms = [DLS_L_Gyr_norm_rms; DLS_L_Gyr_norm_rms_];
    
    % Standard Deviation
    SC_Gyr_x_std = [SC_Gyr_x_std; SC_Gyr_x_std_]
    SC_Gyr_y_std = [SC_Gyr_y_std; SC_Gyr_y_std_]
    SC_Gyr_z_std = [SC_Gyr_z_std; SC_Gyr_z_std_]
    SC_Gyr_norm_std = [SC_Gyr_norm_std; SC_Gyr_norm_std_]
    DLS_R_Gyr_x_std = [DLS_R_Gyr_x_std; DLS_R_Gyr_x_std_];
    DLS_R_Gyr_y_std = [DLS_R_Gyr_y_std; DLS_R_Gyr_y_std_];
    DLS_R_Gyr_z_std = [DLS_R_Gyr_z_std; DLS_R_Gyr_z_std_];
    DLS_R_Gyr_norm_std = [DLS_R_Gyr_norm_std; DLS_R_Gyr_norm_std_];
    DLS_L_Gyr_x_std = [DLS_L_Gyr_x_std; DLS_L_Gyr_x_std_];
    DLS_L_Gyr_y_std = [DLS_L_Gyr_y_std; DLS_L_Gyr_y_std_];
    DLS_L_Gyr_z_std = [DLS_L_Gyr_z_std; DLS_L_Gyr_z_std_];
    DLS_L_Gyr_norm_std = [DLS_L_Gyr_norm_std; DLS_L_Gyr_norm_std_];

    % Skew
    SC_Gyr_x_skew = [SC_Gyr_x_skew; SC_Gyr_x_skew_]
    SC_Gyr_y_skew = [SC_Gyr_y_skew; SC_Gyr_y_skew_]
    SC_Gyr_z_skew = [SC_Gyr_z_skew; SC_Gyr_z_skew_]
    SC_Gyr_norm_skew = [SC_Gyr_norm_skew; SC_Gyr_norm_skew_]
    DLS_R_Gyr_x_skew = [DLS_R_Gyr_x_skew; DLS_R_Gyr_x_skew_];
    DLS_R_Gyr_y_skew = [DLS_R_Gyr_y_skew; DLS_R_Gyr_y_skew_];
    DLS_R_Gyr_z_skew = [DLS_R_Gyr_z_skew; DLS_R_Gyr_z_skew_];
    DLS_R_Gyr_norm_skew = [DLS_R_Gyr_norm_skew; DLS_R_Gyr_norm_skew_];
    DLS_L_Gyr_x_skew = [DLS_L_Gyr_x_skew; DLS_L_Gyr_x_skew_];
    DLS_L_Gyr_y_skew = [DLS_L_Gyr_y_skew; DLS_L_Gyr_y_skew_];
    DLS_L_Gyr_z_skew = [DLS_L_Gyr_z_skew; DLS_L_Gyr_z_skew_];
    DLS_L_Gyr_norm_skew = [DLS_L_Gyr_norm_skew; DLS_L_Gyr_norm_skew_];

    % Kurtosis
    SC_Gyr_x_kurtosis = [SC_Gyr_x_kurtosis; SC_Gyr_x_kurtosis_]
    SC_Gyr_y_kurtosis = [SC_Gyr_y_kurtosis; SC_Gyr_y_kurtosis_]
    SC_Gyr_z_kurtosis = [SC_Gyr_z_kurtosis; SC_Gyr_z_kurtosis_]
    SC_Gyr_norm_kurtosis = [SC_Gyr_norm_kurtosis; SC_Gyr_norm_kurtosis_]
    DLS_R_Gyr_x_kurtosis = [DLS_R_Gyr_x_kurtosis; DLS_R_Gyr_x_kurtosis_];
    DLS_R_Gyr_y_kurtosis = [DLS_R_Gyr_y_kurtosis; DLS_R_Gyr_y_kurtosis_];
    DLS_R_Gyr_z_kurtosis = [DLS_R_Gyr_z_kurtosis; DLS_R_Gyr_z_kurtosis_];
    DLS_R_Gyr_norm_kurtosis = [DLS_R_Gyr_norm_kurtosis; DLS_R_Gyr_norm_kurtosis_];
    DLS_L_Gyr_x_kurtosis = [DLS_L_Gyr_x_kurtosis; DLS_L_Gyr_x_kurtosis_];
    DLS_L_Gyr_y_kurtosis = [DLS_L_Gyr_y_kurtosis; DLS_L_Gyr_y_kurtosis_];
    DLS_L_Gyr_z_kurtosis = [DLS_L_Gyr_z_kurtosis; DLS_L_Gyr_z_kurtosis_];
    DLS_L_Gyr_norm_kurtosis = [DLS_L_Gyr_norm_kurtosis; DLS_L_Gyr_norm_kurtosis_];

    % Pearson correlation coefficient
    SC_Gyr_corr_xy = [SC_Gyr_corr_xy; SC_Gyr_corr_xy_];
    SC_Gyr_corr_xz = [SC_Gyr_corr_xz; SC_Gyr_corr_xz_];
    SC_Gyr_corr_yz = [SC_Gyr_corr_yz; SC_Gyr_corr_yz_];
    DLS_R_Gyr_corr_xy = [DLS_R_Gyr_corr_xy; DLS_R_Gyr_corr_xy_];
    DLS_R_Gyr_corr_xz = [DLS_R_Gyr_corr_xz; DLS_R_Gyr_corr_xz_];
    DLS_R_Gyr_corr_yz = [DLS_R_Gyr_corr_yz; DLS_R_Gyr_corr_yz_];
    DLS_L_Gyr_corr_xy = [DLS_L_Gyr_corr_xy; DLS_L_Gyr_corr_xy_];
    DLS_L_Gyr_corr_xz = [DLS_L_Gyr_corr_xz; DLS_L_Gyr_corr_xz_];
    DLS_L_Gyr_corr_yz = [DLS_L_Gyr_corr_yz; DLS_L_Gyr_corr_yz_];
    
    % Sample Entropy
    SC_Gyr_x_SamEn = [SC_Gyr_x_SamEn; SC_Gyr_x_SamEn_];
    SC_Gyr_y_SamEn = [SC_Gyr_y_SamEn; SC_Gyr_y_SamEn_]
    SC_Gyr_z_SamEn = [SC_Gyr_z_SamEn; SC_Gyr_z_SamEn_]
    SC_Gyr_norm_SamEn = [SC_Gyr_norm_SamEn; SC_Gyr_norm_SamEn_]
    DLS_R_Gyr_x_SamEn = [DLS_R_Gyr_x_SamEn; DLS_R_Gyr_x_SamEn_];
    DLS_R_Gyr_y_SamEn = [DLS_R_Gyr_y_SamEn; DLS_R_Gyr_y_SamEn_];
    DLS_R_Gyr_z_SamEn = [DLS_R_Gyr_z_SamEn; DLS_R_Gyr_z_SamEn_];
    DLS_R_Gyr_norm_SamEn = [DLS_R_Gyr_norm_SamEn; DLS_R_Gyr_norm_SamEn_];
    DLS_L_Gyr_x_SamEn = [DLS_L_Gyr_x_SamEn; DLS_L_Gyr_x_SamEn_];
    DLS_L_Gyr_y_SamEn = [DLS_L_Gyr_y_SamEn; DLS_L_Gyr_y_SamEn_];
    DLS_L_Gyr_z_SamEn = [DLS_L_Gyr_z_SamEn; DLS_L_Gyr_z_SamEn_];
    DLS_L_Gyr_norm_SamEn = [DLS_L_Gyr_norm_SamEn; DLS_L_Gyr_norm_SamEn_];
    
    % Frequency domain
    SC_Gyr_x_DAmp = [SC_Gyr_x_DAmp; SC_Gyr_x_DAmp_];
    SC_Gyr_x_DFreq = [SC_Gyr_x_DFreq; SC_Gyr_x_DFreq_];
    SC_Gyr_x_PSD_mean = [SC_Gyr_x_PSD_mean; SC_Gyr_x_PSD_mean_];
    SC_Gyr_x_PSD_std = [SC_Gyr_x_PSD_std; SC_Gyr_x_PSD_std_];
    SC_Gyr_x_PSD_skew = [SC_Gyr_x_PSD_skew; SC_Gyr_x_PSD_skew_];
    SC_Gyr_x_PSD_kurtosis = [SC_Gyr_x_PSD_kurtosis; SC_Gyr_x_PSD_kurtosis_];

    SC_Gyr_y_DAmp = [SC_Gyr_y_DAmp; SC_Gyr_y_DAmp_];
    SC_Gyr_y_DFreq = [SC_Gyr_y_DFreq; SC_Gyr_y_DFreq_];
    SC_Gyr_y_PSD_mean = [SC_Gyr_y_PSD_mean; SC_Gyr_y_PSD_mean_];
    SC_Gyr_y_PSD_std = [SC_Gyr_y_PSD_std; SC_Gyr_y_PSD_std_];
    SC_Gyr_y_PSD_skew = [SC_Gyr_y_PSD_skew; SC_Gyr_y_PSD_skew_];
    SC_Gyr_y_PSD_kurtosis = [SC_Gyr_y_PSD_kurtosis; SC_Gyr_y_PSD_kurtosis_];

    SC_Gyr_z_DAmp = [SC_Gyr_z_DAmp; SC_Gyr_z_DAmp_];
    SC_Gyr_z_DFreq = [SC_Gyr_z_DFreq; SC_Gyr_z_DFreq_];
    SC_Gyr_z_PSD_mean = [SC_Gyr_z_PSD_mean; SC_Gyr_z_PSD_mean_];
    SC_Gyr_z_PSD_std = [SC_Gyr_z_PSD_std; SC_Gyr_z_PSD_std_];
    SC_Gyr_z_PSD_skew = [SC_Gyr_z_PSD_skew; SC_Gyr_z_PSD_skew_];
    SC_Gyr_z_PSD_kurtosis = [SC_Gyr_z_PSD_kurtosis; SC_Gyr_z_PSD_kurtosis_];

    SC_Gyr_norm_DAmp = [SC_Gyr_norm_DAmp; SC_Gyr_norm_DAmp_];
    SC_Gyr_norm_DFreq = [SC_Gyr_norm_DFreq; SC_Gyr_norm_DFreq_];
    SC_Gyr_norm_PSD_mean = [SC_Gyr_norm_PSD_mean; SC_Gyr_norm_PSD_mean_];
    SC_Gyr_norm_PSD_std = [SC_Gyr_norm_PSD_std; SC_Gyr_norm_PSD_std_];
    SC_Gyr_norm_PSD_skew = [SC_Gyr_norm_PSD_skew; SC_Gyr_norm_PSD_skew_];
    SC_Gyr_norm_PSD_kurtosis = [SC_Gyr_norm_PSD_kurtosis; SC_Gyr_norm_PSD_kurtosis_];
    
    DLS_R_Gyr_x_DAmp = [DLS_R_Gyr_x_DAmp; DLS_R_Gyr_x_DAmp_];
    DLS_R_Gyr_x_DFreq = [DLS_R_Gyr_x_DFreq; DLS_R_Gyr_x_DFreq_];
    DLS_R_Gyr_x_PSD_mean = [DLS_R_Gyr_x_PSD_mean; DLS_R_Gyr_x_PSD_mean_];
    DLS_R_Gyr_x_PSD_std = [DLS_R_Gyr_x_PSD_std; DLS_R_Gyr_x_PSD_std_];
    DLS_R_Gyr_x_PSD_skew = [DLS_R_Gyr_x_PSD_skew; DLS_R_Gyr_x_PSD_skew_];
    DLS_R_Gyr_x_PSD_kurtosis = [DLS_R_Gyr_x_PSD_kurtosis; DLS_R_Gyr_x_PSD_kurtosis_];

    DLS_R_Gyr_y_DAmp = [DLS_R_Gyr_y_DAmp; DLS_R_Gyr_y_DAmp_];
    DLS_R_Gyr_y_DFreq = [DLS_R_Gyr_y_DFreq; DLS_R_Gyr_y_DFreq_];
    DLS_R_Gyr_y_PSD_mean = [DLS_R_Gyr_y_PSD_mean; DLS_R_Gyr_y_PSD_mean_];
    DLS_R_Gyr_y_PSD_std = [DLS_R_Gyr_y_PSD_std; DLS_R_Gyr_y_PSD_std_];
    DLS_R_Gyr_y_PSD_skew = [DLS_R_Gyr_y_PSD_skew; DLS_R_Gyr_y_PSD_skew_];
    DLS_R_Gyr_y_PSD_kurtosis = [DLS_R_Gyr_y_PSD_kurtosis; DLS_R_Gyr_y_PSD_kurtosis_];

    DLS_R_Gyr_z_DAmp = [DLS_R_Gyr_z_DAmp; DLS_R_Gyr_z_DAmp_];
    DLS_R_Gyr_z_DFreq = [DLS_R_Gyr_z_DFreq; DLS_R_Gyr_z_DFreq_];
    DLS_R_Gyr_z_PSD_mean = [DLS_R_Gyr_z_PSD_mean; DLS_R_Gyr_z_PSD_mean_];
    DLS_R_Gyr_z_PSD_std = [DLS_R_Gyr_z_PSD_std; DLS_R_Gyr_z_PSD_std_];
    DLS_R_Gyr_z_PSD_skew = [DLS_R_Gyr_z_PSD_skew; DLS_R_Gyr_z_PSD_skew_];
    DLS_R_Gyr_z_PSD_kurtosis = [DLS_R_Gyr_z_PSD_kurtosis; DLS_R_Gyr_z_PSD_kurtosis_];

    DLS_R_Gyr_norm_DAmp = [DLS_R_Gyr_norm_DAmp; DLS_R_Gyr_norm_DAmp_];
    DLS_R_Gyr_norm_DFreq = [DLS_R_Gyr_norm_DFreq; DLS_R_Gyr_norm_DFreq_];
    DLS_R_Gyr_norm_PSD_mean = [DLS_R_Gyr_norm_PSD_mean; DLS_R_Gyr_norm_PSD_mean_];
    DLS_R_Gyr_norm_PSD_std = [DLS_R_Gyr_norm_PSD_std; DLS_R_Gyr_norm_PSD_std_];
    DLS_R_Gyr_norm_PSD_skew = [DLS_R_Gyr_norm_PSD_skew; DLS_R_Gyr_norm_PSD_skew_];
    DLS_R_Gyr_norm_PSD_kurtosis = [DLS_R_Gyr_norm_PSD_kurtosis; DLS_R_Gyr_norm_PSD_kurtosis_];

    DLS_L_Gyr_x_DAmp = [DLS_L_Gyr_x_DAmp; DLS_L_Gyr_x_DAmp_];
    DLS_L_Gyr_x_DFreq = [DLS_L_Gyr_x_DFreq; DLS_L_Gyr_x_DFreq_];
    DLS_L_Gyr_x_PSD_mean = [DLS_L_Gyr_x_PSD_mean; DLS_L_Gyr_x_PSD_mean_];
    DLS_L_Gyr_x_PSD_std = [DLS_L_Gyr_x_PSD_std; DLS_L_Gyr_x_PSD_std_];
    DLS_L_Gyr_x_PSD_skew = [DLS_L_Gyr_x_PSD_skew; DLS_L_Gyr_x_PSD_skew_];
    DLS_L_Gyr_x_PSD_kurtosis = [DLS_L_Gyr_x_PSD_kurtosis; DLS_L_Gyr_x_PSD_kurtosis_];

    DLS_L_Gyr_y_DAmp = [DLS_L_Gyr_y_DAmp; DLS_L_Gyr_y_DAmp_];
    DLS_L_Gyr_y_DFreq = [DLS_L_Gyr_y_DFreq; DLS_L_Gyr_y_DFreq_];
    DLS_L_Gyr_y_PSD_mean = [DLS_L_Gyr_y_PSD_mean; DLS_L_Gyr_y_PSD_mean_];
    DLS_L_Gyr_y_PSD_std = [DLS_L_Gyr_y_PSD_std; DLS_L_Gyr_y_PSD_std_];
    DLS_L_Gyr_y_PSD_skew = [DLS_L_Gyr_y_PSD_skew; DLS_L_Gyr_y_PSD_skew_];
    DLS_L_Gyr_y_PSD_kurtosis = [DLS_L_Gyr_y_PSD_kurtosis; DLS_L_Gyr_y_PSD_kurtosis_];

    DLS_L_Gyr_z_DAmp = [DLS_L_Gyr_z_DAmp; DLS_L_Gyr_z_DAmp_];
    DLS_L_Gyr_z_DFreq = [DLS_L_Gyr_z_DFreq; DLS_L_Gyr_z_DFreq_];
    DLS_L_Gyr_z_PSD_mean = [DLS_L_Gyr_z_PSD_mean; DLS_L_Gyr_z_PSD_mean_];
    DLS_L_Gyr_z_PSD_std = [DLS_L_Gyr_z_PSD_std; DLS_L_Gyr_z_PSD_std_];
    DLS_L_Gyr_z_PSD_skew = [DLS_L_Gyr_z_PSD_skew; DLS_L_Gyr_z_PSD_skew_];
    DLS_L_Gyr_z_PSD_kurtosis = [DLS_L_Gyr_z_PSD_kurtosis; DLS_L_Gyr_z_PSD_kurtosis_];

    DLS_L_Gyr_norm_DAmp = [DLS_L_Gyr_norm_DAmp; DLS_L_Gyr_norm_DAmp_];
    DLS_L_Gyr_norm_DFreq = [DLS_L_Gyr_norm_DFreq; DLS_L_Gyr_norm_DFreq_];
    DLS_L_Gyr_norm_PSD_mean = [DLS_L_Gyr_norm_PSD_mean; DLS_L_Gyr_norm_PSD_mean_];
    DLS_L_Gyr_norm_PSD_std = [DLS_L_Gyr_norm_PSD_std; DLS_L_Gyr_norm_PSD_std_];
    DLS_L_Gyr_norm_PSD_skew = [DLS_L_Gyr_norm_PSD_skew; DLS_L_Gyr_norm_PSD_skew_];
    DLS_L_Gyr_norm_PSD_kurtosis = [DLS_L_Gyr_norm_PSD_kurtosis; DLS_L_Gyr_norm_PSD_kurtosis_];


    % Acc mean
    % General Features
    SC_Acc_x_mean = [SC_Acc_x_mean; SC_Acc_x_mean_]
    SC_Acc_y_mean = [SC_Acc_y_mean; SC_Acc_y_mean_]
    SC_Acc_z_mean = [SC_Acc_z_mean; SC_Acc_z_mean_]
    SC_Acc_norm_mean = [SC_Acc_norm_mean; SC_Acc_norm_mean_]
    DLS_R_Acc_x_mean = [DLS_R_Acc_x_mean; DLS_R_Acc_x_mean_];
    DLS_R_Acc_y_mean = [DLS_R_Acc_y_mean; DLS_R_Acc_y_mean_];
    DLS_R_Acc_z_mean = [DLS_R_Acc_z_mean; DLS_R_Acc_z_mean_];
    DLS_R_Acc_norm_mean = [DLS_R_Acc_norm_mean; DLS_R_Acc_norm_mean_];
    DLS_L_Acc_x_mean = [DLS_L_Acc_x_mean; DLS_L_Acc_x_mean_];
    DLS_L_Acc_y_mean = [DLS_L_Acc_y_mean; DLS_L_Acc_y_mean_];
    DLS_L_Acc_z_mean = [DLS_L_Acc_z_mean; DLS_L_Acc_z_mean_];
    DLS_L_Acc_norm_mean = [DLS_L_Acc_norm_mean; DLS_L_Acc_norm_mean_];
    
    % Range
    SC_Acc_x_range = [SC_Acc_x_range; SC_Acc_x_range_]
    SC_Acc_y_range = [SC_Acc_y_range; SC_Acc_y_range_]
    SC_Acc_z_range = [SC_Acc_z_range; SC_Acc_z_range_]
    SC_Acc_norm_range = [SC_Acc_norm_range; SC_Acc_norm_range_]
    DLS_R_Acc_x_range = [DLS_R_Acc_x_range; DLS_R_Acc_x_range_];
    DLS_R_Acc_y_range = [DLS_R_Acc_y_range; DLS_R_Acc_y_range_];
    DLS_R_Acc_z_range = [DLS_R_Acc_z_range; DLS_R_Acc_z_range_];
    DLS_R_Acc_norm_range = [DLS_R_Acc_norm_range; DLS_R_Acc_norm_range_];
    DLS_L_Acc_x_range = [DLS_L_Acc_x_range; DLS_L_Acc_x_range_];
    DLS_L_Acc_y_range = [DLS_L_Acc_y_range; DLS_L_Acc_y_range_];
    DLS_L_Acc_z_range = [DLS_L_Acc_z_range; DLS_L_Acc_z_range_];
    DLS_L_Acc_norm_range = [DLS_L_Acc_norm_range; DLS_L_Acc_norm_range_];

    % RMS
    SC_Acc_x_rms = [SC_Acc_x_rms; SC_Acc_x_rms_]
    SC_Acc_y_rms = [SC_Acc_y_rms; SC_Acc_y_rms_]
    SC_Acc_z_rms = [SC_Acc_z_rms; SC_Acc_z_rms_]
    SC_Acc_norm_rms = [SC_Acc_norm_rms; SC_Acc_norm_rms_]
    DLS_R_Acc_x_rms = [DLS_R_Acc_x_rms; DLS_R_Acc_x_rms_];
    DLS_R_Acc_y_rms = [DLS_R_Acc_y_rms; DLS_R_Acc_y_rms_];
    DLS_R_Acc_z_rms = [DLS_R_Acc_z_rms; DLS_R_Acc_z_rms_];
    DLS_R_Acc_norm_rms = [DLS_R_Acc_norm_rms; DLS_R_Acc_norm_rms_];
    DLS_L_Acc_x_rms = [DLS_L_Acc_x_rms; DLS_L_Acc_x_rms_];
    DLS_L_Acc_y_rms = [DLS_L_Acc_y_rms; DLS_L_Acc_y_rms_];
    DLS_L_Acc_z_rms = [DLS_L_Acc_z_rms; DLS_L_Acc_z_rms_];
    DLS_L_Acc_norm_rms = [DLS_L_Acc_norm_rms; DLS_L_Acc_norm_rms_];
    
    % Standard Deviation
    SC_Acc_x_std = [SC_Acc_x_std; SC_Acc_x_std_]
    SC_Acc_y_std = [SC_Acc_y_std; SC_Acc_y_std_]
    SC_Acc_z_std = [SC_Acc_z_std; SC_Acc_z_std_]
    SC_Acc_norm_std = [SC_Acc_norm_std; SC_Acc_norm_std_]
    DLS_R_Acc_x_std = [DLS_R_Acc_x_std; DLS_R_Acc_x_std_];
    DLS_R_Acc_y_std = [DLS_R_Acc_y_std; DLS_R_Acc_y_std_];
    DLS_R_Acc_z_std = [DLS_R_Acc_z_std; DLS_R_Acc_z_std_];
    DLS_R_Acc_norm_std = [DLS_R_Acc_norm_std; DLS_R_Acc_norm_std_];
    DLS_L_Acc_x_std = [DLS_L_Acc_x_std; DLS_L_Acc_x_std_];
    DLS_L_Acc_y_std = [DLS_L_Acc_y_std; DLS_L_Acc_y_std_];
    DLS_L_Acc_z_std = [DLS_L_Acc_z_std; DLS_L_Acc_z_std_];
    DLS_L_Acc_norm_std = [DLS_L_Acc_norm_std; DLS_L_Acc_norm_std_];

    % Skew
    SC_Acc_x_skew = [SC_Acc_x_skew; SC_Acc_x_skew_]
    SC_Acc_y_skew = [SC_Acc_y_skew; SC_Acc_y_skew_]
    SC_Acc_z_skew = [SC_Acc_z_skew; SC_Acc_z_skew_]
    SC_Acc_norm_skew = [SC_Acc_norm_skew; SC_Acc_norm_skew_]
    DLS_R_Acc_x_skew = [DLS_R_Acc_x_skew; DLS_R_Acc_x_skew_];
    DLS_R_Acc_y_skew = [DLS_R_Acc_y_skew; DLS_R_Acc_y_skew_];
    DLS_R_Acc_z_skew = [DLS_R_Acc_z_skew; DLS_R_Acc_z_skew_];
    DLS_R_Acc_norm_skew = [DLS_R_Acc_norm_skew; DLS_R_Acc_norm_skew_];
    DLS_L_Acc_x_skew = [DLS_L_Acc_x_skew; DLS_L_Acc_x_skew_];
    DLS_L_Acc_y_skew = [DLS_L_Acc_y_skew; DLS_L_Acc_y_skew_];
    DLS_L_Acc_z_skew = [DLS_L_Acc_z_skew; DLS_L_Acc_z_skew_];
    DLS_L_Acc_norm_skew = [DLS_L_Acc_norm_skew; DLS_L_Acc_norm_skew_];

    % Kurtosis
    SC_Acc_x_kurtosis = [SC_Acc_x_kurtosis; SC_Acc_x_kurtosis_]
    SC_Acc_y_kurtosis = [SC_Acc_y_kurtosis; SC_Acc_y_kurtosis_]
    SC_Acc_z_kurtosis = [SC_Acc_z_kurtosis; SC_Acc_z_kurtosis_]
    SC_Acc_norm_kurtosis = [SC_Acc_norm_kurtosis; SC_Acc_norm_kurtosis_]
    DLS_R_Acc_x_kurtosis = [DLS_R_Acc_x_kurtosis; DLS_R_Acc_x_kurtosis_];
    DLS_R_Acc_y_kurtosis = [DLS_R_Acc_y_kurtosis; DLS_R_Acc_y_kurtosis_];
    DLS_R_Acc_z_kurtosis = [DLS_R_Acc_z_kurtosis; DLS_R_Acc_z_kurtosis_];
    DLS_R_Acc_norm_kurtosis = [DLS_R_Acc_norm_kurtosis; DLS_R_Acc_norm_kurtosis_];
    DLS_L_Acc_x_kurtosis = [DLS_L_Acc_x_kurtosis; DLS_L_Acc_x_kurtosis_];
    DLS_L_Acc_y_kurtosis = [DLS_L_Acc_y_kurtosis; DLS_L_Acc_y_kurtosis_];
    DLS_L_Acc_z_kurtosis = [DLS_L_Acc_z_kurtosis; DLS_L_Acc_z_kurtosis_];
    DLS_L_Acc_norm_kurtosis = [DLS_L_Acc_norm_kurtosis; DLS_L_Acc_norm_kurtosis_];

    % Pearson correlation coefficient
    SC_Acc_corr_xy = [SC_Acc_corr_xy; SC_Acc_corr_xy_];
    SC_Acc_corr_xz = [SC_Acc_corr_xz; SC_Acc_corr_xz_];
    SC_Acc_corr_yz = [SC_Acc_corr_yz; SC_Acc_corr_yz_];
    DLS_R_Acc_corr_xy = [DLS_R_Acc_corr_xy; DLS_R_Acc_corr_xy_];
    DLS_R_Acc_corr_xz = [DLS_R_Acc_corr_xz; DLS_R_Acc_corr_xz_];
    DLS_R_Acc_corr_yz = [DLS_R_Acc_corr_yz; DLS_R_Acc_corr_yz_];
    DLS_L_Acc_corr_xy = [DLS_L_Acc_corr_xy; DLS_L_Acc_corr_xy_];
    DLS_L_Acc_corr_xz = [DLS_L_Acc_corr_xz; DLS_L_Acc_corr_xz_];
    DLS_L_Acc_corr_yz = [DLS_L_Acc_corr_yz; DLS_L_Acc_corr_yz_];
    
    % Sample Entropy
    SC_Acc_x_SamEn = [SC_Acc_x_SamEn; SC_Acc_x_SamEn_];
    SC_Acc_y_SamEn = [SC_Acc_y_SamEn; SC_Acc_y_SamEn_]
    SC_Acc_z_SamEn = [SC_Acc_z_SamEn; SC_Acc_z_SamEn_]
    SC_Acc_norm_SamEn = [SC_Acc_norm_SamEn; SC_Acc_norm_SamEn_]
    DLS_R_Acc_x_SamEn = [DLS_R_Acc_x_SamEn; DLS_R_Acc_x_SamEn_];
    DLS_R_Acc_y_SamEn = [DLS_R_Acc_y_SamEn; DLS_R_Acc_y_SamEn_];
    DLS_R_Acc_z_SamEn = [DLS_R_Acc_z_SamEn; DLS_R_Acc_z_SamEn_];
    DLS_R_Acc_norm_SamEn = [DLS_R_Acc_norm_SamEn; DLS_R_Acc_norm_SamEn_];
    DLS_L_Acc_x_SamEn = [DLS_L_Acc_x_SamEn; DLS_L_Acc_x_SamEn_];
    DLS_L_Acc_y_SamEn = [DLS_L_Acc_y_SamEn; DLS_L_Acc_y_SamEn_];
    DLS_L_Acc_z_SamEn = [DLS_L_Acc_z_SamEn; DLS_L_Acc_z_SamEn_];
    DLS_L_Acc_norm_SamEn = [DLS_L_Acc_norm_SamEn; DLS_L_Acc_norm_SamEn_];
    
    % Frequency domain
    SC_Acc_x_DAmp = [SC_Acc_x_DAmp; SC_Acc_x_DAmp_];
    SC_Acc_x_DFreq = [SC_Acc_x_DFreq; SC_Acc_x_DFreq_];
    SC_Acc_x_PSD_mean = [SC_Acc_x_PSD_mean; SC_Acc_x_PSD_mean_];
    SC_Acc_x_PSD_std = [SC_Acc_x_PSD_std; SC_Acc_x_PSD_std_];
    SC_Acc_x_PSD_skew = [SC_Acc_x_PSD_skew; SC_Acc_x_PSD_skew_];
    SC_Acc_x_PSD_kurtosis = [SC_Acc_x_PSD_kurtosis; SC_Acc_x_PSD_kurtosis_];

    SC_Acc_y_DAmp = [SC_Acc_y_DAmp; SC_Acc_y_DAmp_];
    SC_Acc_y_DFreq = [SC_Acc_y_DFreq; SC_Acc_y_DFreq_];
    SC_Acc_y_PSD_mean = [SC_Acc_y_PSD_mean; SC_Acc_y_PSD_mean_];
    SC_Acc_y_PSD_std = [SC_Acc_y_PSD_std; SC_Acc_y_PSD_std_];
    SC_Acc_y_PSD_skew = [SC_Acc_y_PSD_skew; SC_Acc_y_PSD_skew_];
    SC_Acc_y_PSD_kurtosis = [SC_Acc_y_PSD_kurtosis; SC_Acc_y_PSD_kurtosis_];

    SC_Acc_z_DAmp = [SC_Acc_z_DAmp; SC_Acc_z_DAmp_];
    SC_Acc_z_DFreq = [SC_Acc_z_DFreq; SC_Acc_z_DFreq_];
    SC_Acc_z_PSD_mean = [SC_Acc_z_PSD_mean; SC_Acc_z_PSD_mean_];
    SC_Acc_z_PSD_std = [SC_Acc_z_PSD_std; SC_Acc_z_PSD_std_];
    SC_Acc_z_PSD_skew = [SC_Acc_z_PSD_skew; SC_Acc_z_PSD_skew_];
    SC_Acc_z_PSD_kurtosis = [SC_Acc_z_PSD_kurtosis; SC_Acc_z_PSD_kurtosis_];

    SC_Acc_norm_DAmp = [SC_Acc_norm_DAmp; SC_Acc_norm_DAmp_];
    SC_Acc_norm_DFreq = [SC_Acc_norm_DFreq; SC_Acc_norm_DFreq_];
    SC_Acc_norm_PSD_mean = [SC_Acc_norm_PSD_mean; SC_Acc_norm_PSD_mean_];
    SC_Acc_norm_PSD_std = [SC_Acc_norm_PSD_std; SC_Acc_norm_PSD_std_];
    SC_Acc_norm_PSD_skew = [SC_Acc_norm_PSD_skew; SC_Acc_norm_PSD_skew_];
    SC_Acc_norm_PSD_kurtosis = [SC_Acc_norm_PSD_kurtosis; SC_Acc_norm_PSD_kurtosis_];
    
    DLS_R_Acc_x_DAmp = [DLS_R_Acc_x_DAmp; DLS_R_Acc_x_DAmp_];
    DLS_R_Acc_x_DFreq = [DLS_R_Acc_x_DFreq; DLS_R_Acc_x_DFreq_];
    DLS_R_Acc_x_PSD_mean = [DLS_R_Acc_x_PSD_mean; DLS_R_Acc_x_PSD_mean_];
    DLS_R_Acc_x_PSD_std = [DLS_R_Acc_x_PSD_std; DLS_R_Acc_x_PSD_std_];
    DLS_R_Acc_x_PSD_skew = [DLS_R_Acc_x_PSD_skew; DLS_R_Acc_x_PSD_skew_];
    DLS_R_Acc_x_PSD_kurtosis = [DLS_R_Acc_x_PSD_kurtosis; DLS_R_Acc_x_PSD_kurtosis_];

    DLS_R_Acc_y_DAmp = [DLS_R_Acc_y_DAmp; DLS_R_Acc_y_DAmp_];
    DLS_R_Acc_y_DFreq = [DLS_R_Acc_y_DFreq; DLS_R_Acc_y_DFreq_];
    DLS_R_Acc_y_PSD_mean = [DLS_R_Acc_y_PSD_mean; DLS_R_Acc_y_PSD_mean_];
    DLS_R_Acc_y_PSD_std = [DLS_R_Acc_y_PSD_std; DLS_R_Acc_y_PSD_std_];
    DLS_R_Acc_y_PSD_skew = [DLS_R_Acc_y_PSD_skew; DLS_R_Acc_y_PSD_skew_];
    DLS_R_Acc_y_PSD_kurtosis = [DLS_R_Acc_y_PSD_kurtosis; DLS_R_Acc_y_PSD_kurtosis_];

    DLS_R_Acc_z_DAmp = [DLS_R_Acc_z_DAmp; DLS_R_Acc_z_DAmp_];
    DLS_R_Acc_z_DFreq = [DLS_R_Acc_z_DFreq; DLS_R_Acc_z_DFreq_];
    DLS_R_Acc_z_PSD_mean = [DLS_R_Acc_z_PSD_mean; DLS_R_Acc_z_PSD_mean_];
    DLS_R_Acc_z_PSD_std = [DLS_R_Acc_z_PSD_std; DLS_R_Acc_z_PSD_std_];
    DLS_R_Acc_z_PSD_skew = [DLS_R_Acc_z_PSD_skew; DLS_R_Acc_z_PSD_skew_];
    DLS_R_Acc_z_PSD_kurtosis = [DLS_R_Acc_z_PSD_kurtosis; DLS_R_Acc_z_PSD_kurtosis_];

    DLS_R_Acc_norm_DAmp = [DLS_R_Acc_norm_DAmp; DLS_R_Acc_norm_DAmp_];
    DLS_R_Acc_norm_DFreq = [DLS_R_Acc_norm_DFreq; DLS_R_Acc_norm_DFreq_];
    DLS_R_Acc_norm_PSD_mean = [DLS_R_Acc_norm_PSD_mean; DLS_R_Acc_norm_PSD_mean_];
    DLS_R_Acc_norm_PSD_std = [DLS_R_Acc_norm_PSD_std; DLS_R_Acc_norm_PSD_std_];
    DLS_R_Acc_norm_PSD_skew = [DLS_R_Acc_norm_PSD_skew; DLS_R_Acc_norm_PSD_skew_];
    DLS_R_Acc_norm_PSD_kurtosis = [DLS_R_Acc_norm_PSD_kurtosis; DLS_R_Acc_norm_PSD_kurtosis_];

    DLS_L_Acc_x_DAmp = [DLS_L_Acc_x_DAmp; DLS_L_Acc_x_DAmp_];
    DLS_L_Acc_x_DFreq = [DLS_L_Acc_x_DFreq; DLS_L_Acc_x_DFreq_];
    DLS_L_Acc_x_PSD_mean = [DLS_L_Acc_x_PSD_mean; DLS_L_Acc_x_PSD_mean_];
    DLS_L_Acc_x_PSD_std = [DLS_L_Acc_x_PSD_std; DLS_L_Acc_x_PSD_std_];
    DLS_L_Acc_x_PSD_skew = [DLS_L_Acc_x_PSD_skew; DLS_L_Acc_x_PSD_skew_];
    DLS_L_Acc_x_PSD_kurtosis = [DLS_L_Acc_x_PSD_kurtosis; DLS_L_Acc_x_PSD_kurtosis_];

    DLS_L_Acc_y_DAmp = [DLS_L_Acc_y_DAmp; DLS_L_Acc_y_DAmp_];
    DLS_L_Acc_y_DFreq = [DLS_L_Acc_y_DFreq; DLS_L_Acc_y_DFreq_];
    DLS_L_Acc_y_PSD_mean = [DLS_L_Acc_y_PSD_mean; DLS_L_Acc_y_PSD_mean_];
    DLS_L_Acc_y_PSD_std = [DLS_L_Acc_y_PSD_std; DLS_L_Acc_y_PSD_std_];
    DLS_L_Acc_y_PSD_skew = [DLS_L_Acc_y_PSD_skew; DLS_L_Acc_y_PSD_skew_];
    DLS_L_Acc_y_PSD_kurtosis = [DLS_L_Acc_y_PSD_kurtosis; DLS_L_Acc_y_PSD_kurtosis_];

    DLS_L_Acc_z_DAmp = [DLS_L_Acc_z_DAmp; DLS_L_Acc_z_DAmp_];
    DLS_L_Acc_z_DFreq = [DLS_L_Acc_z_DFreq; DLS_L_Acc_z_DFreq_];
    DLS_L_Acc_z_PSD_mean = [DLS_L_Acc_z_PSD_mean; DLS_L_Acc_z_PSD_mean_];
    DLS_L_Acc_z_PSD_std = [DLS_L_Acc_z_PSD_std; DLS_L_Acc_z_PSD_std_];
    DLS_L_Acc_z_PSD_skew = [DLS_L_Acc_z_PSD_skew; DLS_L_Acc_z_PSD_skew_];
    DLS_L_Acc_z_PSD_kurtosis = [DLS_L_Acc_z_PSD_kurtosis; DLS_L_Acc_z_PSD_kurtosis_];

    DLS_L_Acc_norm_DAmp = [DLS_L_Acc_norm_DAmp; DLS_L_Acc_norm_DAmp_];
    DLS_L_Acc_norm_DFreq = [DLS_L_Acc_norm_DFreq; DLS_L_Acc_norm_DFreq_];
    DLS_L_Acc_norm_PSD_mean = [DLS_L_Acc_norm_PSD_mean; DLS_L_Acc_norm_PSD_mean_];
    DLS_L_Acc_norm_PSD_std = [DLS_L_Acc_norm_PSD_std; DLS_L_Acc_norm_PSD_std_];
    DLS_L_Acc_norm_PSD_skew = [DLS_L_Acc_norm_PSD_skew; DLS_L_Acc_norm_PSD_skew_];
    DLS_L_Acc_norm_PSD_kurtosis = [DLS_L_Acc_norm_PSD_kurtosis; DLS_L_Acc_norm_PSD_kurtosis_];

  

    
    
end
            
tbl_Features = table(Sub_ID, Sub_Type, Cut_Off_Time, FIM, FIM_M, BBS, MWT10_SSV, MWT10_FV, MWT6, TUG, ...
    AoM_Pel_tilt, AoM_Pel_ro, AoM_Pel_oblq, AoM_Pel_norm, ...
    AoM_Ankle_US_x, AoM_Ankle_US_y, AoM_Ankle_US_z, AoM_Ankle_US_norm, ...
    AoM_Ankle_AS_x, AoM_Ankle_AS_y, AoM_Ankle_AS_z, AoM_Ankle_AS_norm, ...
    SC_Gyr_x_mean, SC_Gyr_y_mean, SC_Gyr_z_mean, SC_Gyr_norm_mean, ...
    DLS_R_Gyr_x_mean, DLS_R_Gyr_y_mean, DLS_R_Gyr_z_mean, DLS_R_Gyr_norm_mean, ...
    DLS_L_Gyr_x_mean, DLS_L_Gyr_y_mean, DLS_L_Gyr_z_mean, DLS_L_Gyr_norm_mean, ...
    SC_Gyr_x_range, SC_Gyr_y_range, SC_Gyr_z_range, SC_Gyr_norm_range, ...
    DLS_R_Gyr_x_range, DLS_R_Gyr_y_range, DLS_R_Gyr_z_range, DLS_R_Gyr_norm_range, ...
    DLS_L_Gyr_x_range, DLS_L_Gyr_y_range, DLS_L_Gyr_z_range, DLS_L_Gyr_norm_range, ...
    SC_Gyr_x_rms, SC_Gyr_y_rms, SC_Gyr_z_rms, SC_Gyr_norm_rms, ...
    DLS_R_Gyr_x_rms, DLS_R_Gyr_y_rms, DLS_R_Gyr_z_rms, DLS_R_Gyr_norm_rms, ...
    DLS_L_Gyr_x_rms, DLS_L_Gyr_y_rms, DLS_L_Gyr_z_rms, DLS_L_Gyr_norm_rms, ...
    SC_Gyr_x_std, SC_Gyr_y_std, SC_Gyr_z_std, SC_Gyr_norm_std, ...
    DLS_R_Gyr_x_std, DLS_R_Gyr_y_std, DLS_R_Gyr_z_std, DLS_R_Gyr_norm_std, ...
    DLS_L_Gyr_x_std, DLS_L_Gyr_y_std, DLS_L_Gyr_z_std, DLS_L_Gyr_norm_std, ...
    SC_Gyr_x_skew, SC_Gyr_y_skew, SC_Gyr_z_skew, SC_Gyr_norm_skew, ...
    DLS_R_Gyr_x_skew, DLS_R_Gyr_y_skew, DLS_R_Gyr_z_skew, DLS_R_Gyr_norm_skew, ...
    DLS_L_Gyr_x_skew, DLS_L_Gyr_y_skew, DLS_L_Gyr_z_skew, DLS_L_Gyr_norm_skew, ...
    SC_Gyr_x_kurtosis, SC_Gyr_y_kurtosis, SC_Gyr_z_kurtosis, SC_Gyr_norm_kurtosis, ...
    DLS_R_Gyr_x_kurtosis, DLS_R_Gyr_y_kurtosis, DLS_R_Gyr_z_kurtosis, DLS_R_Gyr_norm_kurtosis, ...
    DLS_L_Gyr_x_kurtosis, DLS_L_Gyr_y_kurtosis, DLS_L_Gyr_z_kurtosis, DLS_L_Gyr_norm_kurtosis, ...
    SC_Gyr_corr_xy, SC_Gyr_corr_xz, SC_Gyr_corr_yz, ...
    DLS_R_Gyr_corr_xy, DLS_R_Gyr_corr_xz, DLS_R_Gyr_corr_yz, ...
    DLS_L_Gyr_corr_xy, DLS_L_Gyr_corr_xz, DLS_L_Gyr_corr_yz, ...
    SC_Gyr_x_SamEn, SC_Gyr_y_SamEn, SC_Gyr_z_SamEn, SC_Gyr_norm_SamEn, ...
    DLS_R_Gyr_x_SamEn, DLS_R_Gyr_y_SamEn, DLS_R_Gyr_z_SamEn, DLS_R_Gyr_norm_SamEn, ... 
    DLS_L_Gyr_x_SamEn, DLS_L_Gyr_y_SamEn, DLS_L_Gyr_z_SamEn, DLS_L_Gyr_norm_SamEn, ...
    SC_Gyr_x_DAmp, SC_Gyr_x_DFreq, SC_Gyr_x_PSD_mean, SC_Gyr_x_PSD_std, SC_Gyr_x_PSD_skew, SC_Gyr_x_PSD_kurtosis, ...
    SC_Gyr_y_DAmp, SC_Gyr_y_DFreq, SC_Gyr_y_PSD_mean, SC_Gyr_y_PSD_std, SC_Gyr_y_PSD_skew, SC_Gyr_y_PSD_kurtosis, ...
    SC_Gyr_z_DAmp, SC_Gyr_z_DFreq, SC_Gyr_z_PSD_mean, SC_Gyr_z_PSD_std, SC_Gyr_z_PSD_skew, SC_Gyr_z_PSD_kurtosis, ...
    SC_Gyr_norm_DAmp, SC_Gyr_norm_DFreq, SC_Gyr_norm_PSD_mean, SC_Gyr_norm_PSD_std, SC_Gyr_norm_PSD_skew, SC_Gyr_norm_PSD_kurtosis, ...        
    DLS_R_Gyr_x_DAmp, DLS_R_Gyr_x_DFreq, DLS_R_Gyr_x_PSD_mean, DLS_R_Gyr_x_PSD_std, DLS_R_Gyr_x_PSD_skew, DLS_R_Gyr_x_PSD_kurtosis, ...
    DLS_R_Gyr_y_DAmp, DLS_R_Gyr_y_DFreq, DLS_R_Gyr_y_PSD_mean, DLS_R_Gyr_y_PSD_std, DLS_R_Gyr_y_PSD_skew, DLS_R_Gyr_y_PSD_kurtosis, ...
    DLS_R_Gyr_z_DAmp, DLS_R_Gyr_z_DFreq, DLS_R_Gyr_z_PSD_mean, DLS_R_Gyr_z_PSD_std, DLS_R_Gyr_z_PSD_skew, DLS_R_Gyr_z_PSD_kurtosis, ...
    DLS_R_Gyr_norm_DAmp, DLS_R_Gyr_norm_DFreq, DLS_R_Gyr_norm_PSD_mean, DLS_R_Gyr_norm_PSD_std, DLS_R_Gyr_norm_PSD_skew, DLS_R_Gyr_norm_PSD_kurtosis, ...
    DLS_L_Gyr_x_DAmp, DLS_L_Gyr_x_DFreq, DLS_L_Gyr_x_PSD_mean, DLS_L_Gyr_x_PSD_std, DLS_L_Gyr_x_PSD_skew, DLS_L_Gyr_x_PSD_kurtosis, ...
    DLS_L_Gyr_y_DAmp, DLS_L_Gyr_y_DFreq, DLS_L_Gyr_y_PSD_mean, DLS_L_Gyr_y_PSD_std, DLS_L_Gyr_y_PSD_skew, DLS_L_Gyr_y_PSD_kurtosis, ...
    DLS_L_Gyr_z_DAmp, DLS_L_Gyr_z_DFreq, DLS_L_Gyr_z_PSD_mean, DLS_L_Gyr_z_PSD_std, DLS_L_Gyr_z_PSD_skew, DLS_L_Gyr_z_PSD_kurtosis, ...
    DLS_L_Gyr_norm_DAmp, DLS_L_Gyr_norm_DFreq, DLS_L_Gyr_norm_PSD_mean, DLS_L_Gyr_norm_PSD_std, DLS_L_Gyr_norm_PSD_skew, DLS_L_Gyr_norm_PSD_kurtosis, ...
    SC_Acc_x_mean, SC_Acc_y_mean, SC_Acc_z_mean, SC_Acc_norm_mean, ...
    DLS_R_Acc_x_mean, DLS_R_Acc_y_mean, DLS_R_Acc_z_mean, DLS_R_Acc_norm_mean, ...
    DLS_L_Acc_x_mean, DLS_L_Acc_y_mean, DLS_L_Acc_z_mean, DLS_L_Acc_norm_mean, ...
    SC_Acc_x_range, SC_Acc_y_range, SC_Acc_z_range, SC_Acc_norm_range, ...
    DLS_R_Acc_x_range, DLS_R_Acc_y_range, DLS_R_Acc_z_range, DLS_R_Acc_norm_range, ...
    DLS_L_Acc_x_range, DLS_L_Acc_y_range, DLS_L_Acc_z_range, DLS_L_Acc_norm_range, ...
    SC_Acc_x_rms, SC_Acc_y_rms, SC_Acc_z_rms, SC_Acc_norm_rms, ...
    DLS_R_Acc_x_rms, DLS_R_Acc_y_rms, DLS_R_Acc_z_rms, DLS_R_Acc_norm_rms, ...
    DLS_L_Acc_x_rms, DLS_L_Acc_y_rms, DLS_L_Acc_z_rms, DLS_L_Acc_norm_rms, ...
    SC_Acc_x_std, SC_Acc_y_std, SC_Acc_z_std, SC_Acc_norm_std, ...
    DLS_R_Acc_x_std, DLS_R_Acc_y_std, DLS_R_Acc_z_std, DLS_R_Acc_norm_std, ...
    DLS_L_Acc_x_std, DLS_L_Acc_y_std, DLS_L_Acc_z_std, DLS_L_Acc_norm_std, ...
    SC_Acc_x_skew, SC_Acc_y_skew, SC_Acc_z_skew, SC_Acc_norm_skew, ...
    DLS_R_Acc_x_skew, DLS_R_Acc_y_skew, DLS_R_Acc_z_skew, DLS_R_Acc_norm_skew, ...
    DLS_L_Acc_x_skew, DLS_L_Acc_y_skew, DLS_L_Acc_z_skew, DLS_L_Acc_norm_skew, ...
    SC_Acc_x_kurtosis, SC_Acc_y_kurtosis, SC_Acc_z_kurtosis, SC_Acc_norm_kurtosis, ...
    DLS_R_Acc_x_kurtosis, DLS_R_Acc_y_kurtosis, DLS_R_Acc_z_kurtosis, DLS_R_Acc_norm_kurtosis, ...
    DLS_L_Acc_x_kurtosis, DLS_L_Acc_y_kurtosis, DLS_L_Acc_z_kurtosis, DLS_L_Acc_norm_kurtosis, ...
    SC_Acc_corr_xy, SC_Acc_corr_xz, SC_Acc_corr_yz, ...
    DLS_R_Acc_corr_xy, DLS_R_Acc_corr_xz, DLS_R_Acc_corr_yz, ...
    DLS_L_Acc_corr_xy, DLS_L_Acc_corr_xz, DLS_L_Acc_corr_yz, ...
    SC_Acc_x_SamEn, SC_Acc_y_SamEn, SC_Acc_z_SamEn, SC_Acc_norm_SamEn, ...
    DLS_R_Acc_x_SamEn, DLS_R_Acc_y_SamEn, DLS_R_Acc_z_SamEn, DLS_R_Acc_norm_SamEn, ... 
    DLS_L_Acc_x_SamEn, DLS_L_Acc_y_SamEn, DLS_L_Acc_z_SamEn, DLS_L_Acc_norm_SamEn, ...
    SC_Acc_x_DAmp, SC_Acc_x_DFreq, SC_Acc_x_PSD_mean, SC_Acc_x_PSD_std, SC_Acc_x_PSD_skew, SC_Acc_x_PSD_kurtosis, ...
    SC_Acc_y_DAmp, SC_Acc_y_DFreq, SC_Acc_y_PSD_mean, SC_Acc_y_PSD_std, SC_Acc_y_PSD_skew, SC_Acc_y_PSD_kurtosis, ...
    SC_Acc_z_DAmp, SC_Acc_z_DFreq, SC_Acc_z_PSD_mean, SC_Acc_z_PSD_std, SC_Acc_z_PSD_skew, SC_Acc_z_PSD_kurtosis, ...
    SC_Acc_norm_DAmp, SC_Acc_norm_DFreq, SC_Acc_norm_PSD_mean, SC_Acc_norm_PSD_std, SC_Acc_norm_PSD_skew, SC_Acc_norm_PSD_kurtosis, ...        
    DLS_R_Acc_x_DAmp, DLS_R_Acc_x_DFreq, DLS_R_Acc_x_PSD_mean, DLS_R_Acc_x_PSD_std, DLS_R_Acc_x_PSD_skew, DLS_R_Acc_x_PSD_kurtosis, ...
    DLS_R_Acc_y_DAmp, DLS_R_Acc_y_DFreq, DLS_R_Acc_y_PSD_mean, DLS_R_Acc_y_PSD_std, DLS_R_Acc_y_PSD_skew, DLS_R_Acc_y_PSD_kurtosis, ...
    DLS_R_Acc_z_DAmp, DLS_R_Acc_z_DFreq, DLS_R_Acc_z_PSD_mean, DLS_R_Acc_z_PSD_std, DLS_R_Acc_z_PSD_skew, DLS_R_Acc_z_PSD_kurtosis, ...
    DLS_R_Acc_norm_DAmp, DLS_R_Acc_norm_DFreq, DLS_R_Acc_norm_PSD_mean, DLS_R_Acc_norm_PSD_std, DLS_R_Acc_norm_PSD_skew, DLS_R_Acc_norm_PSD_kurtosis, ...
    DLS_L_Acc_x_DAmp, DLS_L_Acc_x_DFreq, DLS_L_Acc_x_PSD_mean, DLS_L_Acc_x_PSD_std, DLS_L_Acc_x_PSD_skew, DLS_L_Acc_x_PSD_kurtosis, ...
    DLS_L_Acc_y_DAmp, DLS_L_Acc_y_DFreq, DLS_L_Acc_y_PSD_mean, DLS_L_Acc_y_PSD_std, DLS_L_Acc_y_PSD_skew, DLS_L_Acc_y_PSD_kurtosis, ...
    DLS_L_Acc_z_DAmp, DLS_L_Acc_z_DFreq, DLS_L_Acc_z_PSD_mean, DLS_L_Acc_z_PSD_std, DLS_L_Acc_z_PSD_skew, DLS_L_Acc_z_PSD_kurtosis, ...
    DLS_L_Acc_norm_DAmp, DLS_L_Acc_norm_DFreq, DLS_L_Acc_norm_PSD_mean, DLS_L_Acc_norm_PSD_std, DLS_L_Acc_norm_PSD_skew, DLS_L_Acc_norm_PSD_kurtosis)


writetable(tbl_Features,'Feature_Matrix_Automation_6MWT.csv','Delimiter',',','QuoteStrings',true)

