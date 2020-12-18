clc
clear all
close all

% CVA
% ID: 1, 3, 22, 25, 30, 32, 42, 44, 49  --> no data
% ID: 7, 40 --> sensor at DLS_L missing
% ID: 8 --> sensor at DLS_R missing
% ID: 35, 43 --> Short data 5595 (179sec, 43), 10167 (325sec, 35)  frames 

% ID: 26, 54 --> Outlier. Sensor data broken % code fixed


ID = [1:55];
Type_of_Subject = 'CVA'

Hz = 31.25;
dt = 1/Hz; 
k=1;

tbl = [];
TS = [5 10 20 30 60 90 120 180 240 300 360]
for k = 1:1:length(TS)
    for n = 1:1:length(ID)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID(n)) '.mat']
        load(file_input);

        Side = get_side('./CVA_Paretic_Side.csv');
        AS = Side.Side(ID(n));

        if n == 1 || n == 3 || n == 7 || n == 8 || n == 22 || n == 25 || n == 30 || n == 32 || n == 35 || n == 40 || n == 42 || n == 43 || n == 44 || n == 49
             SN = nan;     % Session #
             TN = nan;     % Trial #

            % Amount of Motion
            % Pelvis
            AoM_Pel_tilt(n,:) = nan;
            AoM_Pel_ro(n,:) = nan;
            AoM_Pel_oblq(n,:) = nan;
            AoM_Pel_norm(n,:) = nan;

            % Ankle sensor
            AoM_Ankle_AS_x(n,:) = nan;
            AoM_Ankle_US_x(n,:) = nan;
            AoM_Ankle_AS_y(n,:) = nan;
            AoM_Ankle_US_y(n,:) = nan;
            AoM_Ankle_AS_z(n,:) = nan;
            AoM_Ankle_US_z(n,:) = nan;
            AoM_Ankle_AS_norm(n,:) = nan;
            AoM_Ankle_US_norm(n,:) = nan;

            % General Features
            % Gyro mean
            SC_Gyr_x_mean(n,:) = nan;
            SC_Gyr_y_mean(n,:) = nan;
            SC_Gyr_z_mean(n,:) = nan;
            SC_Gyr_norm_mean(n,:) = nan;
            DLS_R_Gyr_x_mean(n,:) = nan;
            DLS_R_Gyr_y_mean(n,:) = nan;
            DLS_R_Gyr_z_mean(n,:) = nan;
            DLS_R_Gyr_norm_mean(n,:) = nan;
            DLS_L_Gyr_x_mean(n,:) = nan;
            DLS_L_Gyr_y_mean(n,:) = nan;
            DLS_L_Gyr_z_mean(n,:) = nan;
            DLS_L_Gyr_norm_mean(n,:) = nan;

            % Range
            SC_Gyr_x_range(n,:) = nan;
            SC_Gyr_y_range(n,:) = nan;
            SC_Gyr_z_range(n,:) = nan;
            SC_Gyr_norm_range(n,:) = nan;
            DLS_R_Gyr_x_range(n,:) = nan;
            DLS_R_Gyr_y_range(n,:) = nan;
            DLS_R_Gyr_z_range(n,:) = nan;
            DLS_R_Gyr_norm_range(n,:) = nan;
            DLS_L_Gyr_x_range(n,:) = nan;
            DLS_L_Gyr_y_range(n,:) = nan;
            DLS_L_Gyr_z_range(n,:) = nan;
            DLS_L_Gyr_norm_range(n,:) = nan;

            % RMS
            SC_Gyr_x_rms(n,:) = nan;
            SC_Gyr_y_rms(n,:) = nan;
            SC_Gyr_z_rms(n,:) = nan;
            SC_Gyr_norm_rms(n,:) = nan;
            DLS_R_Gyr_x_rms(n,:) = nan;
            DLS_R_Gyr_y_rms(n,:) = nan;
            DLS_R_Gyr_z_rms(n,:) = nan;
            DLS_R_Gyr_norm_rms(n,:) = nan;
            DLS_L_Gyr_x_rms(n,:) = nan;
            DLS_L_Gyr_y_rms(n,:) = nan;
            DLS_L_Gyr_z_rms(n,:) = nan;
            DLS_L_Gyr_norm_rms(n,:) = nan;

            % Standard Deviation
            SC_Gyr_x_std(n,:) = nan;
            SC_Gyr_y_std(n,:) = nan;
            SC_Gyr_z_std(n,:) = nan;
            SC_Gyr_norm_std(n,:) = nan;
            DLS_R_Gyr_x_std(n,:) = nan;
            DLS_R_Gyr_y_std(n,:) = nan;
            DLS_R_Gyr_z_std(n,:) = nan;
            DLS_R_Gyr_norm_std(n,:) = nan;
            DLS_L_Gyr_x_std(n,:) = nan;
            DLS_L_Gyr_y_std(n,:) = nan;
            DLS_L_Gyr_z_std(n,:) = nan;
            DLS_L_Gyr_norm_std(n,:) = nan;

            % Skew
            SC_Gyr_x_skew(n,:) = nan;
            SC_Gyr_y_skew(n,:) = nan;
            SC_Gyr_z_skew(n,:) = nan;
            SC_Gyr_norm_skew(n,:) = nan;
            DLS_R_Gyr_x_skew(n,:) = nan;
            DLS_R_Gyr_y_skew(n,:) = nan;
            DLS_R_Gyr_z_skew(n,:) = nan;
            DLS_R_Gyr_norm_skew(n,:) = nan;
            DLS_L_Gyr_x_skew(n,:) = nan;
            DLS_L_Gyr_y_skew(n,:) = nan;
            DLS_L_Gyr_z_skew(n,:) = nan;
            DLS_L_Gyr_norm_skew(n,:) = nan;

            % Kurtosis
            SC_Gyr_x_kurtosis(n,:) = nan;
            SC_Gyr_y_kurtosis(n,:) = nan;
            SC_Gyr_z_kurtosis(n,:) = nan;
            SC_Gyr_norm_kurtosis(n,:) = nan;
            DLS_R_Gyr_x_kurtosis(n,:) = nan;
            DLS_R_Gyr_y_kurtosis(n,:) = nan;
            DLS_R_Gyr_z_kurtosis(n,:) = nan;
            DLS_R_Gyr_norm_kurtosis(n,:) = nan;
            DLS_L_Gyr_x_kurtosis(n,:) = nan;
            DLS_L_Gyr_y_kurtosis(n,:) = nan;
            DLS_L_Gyr_z_kurtosis(n,:) = nan;
            DLS_L_Gyr_norm_kurtosis(n,:) = nan;

            % Pearson correlation coefficient
            SC_Gyr_corr_xy(n,:) = nan;
            SC_Gyr_corr_xz(n,:) = nan;
            SC_Gyr_corr_yz(n,:) = nan;
            DLS_R_Gyr_corr_xy(n,:) = nan;
            DLS_R_Gyr_corr_xz(n,:) = nan;
            DLS_R_Gyr_corr_yz(n,:) = nan;
            DLS_L_Gyr_corr_xy(n,:) = nan;
            DLS_L_Gyr_corr_xz(n,:) = nan;
            DLS_L_Gyr_corr_yz(n,:) = nan;

            % Sample Entropy
            SC_Gyr_x_SamEn(n,:) = nan;
            SC_Gyr_y_SamEn(n,:) = nan;
            SC_Gyr_z_SamEn(n,:) = nan;
            SC_Gyr_norm_SamEn(n,:) = nan;
            DLS_R_Gyr_x_SamEn(n,:) = nan;
            DLS_R_Gyr_y_SamEn(n,:) = nan;
            DLS_R_Gyr_z_SamEn(n,:) = nan;
            DLS_R_Gyr_norm_SamEn(n,:) = nan;
            DLS_L_Gyr_x_SamEn(n,:) = nan;
            DLS_L_Gyr_y_SamEn(n,:) = nan;
            DLS_L_Gyr_z_SamEn(n,:) = nan;
            DLS_L_Gyr_norm_SamEn(n,:) = nan;

            % Frequency Domain
            SC_Gyr_x_DAmp(n,:) = nan;
            SC_Gyr_x_DFreq(n,:) = nan;
            SC_Gyr_x_PSD_mean(n,:) = nan;
            SC_Gyr_x_PSD_std(n,:) = nan;
            SC_Gyr_x_PSD_skew(n,:) = nan;
            SC_Gyr_x_PSD_kurtosis(n,:) = nan;

            SC_Gyr_y_DAmp(n,:) = nan;
            SC_Gyr_y_DFreq(n,:) = nan;
            SC_Gyr_y_PSD_mean(n,:) = nan;
            SC_Gyr_y_PSD_std(n,:) = nan;
            SC_Gyr_y_PSD_skew(n,:) = nan;
            SC_Gyr_y_PSD_kurtosis(n,:) = nan;

            SC_Gyr_z_DAmp(n,:) = nan;
            SC_Gyr_z_DFreq(n,:) = nan;
            SC_Gyr_z_PSD_mean(n,:) = nan;
            SC_Gyr_z_PSD_std(n,:) = nan;
            SC_Gyr_z_PSD_skew(n,:) = nan;
            SC_Gyr_z_PSD_kurtosis(n,:) = nan;

            SC_Gyr_norm_DAmp(n,:) = nan;
            SC_Gyr_norm_DFreq(n,:) = nan;
            SC_Gyr_norm_PSD_mean(n,:) = nan;
            SC_Gyr_norm_PSD_std(n,:) = nan;
            SC_Gyr_norm_PSD_skew(n,:) = nan;
            SC_Gyr_norm_PSD_kurtosis(n,:) = nan;

            DLS_R_Gyr_x_DAmp(n,:) = nan;
            DLS_R_Gyr_x_DFreq(n,:) = nan;
            DLS_R_Gyr_x_PSD_mean(n,:) = nan;
            DLS_R_Gyr_x_PSD_std(n,:) = nan;
            DLS_R_Gyr_x_PSD_skew(n,:) = nan;
            DLS_R_Gyr_x_PSD_kurtosis(n,:) = nan;

            DLS_R_Gyr_y_DAmp(n,:) = nan;
            DLS_R_Gyr_y_DFreq(n,:) = nan;
            DLS_R_Gyr_y_PSD_mean(n,:) = nan;
            DLS_R_Gyr_y_PSD_std(n,:) = nan;
            DLS_R_Gyr_y_PSD_skew(n,:) = nan;
            DLS_R_Gyr_y_PSD_kurtosis(n,:) = nan;

            DLS_R_Gyr_z_DAmp(n,:) = nan;
            DLS_R_Gyr_z_DFreq(n,:) = nan;
            DLS_R_Gyr_z_PSD_mean(n,:) = nan;
            DLS_R_Gyr_z_PSD_std(n,:) = nan;
            DLS_R_Gyr_z_PSD_skew(n,:) = nan;
            DLS_R_Gyr_z_PSD_kurtosis(n,:) = nan;

            DLS_R_Gyr_norm_DAmp(n,:) = nan;
            DLS_R_Gyr_norm_DFreq(n,:) = nan;
            DLS_R_Gyr_norm_PSD_mean(n,:) = nan;
            DLS_R_Gyr_norm_PSD_std(n,:) = nan;
            DLS_R_Gyr_norm_PSD_skew(n,:) = nan;
            DLS_R_Gyr_norm_PSD_kurtosis(n,:) = nan;

            DLS_L_Gyr_x_DAmp(n,:) = nan;
            DLS_L_Gyr_x_DFreq(n,:) = nan;
            DLS_L_Gyr_x_PSD_mean(n,:) = nan;
            DLS_L_Gyr_x_PSD_std(n,:) = nan;
            DLS_L_Gyr_x_PSD_skew(n,:) = nan;
            DLS_L_Gyr_x_PSD_kurtosis(n,:) = nan;

            DLS_L_Gyr_y_DAmp(n,:) = nan;
            DLS_L_Gyr_y_DFreq(n,:) = nan;
            DLS_L_Gyr_y_PSD_mean(n,:) = nan;
            DLS_L_Gyr_y_PSD_std(n,:) = nan;
            DLS_L_Gyr_y_PSD_skew(n,:) = nan;
            DLS_L_Gyr_y_PSD_kurtosis(n,:) = nan;

            DLS_L_Gyr_z_DAmp(n,:) = nan;
            DLS_L_Gyr_z_DFreq(n,:) = nan;
            DLS_L_Gyr_z_PSD_mean(n,:) = nan;
            DLS_L_Gyr_z_PSD_std(n,:) = nan;
            DLS_L_Gyr_z_PSD_skew(n,:) = nan;
            DLS_L_Gyr_z_PSD_kurtosis(n,:) = nan;

            DLS_L_Gyr_norm_DAmp(n,:) = nan;
            DLS_L_Gyr_norm_DFreq(n,:) = nan;
            DLS_L_Gyr_norm_PSD_mean(n,:) = nan;
            DLS_L_Gyr_norm_PSD_std(n,:) = nan;
            DLS_L_Gyr_norm_PSD_skew(n,:) = nan;
            DLS_L_Gyr_norm_PSD_kurtosis(n,:) = nan;



            % Acc mean
            SC_Acc_x_mean(n,:) = nan;
            SC_Acc_y_mean(n,:) = nan;
            SC_Acc_z_mean(n,:) = nan;
            SC_Acc_norm_mean(n,:) = nan;
            DLS_R_Acc_x_mean(n,:) = nan;
            DLS_R_Acc_y_mean(n,:) = nan;
            DLS_R_Acc_z_mean(n,:) = nan;
            DLS_R_Acc_norm_mean(n,:) = nan;
            DLS_L_Acc_x_mean(n,:) = nan;
            DLS_L_Acc_y_mean(n,:) = nan;
            DLS_L_Acc_z_mean(n,:) = nan;
            DLS_L_Acc_norm_mean(n,:) = nan;

            % Range
            SC_Acc_x_range(n,:) = nan;
            SC_Acc_y_range(n,:) = nan;
            SC_Acc_z_range(n,:) = nan;
            SC_Acc_norm_range(n,:) = nan;
            DLS_R_Acc_x_range(n,:) = nan;
            DLS_R_Acc_y_range(n,:) = nan;
            DLS_R_Acc_z_range(n,:) = nan;
            DLS_R_Acc_norm_range(n,:) = nan;
            DLS_L_Acc_x_range(n,:) = nan;
            DLS_L_Acc_y_range(n,:) = nan;
            DLS_L_Acc_z_range(n,:) = nan;
            DLS_L_Acc_norm_range(n,:) = nan;

            % RMS
            SC_Acc_x_rms(n,:) = nan;
            SC_Acc_y_rms(n,:) = nan;
            SC_Acc_z_rms(n,:) = nan;
            SC_Acc_norm_rms(n,:) = nan;
            DLS_R_Acc_x_rms(n,:) = nan;
            DLS_R_Acc_y_rms(n,:) = nan;
            DLS_R_Acc_z_rms(n,:) = nan;
            DLS_R_Acc_norm_rms(n,:) = nan;
            DLS_L_Acc_x_rms(n,:) = nan;
            DLS_L_Acc_y_rms(n,:) = nan;
            DLS_L_Acc_z_rms(n,:) = nan;
            DLS_L_Acc_norm_rms(n,:) = nan;

            % Standard Deviation
            SC_Acc_x_std(n,:) = nan;
            SC_Acc_y_std(n,:) = nan;
            SC_Acc_z_std(n,:) = nan;
            SC_Acc_norm_std(n,:) = nan;
            DLS_R_Acc_x_std(n,:) = nan;
            DLS_R_Acc_y_std(n,:) = nan;
            DLS_R_Acc_z_std(n,:) = nan;
            DLS_R_Acc_norm_std(n,:) = nan;
            DLS_L_Acc_x_std(n,:) = nan;
            DLS_L_Acc_y_std(n,:) = nan;
            DLS_L_Acc_z_std(n,:) = nan;
            DLS_L_Acc_norm_std(n,:) = nan;

            % Skew
            SC_Acc_x_skew(n,:) = nan;
            SC_Acc_y_skew(n,:) = nan;
            SC_Acc_z_skew(n,:) = nan;
            SC_Acc_norm_skew(n,:) = nan;
            DLS_R_Acc_x_skew(n,:) = nan;
            DLS_R_Acc_y_skew(n,:) = nan;
            DLS_R_Acc_z_skew(n,:) = nan;
            DLS_R_Acc_norm_skew(n,:) = nan;
            DLS_L_Acc_x_skew(n,:) = nan;
            DLS_L_Acc_y_skew(n,:) = nan;
            DLS_L_Acc_z_skew(n,:) = nan;
            DLS_L_Acc_norm_skew(n,:) = nan;

            % Kurtosis
            SC_Acc_x_kurtosis(n,:) = nan;
            SC_Acc_y_kurtosis(n,:) = nan;
            SC_Acc_z_kurtosis(n,:) = nan;
            SC_Acc_norm_kurtosis(n,:) = nan;
            DLS_R_Acc_x_kurtosis(n,:) = nan;
            DLS_R_Acc_y_kurtosis(n,:) = nan;
            DLS_R_Acc_z_kurtosis(n,:) = nan;
            DLS_R_Acc_norm_kurtosis(n,:) = nan;
            DLS_L_Acc_x_kurtosis(n,:) = nan;
            DLS_L_Acc_y_kurtosis(n,:) = nan;
            DLS_L_Acc_z_kurtosis(n,:) = nan;
            DLS_L_Acc_norm_kurtosis(n,:) = nan;

            % Pearson correlation coefficient
            SC_Acc_corr_xy(n,:) = nan;
            SC_Acc_corr_xz(n,:) = nan;
            SC_Acc_corr_yz(n,:) = nan;
            DLS_R_Acc_corr_xy(n,:) = nan;
            DLS_R_Acc_corr_xz(n,:) = nan;
            DLS_R_Acc_corr_yz(n,:) = nan;
            DLS_L_Acc_corr_xy(n,:) = nan;
            DLS_L_Acc_corr_xz(n,:) = nan;
            DLS_L_Acc_corr_yz(n,:) = nan;

            % Sample Entropy
            SC_Acc_x_SamEn(n,:) = nan;
            SC_Acc_y_SamEn(n,:) = nan;
            SC_Acc_z_SamEn(n,:) = nan;
            SC_Acc_norm_SamEn(n,:) = nan;
            DLS_R_Acc_x_SamEn(n,:) = nan;
            DLS_R_Acc_y_SamEn(n,:) = nan;
            DLS_R_Acc_z_SamEn(n,:) = nan;
            DLS_R_Acc_norm_SamEn(n,:) = nan;
            DLS_L_Acc_x_SamEn(n,:) = nan;
            DLS_L_Acc_y_SamEn(n,:) = nan;
            DLS_L_Acc_z_SamEn(n,:) = nan;
            DLS_L_Acc_norm_SamEn(n,:) = nan;

            % Frequency Domain
            SC_Acc_x_DAmp(n,:) = nan;
            SC_Acc_x_DFreq(n,:) = nan;
            SC_Acc_x_PSD_mean(n,:) = nan;
            SC_Acc_x_PSD_std(n,:) = nan;
            SC_Acc_x_PSD_skew(n,:) = nan;
            SC_Acc_x_PSD_kurtosis(n,:) = nan;

            SC_Acc_y_DAmp(n,:) = nan;
            SC_Acc_y_DFreq(n,:) = nan;
            SC_Acc_y_PSD_mean(n,:) = nan;
            SC_Acc_y_PSD_std(n,:) = nan;
            SC_Acc_y_PSD_skew(n,:) = nan;
            SC_Acc_y_PSD_kurtosis(n,:) = nan;

            SC_Acc_z_DAmp(n,:) = nan;
            SC_Acc_z_DFreq(n,:) = nan;
            SC_Acc_z_PSD_mean(n,:) = nan;
            SC_Acc_z_PSD_std(n,:) = nan;
            SC_Acc_z_PSD_skew(n,:) = nan;
            SC_Acc_z_PSD_kurtosis(n,:) = nan;

            SC_Acc_norm_DAmp(n,:) = nan;
            SC_Acc_norm_DFreq(n,:) = nan;
            SC_Acc_norm_PSD_mean(n,:) = nan;
            SC_Acc_norm_PSD_std(n,:) = nan;
            SC_Acc_norm_PSD_skew(n,:) = nan;
            SC_Acc_norm_PSD_kurtosis(n,:) = nan;

            DLS_R_Acc_x_DAmp(n,:) = nan;
            DLS_R_Acc_x_DFreq(n,:) = nan;
            DLS_R_Acc_x_PSD_mean(n,:) = nan;
            DLS_R_Acc_x_PSD_std(n,:) = nan;
            DLS_R_Acc_x_PSD_skew(n,:) = nan;
            DLS_R_Acc_x_PSD_kurtosis(n,:) = nan;

            DLS_R_Acc_y_DAmp(n,:) = nan;
            DLS_R_Acc_y_DFreq(n,:) = nan;
            DLS_R_Acc_y_PSD_mean(n,:) = nan;
            DLS_R_Acc_y_PSD_std(n,:) = nan;
            DLS_R_Acc_y_PSD_skew(n,:) = nan;
            DLS_R_Acc_y_PSD_kurtosis(n,:) = nan;

            DLS_R_Acc_z_DAmp(n,:) = nan;
            DLS_R_Acc_z_DFreq(n,:) = nan;
            DLS_R_Acc_z_PSD_mean(n,:) = nan;
            DLS_R_Acc_z_PSD_std(n,:) = nan;
            DLS_R_Acc_z_PSD_skew(n,:) = nan;
            DLS_R_Acc_z_PSD_kurtosis(n,:) = nan;

            DLS_R_Acc_norm_DAmp(n,:) = nan;
            DLS_R_Acc_norm_DFreq(n,:) = nan;
            DLS_R_Acc_norm_PSD_mean(n,:) = nan;
            DLS_R_Acc_norm_PSD_std(n,:) = nan;
            DLS_R_Acc_norm_PSD_skew(n,:) = nan;
            DLS_R_Acc_norm_PSD_kurtosis(n,:) = nan;

            DLS_L_Acc_x_DAmp(n,:) = nan;
            DLS_L_Acc_x_DFreq(n,:) = nan;
            DLS_L_Acc_x_PSD_mean(n,:) = nan;
            DLS_L_Acc_x_PSD_std(n,:) = nan;
            DLS_L_Acc_x_PSD_skew(n,:) = nan;
            DLS_L_Acc_x_PSD_kurtosis(n,:) = nan;

            DLS_L_Acc_y_DAmp(n,:) = nan;
            DLS_L_Acc_y_DFreq(n,:) = nan;
            DLS_L_Acc_y_PSD_mean(n,:) = nan;
            DLS_L_Acc_y_PSD_std(n,:) = nan;
            DLS_L_Acc_y_PSD_skew(n,:) = nan;
            DLS_L_Acc_y_PSD_kurtosis(n,:) = nan;

            DLS_L_Acc_z_DAmp(n,:) = nan;
            DLS_L_Acc_z_DFreq(n,:) = nan;
            DLS_L_Acc_z_PSD_mean(n,:) = nan;
            DLS_L_Acc_z_PSD_std(n,:) = nan;
            DLS_L_Acc_z_PSD_skew(n,:) = nan;
            DLS_L_Acc_z_PSD_kurtosis(n,:) = nan;

            DLS_L_Acc_norm_DAmp(n,:) = nan;
            DLS_L_Acc_norm_DFreq(n,:) = nan;
            DLS_L_Acc_norm_PSD_mean(n,:) = nan;
            DLS_L_Acc_norm_PSD_std(n,:) = nan;
            DLS_L_Acc_norm_PSD_skew(n,:) = nan;
            DLS_L_Acc_norm_PSD_kurtosis(n,:) = nan;

        else
            SN = 1;
            TN = 1;
            
            if n == 26 || n == 54
                SC_Gyr = [data.Session{SN}.Motion.SC.Gyr{1}; data.Session{SN}.Motion.SC.Gyr{2}];
                DLS_R_Gyr = [data.Session{SN}.Motion.DLS_R.Gyr{1}; data.Session{SN}.Motion.DLS_R.Gyr{2}];
                DLS_L_Gyr = [data.Session{SN}.Motion.DLS_L.Gyr{1}; data.Session{SN}.Motion.DLS_L.Gyr{2}];

                % Acceleration data
                SC_Acc = [data.Session{SN}.Motion.SC.Acc{1}; data.Session{SN}.Motion.SC.Acc{2}];
                DLS_R_Acc = [data.Session{SN}.Motion.DLS_R.Acc{1}; data.Session{SN}.Motion.DLS_R.Acc{2}]
                DLS_L_Acc = [data.Session{SN}.Motion.DLS_L.Acc{1}; data.Session{SN}.Motion.DLS_L.Acc{2}]
                
                 % To match final time
    %              final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
                final = Hz*TS(k);  
                Time_ = [data.Session{SN}.Motion.Time{1}; data.Session{SN}.Motion.Time{2}]
                Time = Time_(1:final,:);
                
                SC_Gyr = SC_Gyr(1:final,:);
                PLS_R_Gyr = DLS_R_Gyr(1:final,:);
                DLS_L_Gyr = DLS_L_Gyr(1:final,:);

                SC_Acc = SC_Acc(1:final,:);
                DLS_R_Acc = DLS_R_Acc(1:final,:);
                DLS_L_Acc = DLS_L_Acc(1:final,:);
                
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
                
                SC_Acc = SC_Acc(1:final,:);
                DLS_R_Acc = DLS_R_Acc(1:final,:);
                DLS_L_Acc = DLS_L_Acc(1:final,:);
                
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

            % Amount of motion
            % Pelvic Amount of motion 
            AoM_Pel_tilt(n,:) = sum(abs(SC_Gyr(:,1))) * dt;
            AoM_Pel_ro(n,:) = sum(abs(SC_Gyr(:,2))) * dt;
            AoM_Pel_oblq(n,:) = sum(abs(SC_Gyr(:,3))) * dt;
            AoM_Pel_norm(n,:) = sum(abs(SC_Gyr_norm)) * dt;     

           % Ankle Amount of motion
            if strcmp(AS,'L') == 1
                AoM_Ankle_AS_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;  
                AoM_Ankle_US_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;
                AoM_Ankle_AS_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;  
                AoM_Ankle_US_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;
                AoM_Ankle_AS_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;  
                AoM_Ankle_US_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;
                AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt; 
                AoM_Ankle_US_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt;

            elseif strcmp(AS,'R') == 1
                AoM_Ankle_AS_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;  
                AoM_Ankle_US_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;
                AoM_Ankle_AS_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;  
                AoM_Ankle_US_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;
                AoM_Ankle_AS_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;  
                AoM_Ankle_US_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;
                AoM_Ankle_AS_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt; 
                AoM_Ankle_US_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt;
            end  

            % General Features
            % Gyro mean
            SC_Gyr_x_mean(n,:) = mean(SC_Gyr(:,1));
            SC_Gyr_y_mean(n,:) = mean(SC_Gyr(:,2));
            SC_Gyr_z_mean(n,:) = mean(SC_Gyr(:,3));
            SC_Gyr_norm_mean(n,:) = mean(SC_Gyr_norm);
            DLS_R_Gyr_x_mean(n,:) = mean(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_mean(n,:) = mean(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_mean(n,:) = mean(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_mean(n,:) = mean(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_mean(n,:) = mean(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_mean(n,:) = mean(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_mean(n,:) = mean(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_mean(n,:) = mean(DLS_L_Gyr_norm);

            % Range
            SC_Gyr_x_range(n,:) = range(SC_Gyr(:,1));
            SC_Gyr_y_range(n,:) = range(SC_Gyr(:,2));
            SC_Gyr_z_range(n,:) = range(SC_Gyr(:,3));
            SC_Gyr_norm_range(n,:) = range(SC_Gyr_norm);
            DLS_R_Gyr_x_range(n,:) = range(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_range(n,:) = range(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_range(n,:) = range(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_range(n,:) = range(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_range(n,:) = range(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_range(n,:) = range(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_range(n,:) = range(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_range(n,:) = range(DLS_L_Gyr_norm);

            % RMS
            SC_Gyr_x_rms(n,:) = rms(SC_Gyr(:,1));
            SC_Gyr_y_rms(n,:) = rms(SC_Gyr(:,2));
            SC_Gyr_z_rms(n,:) = rms(SC_Gyr(:,3));
            SC_Gyr_norm_rms(n,:) = rms(SC_Gyr_norm);
            DLS_R_Gyr_x_rms(n,:) = rms(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_rms(n,:) = rms(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_rms(n,:) = rms(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_rms(n,:) = rms(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_rms(n,:) = rms(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_rms(n,:) = rms(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_rms(n,:) = rms(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_rms(n,:) = rms(DLS_L_Gyr_norm);

            % Standard Deviation
            SC_Gyr_x_std(n,:) = std(SC_Gyr(:,1));
            SC_Gyr_y_std(n,:) = std(SC_Gyr(:,2));
            SC_Gyr_z_std(n,:) = std(SC_Gyr(:,3));
            SC_Gyr_norm_std(n,:) = std(SC_Gyr_norm);
            DLS_R_Gyr_x_std(n,:) = std(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_std(n,:) = std(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_std(n,:) = std(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_std(n,:) = std(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_std(n,:) = std(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_std(n,:) = std(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_std(n,:) = std(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_std(n,:) = std(DLS_L_Gyr_norm);

            % Skew
            SC_Gyr_x_skew(n,:) = skewness(SC_Gyr(:,1));
            SC_Gyr_y_skew(n,:) = skewness(SC_Gyr(:,2));
            SC_Gyr_z_skew(n,:) = skewness(SC_Gyr(:,3));
            SC_Gyr_norm_skew(n,:) = skewness(SC_Gyr_norm);
            DLS_R_Gyr_x_skew(n,:) = skewness(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_skew(n,:) = skewness(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_skew(n,:) = skewness(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_skew(n,:) = skewness(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_skew(n,:) = skewness(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_skew(n,:) = skewness(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_skew(n,:) = skewness(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_skew(n,:) = skewness(DLS_L_Gyr_norm);

            % Kurtosis
            SC_Gyr_x_kurtosis(n,:) = kurtosis(SC_Gyr(:,1));
            SC_Gyr_y_kurtosis(n,:) = kurtosis(SC_Gyr(:,2));
            SC_Gyr_z_kurtosis(n,:) = kurtosis(SC_Gyr(:,3));
            SC_Gyr_norm_kurtosis(n,:) = kurtosis(SC_Gyr_norm);
            DLS_R_Gyr_x_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,1));
            DLS_R_Gyr_y_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,2));
            DLS_R_Gyr_z_kurtosis(n,:) = kurtosis(DLS_R_Gyr(:,3));
            DLS_R_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_R_Gyr_norm);
            DLS_L_Gyr_x_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,1));
            DLS_L_Gyr_y_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,2));
            DLS_L_Gyr_z_kurtosis(n,:) = kurtosis(DLS_L_Gyr(:,3));
            DLS_L_Gyr_norm_kurtosis(n,:) = kurtosis(DLS_L_Gyr_norm);

            % Pearson correlation coefficient
            SC_Gyr_corr_xy(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,2));
            SC_Gyr_corr_xz(n,:) = corr(SC_Gyr(:,1),SC_Gyr(:,3));
            SC_Gyr_corr_yz(n,:) = corr(SC_Gyr(:,2),SC_Gyr(:,3));
            DLS_R_Gyr_corr_xy(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,2));
            DLS_R_Gyr_corr_xz(n,:) = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,3));
            DLS_R_Gyr_corr_yz(n,:) = corr(DLS_R_Gyr(:,2),DLS_R_Gyr(:,3));
            DLS_L_Gyr_corr_xy(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,2));
            DLS_L_Gyr_corr_xz(n,:) = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,3));
            DLS_L_Gyr_corr_yz(n,:) = corr(DLS_L_Gyr(:,2),DLS_L_Gyr(:,3));

            % Sample Entropy
            r = 0.2;
            SC_Gyr_x_SamEn(n,:) = sampen(SC_Gyr(:,1),1,r);
            SC_Gyr_y_SamEn(n,:) = sampen(SC_Gyr(:,2),1,r);
            SC_Gyr_z_SamEn(n,:) = sampen(SC_Gyr(:,3),1,r);
            SC_Gyr_norm_SamEn(n,:) = sampen(SC_Gyr_norm,1,r);
            DLS_R_Gyr_x_SamEn(n,:) = sampen(DLS_R_Gyr(:,1),1,r);
            DLS_R_Gyr_y_SamEn(n,:) = sampen(DLS_R_Gyr(:,2),1,r);
            DLS_R_Gyr_z_SamEn(n,:) = sampen(DLS_R_Gyr(:,3),1,r);
            DLS_R_Gyr_norm_SamEn(n,:) = sampen(DLS_R_Gyr_norm,1,r);
            DLS_L_Gyr_x_SamEn(n,:) = sampen(DLS_L_Gyr(:,1),1,r);
            DLS_L_Gyr_y_SamEn(n,:) = sampen(DLS_L_Gyr(:,2),1,r);
            DLS_L_Gyr_z_SamEn(n,:) = sampen(DLS_L_Gyr(:,3),1,r);
            DLS_L_Gyr_norm_SamEn(n,:) = sampen(DLS_L_Gyr_norm,1,r);
            
            
            % Frequency Domain
            ff = FFeatures(SC_Gyr(:,1), Hz);
            SC_Gyr_x_DAmp(n,:) = ff(1);
            SC_Gyr_x_DFreq(n,:) = ff(2);
            SC_Gyr_x_PSD_mean(n,:) = ff(3);
            SC_Gyr_x_PSD_std(n,:) = ff(4);
            SC_Gyr_x_PSD_skew(n,:) = ff(5);
            SC_Gyr_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(SC_Gyr(:,2), Hz);
            SC_Gyr_y_DAmp(n,:) = ff(1);
            SC_Gyr_y_DFreq(n,:) = ff(2);
            SC_Gyr_y_PSD_mean(n,:) = ff(3);
            SC_Gyr_y_PSD_std(n,:) = ff(4);
            SC_Gyr_y_PSD_skew(n,:) = ff(5);
            SC_Gyr_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(SC_Gyr(:,3), Hz);
            SC_Gyr_z_DAmp(n,:) = ff(1);
            SC_Gyr_z_DFreq(n,:) = ff(2);
            SC_Gyr_z_PSD_mean(n,:) = ff(3);
            SC_Gyr_z_PSD_std(n,:) = ff(4);
            SC_Gyr_z_PSD_skew(n,:) = ff(5);
            SC_Gyr_z_PSD_kurtosis(n,:) =ff(6);

            ff = FFeatures(SC_Gyr_norm, Hz);
            SC_Gyr_norm_DAmp(n,:) = ff(1);
            SC_Gyr_norm_DFreq(n,:) = ff(2);
            SC_Gyr_norm_PSD_mean(n,:) = ff(3);
            SC_Gyr_norm_PSD_std(n,:) = ff(4);
            SC_Gyr_norm_PSD_skew(n,:) = ff(5);
            SC_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Gyr(:,1), Hz);
            DLS_R_Gyr_x_DAmp(n,:) = ff(1);
            DLS_R_Gyr_x_DFreq(n,:) = ff(2);
            DLS_R_Gyr_x_PSD_mean(n,:) = ff(3);
            DLS_R_Gyr_x_PSD_std(n,:) = ff(4);
            DLS_R_Gyr_x_PSD_skew(n,:) = ff(5);
            DLS_R_Gyr_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Gyr(:,2), Hz);
            DLS_R_Gyr_y_DAmp(n,:) = ff(1);
            DLS_R_Gyr_y_DFreq(n,:) = ff(2);
            DLS_R_Gyr_y_PSD_mean(n,:) = ff(3);
            DLS_R_Gyr_y_PSD_std(n,:) = ff(4);
            DLS_R_Gyr_y_PSD_skew(n,:) = ff(5);
            DLS_R_Gyr_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Gyr(:,3), Hz);
            DLS_R_Gyr_z_DAmp(n,:) = ff(1);
            DLS_R_Gyr_z_DFreq(n,:) = ff(2);
            DLS_R_Gyr_z_PSD_mean(n,:) = ff(3);
            DLS_R_Gyr_z_PSD_std(n,:) = ff(4);
            DLS_R_Gyr_z_PSD_skew(n,:) = ff(5);
            DLS_R_Gyr_z_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Gyr_norm, Hz);
            DLS_R_Gyr_norm_DAmp(n,:) = ff(1);
            DLS_R_Gyr_norm_DFreq(n,:) = ff(2);
            DLS_R_Gyr_norm_PSD_mean(n,:) = ff(3);
            DLS_R_Gyr_norm_PSD_std(n,:) = ff(4);
            DLS_R_Gyr_norm_PSD_skew(n,:) = ff(5);
            DLS_R_Gyr_norm_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Gyr(:,1), Hz);
            DLS_L_Gyr_x_DAmp(n,:) = ff(1);
            DLS_L_Gyr_x_DFreq(n,:) = ff(2);
            DLS_L_Gyr_x_PSD_mean(n,:) = ff(3);
            DLS_L_Gyr_x_PSD_std(n,:) = ff(4);
            DLS_L_Gyr_x_PSD_skew(n,:) = ff(5);
            DLS_L_Gyr_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Gyr(:,2), Hz);
            DLS_L_Gyr_y_DAmp(n,:) = ff(1);
            DLS_L_Gyr_y_DFreq(n,:) = ff(2);
            DLS_L_Gyr_y_PSD_mean(n,:) = ff(3);
            DLS_L_Gyr_y_PSD_std(n,:) = ff(4);
            DLS_L_Gyr_y_PSD_skew(n,:) = ff(5);
            DLS_L_Gyr_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Gyr(:,3), Hz);
            DLS_L_Gyr_z_DAmp(n,:) = ff(1);
            DLS_L_Gyr_z_DFreq(n,:) = ff(2);
            DLS_L_Gyr_z_PSD_mean(n,:) = ff(3);
            DLS_L_Gyr_z_PSD_std(n,:) = ff(4);
            DLS_L_Gyr_z_PSD_skew(n,:) = ff(5);
            DLS_L_Gyr_z_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Gyr_norm, Hz);
            DLS_L_Gyr_norm_DAmp(n,:) = ff(1);
            DLS_L_Gyr_norm_DFreq(n,:) = ff(2);
            DLS_L_Gyr_norm_PSD_mean(n,:) = ff(3);
            DLS_L_Gyr_norm_PSD_std(n,:) = ff(4);
            DLS_L_Gyr_norm_PSD_skew(n,:) = ff(5);
            DLS_L_Gyr_norm_PSD_kurtosis(n,:) = ff(6);


            % Acc mean
            SC_Acc_x_mean(n,:) = mean(SC_Acc(:,1));
            SC_Acc_y_mean(n,:) = mean(SC_Acc(:,2));
            SC_Acc_z_mean(n,:) = mean(SC_Acc(:,3));
            SC_Acc_norm_mean(n,:) = mean(SC_Acc_norm);
            DLS_R_Acc_x_mean(n,:) = mean(DLS_R_Acc(:,1));
            DLS_R_Acc_y_mean(n,:) = mean(DLS_R_Acc(:,2));
            DLS_R_Acc_z_mean(n,:) = mean(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_mean(n,:) = mean(DLS_R_Acc_norm);
            DLS_L_Acc_x_mean(n,:) = mean(DLS_L_Acc(:,1));
            DLS_L_Acc_y_mean(n,:) = mean(DLS_L_Acc(:,2));
            DLS_L_Acc_z_mean(n,:) = mean(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_mean(n,:) = mean(DLS_L_Acc_norm);

            % Range
            SC_Acc_x_range(n,:) = range(SC_Acc(:,1));
            SC_Acc_y_range(n,:) = range(SC_Acc(:,2));
            SC_Acc_z_range(n,:) = range(SC_Acc(:,3));
            SC_Acc_norm_range(n,:) = range(SC_Acc_norm);
            DLS_R_Acc_x_range(n,:) = range(DLS_R_Acc(:,1));
            DLS_R_Acc_y_range(n,:) = range(DLS_R_Acc(:,2));
            DLS_R_Acc_z_range(n,:) = range(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_range(n,:) = range(DLS_R_Acc_norm);
            DLS_L_Acc_x_range(n,:) = range(DLS_L_Acc(:,1));
            DLS_L_Acc_y_range(n,:) = range(DLS_L_Acc(:,2));
            DLS_L_Acc_z_range(n,:) = range(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_range(n,:) = range(DLS_L_Acc_norm);

            % RMS
            SC_Acc_x_rms(n,:) = rms(SC_Acc(:,1));
            SC_Acc_y_rms(n,:) = rms(SC_Acc(:,2));
            SC_Acc_z_rms(n,:) = rms(SC_Acc(:,3));
            SC_Acc_norm_rms(n,:) = rms(SC_Acc_norm);
            DLS_R_Acc_x_rms(n,:) = rms(DLS_R_Acc(:,1));
            DLS_R_Acc_y_rms(n,:) = rms(DLS_R_Acc(:,2));
            DLS_R_Acc_z_rms(n,:) = rms(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_rms(n,:) = rms(DLS_R_Acc_norm);
            DLS_L_Acc_x_rms(n,:) = rms(DLS_L_Acc(:,1));
            DLS_L_Acc_y_rms(n,:) = rms(DLS_L_Acc(:,2));
            DLS_L_Acc_z_rms(n,:) = rms(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_rms(n,:) = rms(DLS_L_Acc_norm);

            % Standard Deviation
            SC_Acc_x_std(n,:) = std(SC_Acc(:,1));
            SC_Acc_y_std(n,:) = std(SC_Acc(:,2));
            SC_Acc_z_std(n,:) = std(SC_Acc(:,3));
            SC_Acc_norm_std(n,:) = std(SC_Acc_norm);
            DLS_R_Acc_x_std(n,:) = std(DLS_R_Acc(:,1));
            DLS_R_Acc_y_std(n,:) = std(DLS_R_Acc(:,2));
            DLS_R_Acc_z_std(n,:) = std(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_std(n,:) = std(DLS_R_Acc_norm);
            DLS_L_Acc_x_std(n,:) = std(DLS_L_Acc(:,1));
            DLS_L_Acc_y_std(n,:) = std(DLS_L_Acc(:,2));
            DLS_L_Acc_z_std(n,:) = std(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_std(n,:) = std(DLS_L_Acc_norm);

            % Skew
            SC_Acc_x_skew(n,:) = skewness(SC_Acc(:,1));
            SC_Acc_y_skew(n,:) = skewness(SC_Acc(:,2));
            SC_Acc_z_skew(n,:) = skewness(SC_Acc(:,3));
            SC_Acc_norm_skew(n,:) = skewness(SC_Acc_norm);
            DLS_R_Acc_x_skew(n,:) = skewness(DLS_R_Acc(:,1));
            DLS_R_Acc_y_skew(n,:) = skewness(DLS_R_Acc(:,2));
            DLS_R_Acc_z_skew(n,:) = skewness(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_skew(n,:) = skewness(DLS_R_Acc_norm);
            DLS_L_Acc_x_skew(n,:) = skewness(DLS_L_Acc(:,1));
            DLS_L_Acc_y_skew(n,:) = skewness(DLS_L_Acc(:,2));
            DLS_L_Acc_z_skew(n,:) = skewness(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_skew(n,:) = skewness(DLS_L_Acc_norm);

            % Kurtosis
            SC_Acc_x_kurtosis(n,:) = kurtosis(SC_Acc(:,1));
            SC_Acc_y_kurtosis(n,:) = kurtosis(SC_Acc(:,2));
            SC_Acc_z_kurtosis(n,:) = kurtosis(SC_Acc(:,3));
            SC_Acc_norm_kurtosis(n,:) = kurtosis(SC_Acc_norm);
            DLS_R_Acc_x_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,1));
            DLS_R_Acc_y_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,2));
            DLS_R_Acc_z_kurtosis(n,:) = kurtosis(DLS_R_Acc(:,3));
            DLS_R_Acc_norm_kurtosis(n,:) = kurtosis(DLS_R_Acc_norm);
            DLS_L_Acc_x_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,1));
            DLS_L_Acc_y_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,2));
            DLS_L_Acc_z_kurtosis(n,:) = kurtosis(DLS_L_Acc(:,3));
            DLS_L_Acc_norm_kurtosis(n,:) = kurtosis(DLS_L_Acc_norm);

            % Pearson correlation coefficient
            SC_Acc_corr_xy(n,:) = corr(SC_Acc(:,1),SC_Acc(:,2));
            SC_Acc_corr_xz(n,:) = corr(SC_Acc(:,1),SC_Acc(:,3));
            SC_Acc_corr_yz(n,:) = corr(SC_Acc(:,2),SC_Acc(:,3));
            DLS_R_Acc_corr_xy(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,2));
            DLS_R_Acc_corr_xz(n,:) = corr(DLS_R_Acc(:,1),DLS_R_Acc(:,3));
            DLS_R_Acc_corr_yz(n,:) = corr(DLS_R_Acc(:,2),DLS_R_Acc(:,3));
            DLS_L_Acc_corr_xy(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,2));
            DLS_L_Acc_corr_xz(n,:) = corr(DLS_L_Acc(:,1),DLS_L_Acc(:,3));
            DLS_L_Acc_corr_yz(n,:) = corr(DLS_L_Acc(:,2),DLS_L_Acc(:,3));

            % Sample Entropy
            r = 0.2;
            SC_Acc_x_SamEn(n,:) = sampen(SC_Acc(:,1),1,r);
            SC_Acc_y_SamEn(n,:) = sampen(SC_Acc(:,2),1,r);
            SC_Acc_z_SamEn(n,:) = sampen(SC_Acc(:,3),1,r);
            SC_Acc_norm_SamEn(n,:) = sampen(SC_Acc_norm,1,r);
            DLS_R_Acc_x_SamEn(n,:) = sampen(DLS_R_Acc(:,1),1,r);
            DLS_R_Acc_y_SamEn(n,:) = sampen(DLS_R_Acc(:,2),1,r);
            DLS_R_Acc_z_SamEn(n,:) = sampen(DLS_R_Acc(:,3),1,r);
            DLS_R_Acc_norm_SamEn(n,:) = sampen(DLS_R_Acc_norm,1,r);
            DLS_L_Acc_x_SamEn(n,:) = sampen(DLS_L_Acc(:,1),1,r);
            DLS_L_Acc_y_SamEn(n,:) = sampen(DLS_L_Acc(:,2),1,r);
            DLS_L_Acc_z_SamEn(n,:) = sampen(DLS_L_Acc(:,3),1,r);
            DLS_L_Acc_norm_SamEn(n,:) = sampen(DLS_L_Acc_norm,1,r);
            
            % Frequency Domain
            ff = FFeatures(SC_Acc(:,1), Hz);
            SC_Acc_x_DAmp(n,:) = ff(1);
            SC_Acc_x_DFreq(n,:) = ff(2);
            SC_Acc_x_PSD_mean(n,:) = ff(3);
            SC_Acc_x_PSD_std(n,:) = ff(4);
            SC_Acc_x_PSD_skew(n,:) = ff(5);
            SC_Acc_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(SC_Acc(:,2), Hz);
            SC_Acc_y_DAmp(n,:) = ff(1);
            SC_Acc_y_DFreq(n,:) = ff(2);
            SC_Acc_y_PSD_mean(n,:) = ff(3);
            SC_Acc_y_PSD_std(n,:) = ff(4);
            SC_Acc_y_PSD_skew(n,:) = ff(5);
            SC_Acc_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(SC_Acc(:,3), Hz);
            SC_Acc_z_DAmp(n,:) = ff(1);
            SC_Acc_z_DFreq(n,:) = ff(2);
            SC_Acc_z_PSD_mean(n,:) = ff(3);
            SC_Acc_z_PSD_std(n,:) = ff(4);
            SC_Acc_z_PSD_skew(n,:) = ff(5);
            SC_Acc_z_PSD_kurtosis(n,:) =ff(6);

            ff = FFeatures(SC_Acc_norm, Hz);
            SC_Acc_norm_DAmp(n,:) = ff(1);
            SC_Acc_norm_DFreq(n,:) = ff(2);
            SC_Acc_norm_PSD_mean(n,:) = ff(3);
            SC_Acc_norm_PSD_std(n,:) = ff(4);
            SC_Acc_norm_PSD_skew(n,:) = ff(5);
            SC_Acc_norm_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Acc(:,1), Hz);
            DLS_R_Acc_x_DAmp(n,:) = ff(1);
            DLS_R_Acc_x_DFreq(n,:) = ff(2);
            DLS_R_Acc_x_PSD_mean(n,:) = ff(3);
            DLS_R_Acc_x_PSD_std(n,:) = ff(4);
            DLS_R_Acc_x_PSD_skew(n,:) = ff(5);
            DLS_R_Acc_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Acc(:,2), Hz);
            DLS_R_Acc_y_DAmp(n,:) = ff(1);
            DLS_R_Acc_y_DFreq(n,:) = ff(2);
            DLS_R_Acc_y_PSD_mean(n,:) = ff(3);
            DLS_R_Acc_y_PSD_std(n,:) = ff(4);
            DLS_R_Acc_y_PSD_skew(n,:) = ff(5);
            DLS_R_Acc_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Acc(:,3), Hz);
            DLS_R_Acc_z_DAmp(n,:) = ff(1);
            DLS_R_Acc_z_DFreq(n,:) = ff(2);
            DLS_R_Acc_z_PSD_mean(n,:) = ff(3);
            DLS_R_Acc_z_PSD_std(n,:) = ff(4);
            DLS_R_Acc_z_PSD_skew(n,:) = ff(5);
            DLS_R_Acc_z_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_R_Acc_norm, Hz);
            DLS_R_Acc_norm_DAmp(n,:) = ff(1);
            DLS_R_Acc_norm_DFreq(n,:) = ff(2);
            DLS_R_Acc_norm_PSD_mean(n,:) = ff(3);
            DLS_R_Acc_norm_PSD_std(n,:) = ff(4);
            DLS_R_Acc_norm_PSD_skew(n,:) = ff(5);
            DLS_R_Acc_norm_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Acc(:,1), Hz);
            DLS_L_Acc_x_DAmp(n,:) = ff(1);
            DLS_L_Acc_x_DFreq(n,:) = ff(2);
            DLS_L_Acc_x_PSD_mean(n,:) = ff(3);
            DLS_L_Acc_x_PSD_std(n,:) = ff(4);
            DLS_L_Acc_x_PSD_skew(n,:) = ff(5);
            DLS_L_Acc_x_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Acc(:,2), Hz);
            DLS_L_Acc_y_DAmp(n,:) = ff(1);
            DLS_L_Acc_y_DFreq(n,:) = ff(2);
            DLS_L_Acc_y_PSD_mean(n,:) = ff(3);
            DLS_L_Acc_y_PSD_std(n,:) = ff(4);
            DLS_L_Acc_y_PSD_skew(n,:) = ff(5);
            DLS_L_Acc_y_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Acc(:,3), Hz);
            DLS_L_Acc_z_DAmp(n,:) = ff(1);
            DLS_L_Acc_z_DFreq(n,:) = ff(2);
            DLS_L_Acc_z_PSD_mean(n,:) = ff(3);
            DLS_L_Acc_z_PSD_std(n,:) = ff(4);
            DLS_L_Acc_z_PSD_skew(n,:) = ff(5);
            DLS_L_Acc_z_PSD_kurtosis(n,:) = ff(6);

            ff = FFeatures(DLS_L_Acc_norm, Hz);
            DLS_L_Acc_norm_DAmp(n,:) = ff(1);
            DLS_L_Acc_norm_DFreq(n,:) = ff(2);
            DLS_L_Acc_norm_PSD_mean(n,:) = ff(3);
            DLS_L_Acc_norm_PSD_std(n,:) = ff(4);
            DLS_L_Acc_norm_PSD_skew(n,:) = ff(5);
            DLS_L_Acc_norm_PSD_kurtosis(n,:) = ff(6);


        end
        subject{n,:} = ['CVA' sprintf('%02d',n)]
        activity{n,:} = ['6MWT']
        cut_off_time(n,:) = TS(k)
    end

    CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_Admission.csv');
    FIM_AD = CO(:,2);
    BBS_AD = CO(:,3);
    MWT10_SSV_AD = CO(:,4);
    MWT10_FV_AD = CO(:,5);
    MWT6_AD = CO(:,6);
    TUG_AD = CO(:,7);
    
    for  i = 1:1:length(MWT10_SSV_AD)
        if MWT10_SSV_AD(i) < 0.4
            Ambul_AD{i,:} = 'Severe';
        elseif MWT10_SSV_AD(i) >= 0.4 && MWT10_SSV_AD(i) < 0.8
            Ambul_AD{i,:} = 'Moderate';
        elseif MWT10_SSV_AD(i) >= 0.8
            Ambul_AD{i,:} = 'Mild';
        end
    end

    CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_Discharge.csv');
    FIM_DC = CO(:,2);
    BBS_DC = CO(:,3);
    MWT10_SSV_DC = CO(:,4);
    MWT10_FV_DC = CO(:,5);
    MWT6_DC = CO(:,6);
    TUG_DC = CO(:,7);

    for  i = 1:1:length(MWT10_SSV_DC)
        if MWT10_SSV_DC(i) < 0.4
            Ambul_DC{i,:} = 'Severe';
        elseif MWT10_SSV_DC(i) >= 0.4 && MWT10_SSV_DC(i) < 0.8
            Ambul_DC{i,:} = 'Moderate';
        elseif MWT10_SSV_DC(i) >= 0.8                           
            Ambul_DC{i,:} = 'Mild';
        end
    end
    

    


    tbl_Features = table(subject, activity, cut_off_time, FIM_AD, BBS_AD, MWT10_SSV_AD, MWT10_FV_AD, MWT6_AD, TUG_AD, ...
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
        DLS_L_Acc_norm_DAmp, DLS_L_Acc_norm_DFreq, DLS_L_Acc_norm_PSD_mean, DLS_L_Acc_norm_PSD_std, DLS_L_Acc_norm_PSD_skew, DLS_L_Acc_norm_PSD_kurtosis, ...
        Ambul_AD, Ambul_DC)
    
    
    tbl = [tbl; tbl_Features]
end


writetable(tbl,'Feature_Matrix_Prediction_6MWT.csv','Delimiter',',','QuoteStrings',true)









