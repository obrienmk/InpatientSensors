clc
clear all
close all

% CVA
% ID: 1, 3, 22, 25, 30, 32, 42, 49  --> no data
% ID: 7, 40 --> sensor at DLS_L missing
% ID: 8 --> sensor at DLS_R missing

ID = [1:55];

Type_of_Subject = 'CVA'

Hz = 31.25;
dt = 1/Hz; 
k=1;
for n = 1:1:length(ID)
    file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID(n)) '.mat']
    load(file_input);
    
    Side = get_side('./CVA_Paretic_Side.csv');
    AS = Side.Side(ID(n));

    if n == 1 || n == 3 || n == 7 || n == 8 || n == 22 || n == 25 || n == 30 || n == 32 || n == 40 || n == 42 || n == 49
         SN = nan;     % Session #
         TN = nan;     % Trial #
         
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
        
        
        % Derivative
        % Gyro mean
        dSC_Gyr_x_mean(n,:) = nan;
        dSC_Gyr_y_mean(n,:) = nan;
        dSC_Gyr_z_mean(n,:) = nan;
        dSC_Gyr_norm_mean(n,:) = nan;
        dDLS_R_Gyr_x_mean(n,:) = nan;
        dDLS_R_Gyr_y_mean(n,:) = nan;
        dDLS_R_Gyr_z_mean(n,:) = nan;
        dDLS_R_Gyr_norm_mean(n,:) = nan;
        dDLS_L_Gyr_x_mean(n,:) = nan;
        dDLS_L_Gyr_y_mean(n,:) = nan;
        dDLS_L_Gyr_z_mean(n,:) = nan;
        dDLS_L_Gyr_norm_mean(n,:) = nan;

        % Range
        dSC_Gyr_x_range(n,:) = nan;
        dSC_Gyr_y_range(n,:) = nan;
        dSC_Gyr_z_range(n,:) = nan;
        dSC_Gyr_norm_range(n,:) = nan;
        dDLS_R_Gyr_x_range(n,:) = nan;
        dDLS_R_Gyr_y_range(n,:) = nan;
        dDLS_R_Gyr_z_range(n,:) = nan;
        dDLS_R_Gyr_norm_range(n,:) = nan;
        dDLS_L_Gyr_x_range(n,:) = nan;
        dDLS_L_Gyr_y_range(n,:) = nan;
        dDLS_L_Gyr_z_range(n,:) = nan;
        dDLS_L_Gyr_norm_range(n,:) = nan;

        % RMS
        dSC_Gyr_x_rms(n,:) = nan;
        dSC_Gyr_y_rms(n,:) = nan;
        dSC_Gyr_z_rms(n,:) = nan;
        dSC_Gyr_norm_rms(n,:) = nan;
        dDLS_R_Gyr_x_rms(n,:) = nan;
        dDLS_R_Gyr_y_rms(n,:) = nan;
        dDLS_R_Gyr_z_rms(n,:) = nan;
        dDLS_R_Gyr_norm_rms(n,:) = nan;
        dDLS_L_Gyr_x_rms(n,:) = nan;
        dDLS_L_Gyr_y_rms(n,:) = nan;
        dDLS_L_Gyr_z_rms(n,:) = nan;
        dDLS_L_Gyr_norm_rms(n,:) = nan;

        % Standard Deviation
        dSC_Gyr_x_std(n,:) = nan;
        dSC_Gyr_y_std(n,:) = nan;
        dSC_Gyr_z_std(n,:) = nan;
        dSC_Gyr_norm_std(n,:) = nan;
        dDLS_R_Gyr_x_std(n,:) = nan;
        dDLS_R_Gyr_y_std(n,:) = nan;
        dDLS_R_Gyr_z_std(n,:) = nan;
        dDLS_R_Gyr_norm_std(n,:) = nan;
        dDLS_L_Gyr_x_std(n,:) = nan;
        dDLS_L_Gyr_y_std(n,:) = nan;
        dDLS_L_Gyr_z_std(n,:) = nan;
        dDLS_L_Gyr_norm_std(n,:) = nan;

        % Skew
        dSC_Gyr_x_skew(n,:) = nan;
        dSC_Gyr_y_skew(n,:) = nan;
        dSC_Gyr_z_skew(n,:) = nan;
        dSC_Gyr_norm_skew(n,:) = nan;
        dDLS_R_Gyr_x_skew(n,:) = nan;
        dDLS_R_Gyr_y_skew(n,:) = nan;
        dDLS_R_Gyr_z_skew(n,:) = nan;
        dDLS_R_Gyr_norm_skew(n,:) = nan;
        dDLS_L_Gyr_x_skew(n,:) = nan;
        dDLS_L_Gyr_y_skew(n,:) = nan;
        dDLS_L_Gyr_z_skew(n,:) = nan;
        dDLS_L_Gyr_norm_skew(n,:) = nan;

        % Kurtosis
        dSC_Gyr_x_kurtosis(n,:) = nan;
        dSC_Gyr_y_kurtosis(n,:) = nan;
        dSC_Gyr_z_kurtosis(n,:) = nan;
        dSC_Gyr_norm_kurtosis(n,:) = nan;
        dDLS_R_Gyr_x_kurtosis(n,:) = nan;
        dDLS_R_Gyr_y_kurtosis(n,:) = nan;
        dDLS_R_Gyr_z_kurtosis(n,:) = nan;
        dDLS_R_Gyr_norm_kurtosis(n,:) = nan;
        dDLS_L_Gyr_x_kurtosis(n,:) = nan;
        dDLS_L_Gyr_y_kurtosis(n,:) = nan;
        dDLS_L_Gyr_z_kurtosis(n,:) = nan;
        dDLS_L_Gyr_norm_kurtosis(n,:) = nan;
        
        
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

       
    else

         data
         SN = 1;     % Session #: Admission
         TN = 1;     % Trial #

        % Gyroscope data
        SC_Gyr = data.Session{SN}.Motion.SC.Gyr{TN};
        DLS_R_Gyr = data.Session{SN}.Motion.DLS_R.Gyr{TN};
        DLS_L_Gyr = data.Session{SN}.Motion.DLS_L.Gyr{TN};

        % Acceleration data
        SC_Acc = data.Session{SN}.Motion.SC.Acc{TN};
        DLS_R_Acc = data.Session{SN}.Motion.DLS_R.Acc{TN};
        DLS_L_Acc = data.Session{SN}.Motion.DLS_L.Acc{TN};

        % Demean Acc data
        SC_Acc = SC_Acc - ones(length(SC_Acc),1)*mean(SC_Acc);
        DLS_R_Acc = DLS_R_Acc - ones(length(DLS_R_Acc),1)*mean(DLS_R_Acc);
        DLS_L_Acc = DLS_L_Acc - ones(length(DLS_L_Acc),1)*mean(DLS_L_Acc);

        % To match final time
        final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
%         final = Hz*20;  % Initial 20 sec
        Time = data.Session{SN}.Motion.Time{TN}(1:final,:);
        SC_Gyr = SC_Gyr(1:final,:);
        DLS_R_Gyr = DLS_R_Gyr(1:final,:);
        DLS_L_Gyr = DLS_L_Gyr(1:final,:);
        SC_Acc = SC_Acc(1:final,:);
        DLS_R_Acc = DLS_R_Acc(1:final,:);
        DLS_L_Acc = DLS_L_Acc(1:final,:);
        
        % Norm of Gyro and Acc Data
        for i = 1:1:length(Time)
            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
            DLS_R_Gyr_norm(i,:) = norm(DLS_R_Gyr(i,:));
            DLS_L_Gyr_norm(i,:) = norm(DLS_L_Gyr(i,:));

            SC_Acc_norm(i,:) = norm(SC_Acc(i,:));
            DLS_R_Acc_norm(i,:) = norm(DLS_R_Acc(i,:));
            DLS_L_Acc_norm(i,:) = norm(DLS_L_Acc(i,:));
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
        
        
        % Derivative
        % Gyro mean
        dSC_Gyr_x_mean(n,:) = mean(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_mean(n,:) = mean(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_mean(n,:) = mean(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_mean(n,:) = mean(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_mean(n,:) = mean(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_mean(n,:) = mean(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_mean(n,:) = mean(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_mean(n,:) = mean(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_mean(n,:) = mean(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_mean(n,:) = mean(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_mean(n,:) = mean(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_mean(n,:) = mean(diff(DLS_L_Gyr_norm)/dt);

        % Range
        dSC_Gyr_x_range(n,:) = range(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_range(n,:) = range(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_range(n,:) = range(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_range(n,:) = range(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_range(n,:) = range(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_range(n,:) = range(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_range(n,:) = range(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_range(n,:) = range(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_range(n,:) = range(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_range(n,:) = range(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_range(n,:) = range(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_range(n,:) = range(diff(DLS_L_Gyr_norm)/dt);

        % RMS
        dSC_Gyr_x_rms(n,:) = rms(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_rms(n,:) = rms(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_rms(n,:) = rms(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_rms(n,:) = rms(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_rms(n,:) = rms(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_rms(n,:) = rms(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_rms(n,:) = rms(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_rms(n,:) = rms(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_rms(n,:) = rms(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_rms(n,:) = rms(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_rms(n,:) = rms(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_rms(n,:) = rms(diff(DLS_L_Gyr_norm)/dt);

        % Standard Deviation
        dSC_Gyr_x_std(n,:) = std(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_std(n,:) = std(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_std(n,:) = std(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_std(n,:) = std(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_std(n,:) = std(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_std(n,:) = std(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_std(n,:) = std(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_std(n,:) = std(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_std(n,:) = std(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_std(n,:) = std(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_std(n,:) = std(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_std(n,:) = std(diff(DLS_L_Gyr_norm)/dt);

        % Skew
        dSC_Gyr_x_skew(n,:) = skewness(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_skew(n,:) = skewness(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_skew(n,:) = skewness(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_skew(n,:) = skewness(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_skew(n,:) = skewness(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_skew(n,:) = skewness(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_skew(n,:) = skewness(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_skew(n,:) = skewness(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_skew(n,:) = skewness(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_skew(n,:) = skewness(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_skew(n,:) = skewness(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_skew(n,:) = skewness(diff(DLS_L_Gyr_norm)/dt);

        % Kurtosis
        dSC_Gyr_x_kurtosis(n,:) = kurtosis(diff(SC_Gyr(:,1))/dt);
        dSC_Gyr_y_kurtosis(n,:) = kurtosis(diff(SC_Gyr(:,2))/dt);
        dSC_Gyr_z_kurtosis(n,:) = kurtosis(diff(SC_Gyr(:,3))/dt);
        dSC_Gyr_norm_kurtosis(n,:) = kurtosis(diff(SC_Gyr_norm)/dt);
        dDLS_R_Gyr_x_kurtosis(n,:) = kurtosis(diff(DLS_R_Gyr(:,1))/dt);
        dDLS_R_Gyr_y_kurtosis(n,:) = kurtosis(diff(DLS_R_Gyr(:,2))/dt);
        dDLS_R_Gyr_z_kurtosis(n,:) = kurtosis(diff(DLS_R_Gyr(:,3))/dt);
        dDLS_R_Gyr_norm_kurtosis(n,:) = kurtosis(diff(DLS_R_Gyr_norm)/dt);
        dDLS_L_Gyr_x_kurtosis(n,:) = kurtosis(diff(DLS_L_Gyr(:,1))/dt);
        dDLS_L_Gyr_y_kurtosis(n,:) = kurtosis(diff(DLS_L_Gyr(:,2))/dt);
        dDLS_L_Gyr_z_kurtosis(n,:) = kurtosis(diff(DLS_L_Gyr(:,3))/dt);
        dDLS_L_Gyr_norm_kurtosis(n,:) = kurtosis(diff(DLS_L_Gyr_norm)/dt);
        
        
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
        SC_Gyr_z_PSD_kurtosis(n,:) = ff(6);

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
        

    end
    
end
ID = ID';


tbl_features = table(ID, SC_Gyr_x_mean, SC_Gyr_y_mean, SC_Gyr_z_mean, SC_Gyr_norm_mean, ... 
                    DLS_L_Gyr_x_mean, DLS_L_Gyr_y_mean, DLS_L_Gyr_z_mean, DLS_L_Gyr_norm_mean, ...
                    DLS_R_Gyr_x_mean, DLS_R_Gyr_y_mean, DLS_R_Gyr_z_mean, DLS_R_Gyr_norm_mean, ...
                    SC_Gyr_x_range, SC_Gyr_y_range, SC_Gyr_z_range, SC_Gyr_norm_range, ... 
                    DLS_L_Gyr_x_range, DLS_L_Gyr_y_range, DLS_L_Gyr_z_range, DLS_L_Gyr_norm_range, ...
                    DLS_R_Gyr_x_range, DLS_R_Gyr_y_range, DLS_R_Gyr_z_range, DLS_R_Gyr_norm_range, ...
                    SC_Gyr_x_rms, SC_Gyr_y_rms, SC_Gyr_z_rms, SC_Gyr_norm_rms, ... 
                    DLS_L_Gyr_x_rms, DLS_L_Gyr_y_rms, DLS_L_Gyr_z_rms, DLS_L_Gyr_norm_rms, ...
                    DLS_R_Gyr_x_rms, DLS_R_Gyr_y_rms, DLS_R_Gyr_z_rms, DLS_R_Gyr_norm_rms, ...
                    SC_Gyr_x_std, SC_Gyr_y_std, SC_Gyr_z_std, SC_Gyr_norm_std, ... 
                    DLS_L_Gyr_x_std, DLS_L_Gyr_y_std, DLS_L_Gyr_z_std, DLS_L_Gyr_norm_std, ...
                    DLS_R_Gyr_x_std, DLS_R_Gyr_y_std, DLS_R_Gyr_z_std, DLS_R_Gyr_norm_std, ...
                    SC_Gyr_x_skew, SC_Gyr_y_skew, SC_Gyr_z_skew, SC_Gyr_norm_skew, ... 
                    DLS_L_Gyr_x_skew, DLS_L_Gyr_y_skew, DLS_L_Gyr_z_skew, DLS_L_Gyr_norm_skew, ...
                    DLS_R_Gyr_x_skew, DLS_R_Gyr_y_skew, DLS_R_Gyr_z_skew, DLS_R_Gyr_norm_skew, ...
                    SC_Gyr_x_kurtosis, SC_Gyr_y_kurtosis, SC_Gyr_z_kurtosis, SC_Gyr_norm_kurtosis, ... 
                    DLS_L_Gyr_x_kurtosis, DLS_L_Gyr_y_kurtosis, DLS_L_Gyr_z_kurtosis, DLS_L_Gyr_norm_kurtosis, ...
                    DLS_R_Gyr_x_kurtosis, DLS_R_Gyr_y_kurtosis, DLS_R_Gyr_z_kurtosis, DLS_R_Gyr_norm_kurtosis, ...
                    dSC_Gyr_x_mean, dSC_Gyr_y_mean, dSC_Gyr_z_mean, dSC_Gyr_norm_mean, ... 
                    dDLS_L_Gyr_x_mean, dDLS_L_Gyr_y_mean, dDLS_L_Gyr_z_mean, dDLS_L_Gyr_norm_mean, ...
                    dDLS_R_Gyr_x_mean, dDLS_R_Gyr_y_mean, dDLS_R_Gyr_z_mean, dDLS_R_Gyr_norm_mean, ...
                    dSC_Gyr_x_range, dSC_Gyr_y_range, dSC_Gyr_z_range, dSC_Gyr_norm_range, ... 
                    dDLS_L_Gyr_x_range, dDLS_L_Gyr_y_range, dDLS_L_Gyr_z_range, dDLS_L_Gyr_norm_range, ...
                    dDLS_R_Gyr_x_range, dDLS_R_Gyr_y_range, dDLS_R_Gyr_z_range, dDLS_R_Gyr_norm_range, ...
                    dSC_Gyr_x_rms, dSC_Gyr_y_rms, dSC_Gyr_z_rms, dSC_Gyr_norm_rms, ... 
                    dDLS_L_Gyr_x_rms, dDLS_L_Gyr_y_rms, dDLS_L_Gyr_z_rms, dDLS_L_Gyr_norm_rms, ...
                    dDLS_R_Gyr_x_rms, dDLS_R_Gyr_y_rms, dDLS_R_Gyr_z_rms, dDLS_R_Gyr_norm_rms, ...
                    dSC_Gyr_x_std, dSC_Gyr_y_std, dSC_Gyr_z_std, dSC_Gyr_norm_std, ... 
                    dDLS_L_Gyr_x_std, dDLS_L_Gyr_y_std, dDLS_L_Gyr_z_std, dDLS_L_Gyr_norm_std, ...
                    dDLS_R_Gyr_x_std, dDLS_R_Gyr_y_std, dDLS_R_Gyr_z_std, dDLS_R_Gyr_norm_std, ...
                    dSC_Gyr_x_skew, dSC_Gyr_y_skew, dSC_Gyr_z_skew, dSC_Gyr_norm_skew, ... 
                    dDLS_L_Gyr_x_skew, dDLS_L_Gyr_y_skew, dDLS_L_Gyr_z_skew, dDLS_L_Gyr_norm_skew, ...
                    dDLS_R_Gyr_x_skew, dDLS_R_Gyr_y_skew, dDLS_R_Gyr_z_skew, dDLS_R_Gyr_norm_skew, ...
                    dSC_Gyr_x_kurtosis, dSC_Gyr_y_kurtosis, dSC_Gyr_z_kurtosis, dSC_Gyr_norm_kurtosis, ... 
                    dDLS_L_Gyr_x_kurtosis, dDLS_L_Gyr_y_kurtosis, dDLS_L_Gyr_z_kurtosis, dDLS_L_Gyr_norm_kurtosis, ...
                    dDLS_R_Gyr_x_kurtosis, dDLS_R_Gyr_y_kurtosis, dDLS_R_Gyr_z_kurtosis, dDLS_R_Gyr_norm_kurtosis, ...
                    SC_Gyr_corr_xy, SC_Gyr_corr_xz, SC_Gyr_corr_yz, DLS_R_Gyr_corr_xy, DLS_R_Gyr_corr_xz, DLS_R_Gyr_corr_yz, DLS_L_Gyr_corr_xy, DLS_L_Gyr_corr_xz, DLS_L_Gyr_corr_yz, ...
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
                    DLS_L_Gyr_norm_DAmp, DLS_L_Gyr_norm_DFreq, DLS_L_Gyr_norm_PSD_mean, DLS_L_Gyr_norm_PSD_std, DLS_L_Gyr_norm_PSD_skew, DLS_L_Gyr_norm_PSD_kurtosis ...
                )   

writetable(tbl_features,'General_Feature_Matrix_Admission_6MWT.csv','Delimiter',',','QuoteStrings',true)



