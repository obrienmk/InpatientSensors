clc
clear all
close all

% CVA
% ID: 1, 3, 22, 25, 30, 32, 42, 44, 49  --> no data
% ID: 7, 40 --> sensor at DLS_L missing
% ID: 8 --> sensor at DLS_R missing


Type_of_Subject = 'CVA'

Hz = 31.25;
dt = 1/Hz; 

ID = 2
file_input = ['\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Data analysis\2_Clean_Data_Extracted\6MWT\' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID) '.mat']
load(file_input);
    

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
SC_Gyr_x_mean = mean(SC_Gyr(:,1));
SC_Gyr_y_mean = mean(SC_Gyr(:,2));
SC_Gyr_z_mean = mean(SC_Gyr(:,3));
SC_Gyr_norm_mean = mean(SC_Gyr_norm);
DLS_R_Gyr_x_mean = mean(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_mean = mean(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_mean = mean(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_mean = mean(DLS_R_Gyr_norm);
DLS_L_Gyr_x_mean = mean(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_mean = mean(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_mean = mean(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_mean = mean(DLS_L_Gyr_norm);

% Range
SC_Gyr_x_range = range(SC_Gyr(:,1));
SC_Gyr_y_range = range(SC_Gyr(:,2));
SC_Gyr_z_range = range(SC_Gyr(:,3));
SC_Gyr_norm_range = range(SC_Gyr_norm);
DLS_R_Gyr_x_range = range(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_range = range(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_range = range(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_range = range(DLS_R_Gyr_norm);
DLS_L_Gyr_x_range = range(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_range = range(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_range = range(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_range = range(DLS_L_Gyr_norm);

% RMS
SC_Gyr_x_rms = rms(SC_Gyr(:,1));
SC_Gyr_y_rms = rms(SC_Gyr(:,2));
SC_Gyr_z_rms = rms(SC_Gyr(:,3));
SC_Gyr_norm_rms = rms(SC_Gyr_norm);
DLS_R_Gyr_x_rms = rms(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_rms = rms(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_rms = rms(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_rms = rms(DLS_R_Gyr_norm);
DLS_L_Gyr_x_rms = rms(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_rms = rms(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_rms = rms(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_rms = rms(DLS_L_Gyr_norm);

% Standard Deviation
SC_Gyr_x_std = std(SC_Gyr(:,1));
SC_Gyr_y_std = std(SC_Gyr(:,2));
SC_Gyr_z_std = std(SC_Gyr(:,3));
SC_Gyr_norm_std = std(SC_Gyr_norm);
DLS_R_Gyr_x_std = std(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_std = std(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_std = std(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_std = std(DLS_R_Gyr_norm);
DLS_L_Gyr_x_std = std(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_std = std(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_std = std(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_std = std(DLS_L_Gyr_norm);

% Skew
SC_Gyr_x_skew = skewness(SC_Gyr(:,1));
SC_Gyr_y_skew = skewness(SC_Gyr(:,2));
SC_Gyr_z_skew = skewness(SC_Gyr(:,3));
SC_Gyr_norm_skew = skewness(SC_Gyr_norm);
DLS_R_Gyr_x_skew = skewness(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_skew = skewness(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_skew = skewness(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_skew = skewness(DLS_R_Gyr_norm);
DLS_L_Gyr_x_skew = skewness(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_skew = skewness(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_skew = skewness(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_skew = skewness(DLS_L_Gyr_norm);

% Kurtosis
SC_Gyr_x_kurtosis = kurtosis(SC_Gyr(:,1));
SC_Gyr_y_kurtosis = kurtosis(SC_Gyr(:,2));
SC_Gyr_z_kurtosis = kurtosis(SC_Gyr(:,3));
SC_Gyr_norm_kurtosis = kurtosis(SC_Gyr_norm);
DLS_R_Gyr_x_kurtosis = kurtosis(DLS_R_Gyr(:,1));
DLS_R_Gyr_y_kurtosis = kurtosis(DLS_R_Gyr(:,2));
DLS_R_Gyr_z_kurtosis = kurtosis(DLS_R_Gyr(:,3));
DLS_R_Gyr_norm_kurtosis = kurtosis(DLS_R_Gyr_norm);
DLS_L_Gyr_x_kurtosis = kurtosis(DLS_L_Gyr(:,1));
DLS_L_Gyr_y_kurtosis = kurtosis(DLS_L_Gyr(:,2));
DLS_L_Gyr_z_kurtosis = kurtosis(DLS_L_Gyr(:,3));
DLS_L_Gyr_norm_kurtosis = kurtosis(DLS_L_Gyr_norm);
        

% Derivative
% Gyro mean
dSC_Gyr_x_mean = mean(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_mean = mean(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_mean = mean(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_mean = mean(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_mean = mean(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_mean = mean(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_mean = mean(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_mean = mean(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_mean = mean(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_mean = mean(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_mean = mean(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_mean = mean(diff(DLS_L_Gyr_norm)/dt);

% Range
dSC_Gyr_x_range = range(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_range = range(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_range = range(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_range = range(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_range = range(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_range = range(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_range = range(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_range = range(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_range = range(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_range = range(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_range = range(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_range = range(diff(DLS_L_Gyr_norm)/dt);

% RMS
dSC_Gyr_x_rms = rms(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_rms = rms(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_rms = rms(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_rms = rms(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_rms = rms(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_rms = rms(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_rms = rms(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_rms = rms(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_rms = rms(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_rms = rms(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_rms = rms(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_rms = rms(diff(DLS_L_Gyr_norm)/dt);

% Standard Deviation
dSC_Gyr_x_std = std(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_std = std(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_std = std(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_std = std(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_std = std(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_std = std(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_std = std(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_std = std(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_std = std(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_std = std(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_std = std(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_std = std(diff(DLS_L_Gyr_norm)/dt);

% Skew
dSC_Gyr_x_skew = skewness(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_skew = skewness(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_skew = skewness(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_skew = skewness(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_skew = skewness(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_skew = skewness(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_skew = skewness(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_skew = skewness(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_skew = skewness(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_skew = skewness(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_skew = skewness(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_skew = skewness(diff(DLS_L_Gyr_norm)/dt);

% Kurtosis
dSC_Gyr_x_kurtosis = kurtosis(diff(SC_Gyr(:,1))/dt);
dSC_Gyr_y_kurtosis = kurtosis(diff(SC_Gyr(:,2))/dt);
dSC_Gyr_z_kurtosis = kurtosis(diff(SC_Gyr(:,3))/dt);
dSC_Gyr_norm_kurtosis = kurtosis(diff(SC_Gyr_norm)/dt);
dDLS_R_Gyr_x_kurtosis = kurtosis(diff(DLS_R_Gyr(:,1))/dt);
dDLS_R_Gyr_y_kurtosis = kurtosis(diff(DLS_R_Gyr(:,2))/dt);
dDLS_R_Gyr_z_kurtosis = kurtosis(diff(DLS_R_Gyr(:,3))/dt);
dDLS_R_Gyr_norm_kurtosis = kurtosis(diff(DLS_R_Gyr_norm)/dt);
dDLS_L_Gyr_x_kurtosis = kurtosis(diff(DLS_L_Gyr(:,1))/dt);
dDLS_L_Gyr_y_kurtosis = kurtosis(diff(DLS_L_Gyr(:,2))/dt);
dDLS_L_Gyr_z_kurtosis = kurtosis(diff(DLS_L_Gyr(:,3))/dt);
dDLS_L_Gyr_norm_kurtosis = kurtosis(diff(DLS_L_Gyr_norm)/dt);
        

% Pearson correlation coefficient
SC_Gyr_corr_xy = corr(SC_Gyr(:,1),SC_Gyr(:,2));
SC_Gyr_corr_xz = corr(SC_Gyr(:,1),SC_Gyr(:,3));
SC_Gyr_corr_yz = corr(SC_Gyr(:,2),SC_Gyr(:,3));
DLS_R_Gyr_corr_xy = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,2));
DLS_R_Gyr_corr_xz = corr(DLS_R_Gyr(:,1),DLS_R_Gyr(:,3));
DLS_R_Gyr_corr_yz = corr(DLS_R_Gyr(:,2),DLS_R_Gyr(:,3));
DLS_L_Gyr_corr_xy = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,2));
DLS_L_Gyr_corr_xz = corr(DLS_L_Gyr(:,1),DLS_L_Gyr(:,3));
DLS_L_Gyr_corr_yz = corr(DLS_L_Gyr(:,2),DLS_L_Gyr(:,3));


% Sample Entropy
r = 0.2
SC_Gyr_x_SamEn = sampen(SC_Gyr(:,1),1,r);
SC_Gyr_y_SamEn = sampen(SC_Gyr(:,2),1,r);
SC_Gyr_z_SamEn = sampen(SC_Gyr(:,3),1,r);
SC_Gyr_norm_SamEn = sampen(SC_Gyr_norm,1,r);
DLS_R_Gyr_x_SamEn = sampen(DLS_R_Gyr(:,1),1,r);
DLS_R_Gyr_y_SamEn = sampen(DLS_R_Gyr(:,2),1,r);
DLS_R_Gyr_z_SamEn = sampen(DLS_R_Gyr(:,3),1,r);
DLS_R_Gyr_norm_SamEn = sampen(DLS_R_Gyr_norm,1,r);
DLS_L_Gyr_x_SamEn = sampen(DLS_L_Gyr(:,1),1,r);
DLS_L_Gyr_y_SamEn = sampen(DLS_L_Gyr(:,2),1,r);
DLS_L_Gyr_z_SamEn = sampen(DLS_L_Gyr(:,3),1,r);
DLS_L_Gyr_norm_SamEn = sampen(DLS_L_Gyr_norm,1,r);


% Frequency Domain
ff = FFeatures(SC_Gyr(:,1), Hz);
SC_Gyr_x_DAmp = ff(1);
SC_Gyr_x_DFreq = ff(2);
SC_Gyr_x_PSD_mean = ff(3);
SC_Gyr_x_PSD_std = ff(4);
SC_Gyr_x_PSD_skew = ff(5);
SC_Gyr_x_PSD_kurtosis = ff(6);

ff = FFeatures(SC_Gyr(:,2), Hz);
SC_Gyr_y_DAmp = ff(1);
SC_Gyr_y_DFreq = ff(2);
SC_Gyr_y_PSD_mean = ff(3);
SC_Gyr_y_PSD_std = ff(4);
SC_Gyr_y_PSD_skew = ff(5);
SC_Gyr_y_PSD_kurtosis = ff(6);

ff = FFeatures(SC_Gyr(:,3), Hz);
SC_Gyr_z_DAmp = ff(1);
SC_Gyr_z_DFreq = ff(2);
SC_Gyr_z_PSD_mean = ff(3);
SC_Gyr_z_PSD_std = ff(4);
SC_Gyr_z_PSD_skew = ff(5);
SC_Gyr_z_PSD_kurtosis = ff(6);

ff = FFeatures(SC_Gyr_norm, Hz);
SC_Gyr_norm_DAmp = ff(1);
SC_Gyr_norm_DFreq = ff(2);
SC_Gyr_norm_PSD_mean = ff(3);
SC_Gyr_norm_PSD_std = ff(4);
SC_Gyr_norm_PSD_skew = ff(5);
SC_Gyr_norm_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_R_Gyr(:,1), Hz);
DLS_R_Gyr_x_DAmp = ff(1);
DLS_R_Gyr_x_DFreq = ff(2);
DLS_R_Gyr_x_PSD_mean = ff(3);
DLS_R_Gyr_x_PSD_std = ff(4);
DLS_R_Gyr_x_PSD_skew = ff(5);
DLS_R_Gyr_x_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_R_Gyr(:,2), Hz);
DLS_R_Gyr_y_DAmp = ff(1);
DLS_R_Gyr_y_DFreq = ff(2);
DLS_R_Gyr_y_PSD_mean = ff(3);
DLS_R_Gyr_y_PSD_std = ff(4);
DLS_R_Gyr_y_PSD_skew = ff(5);
DLS_R_Gyr_y_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_R_Gyr(:,3), Hz);
DLS_R_Gyr_z_DAmp = ff(1);
DLS_R_Gyr_z_DFreq = ff(2);
DLS_R_Gyr_z_PSD_mean = ff(3);
DLS_R_Gyr_z_PSD_std = ff(4);
DLS_R_Gyr_z_PSD_skew = ff(5);
DLS_R_Gyr_z_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_R_Gyr_norm, Hz);
DLS_R_Gyr_norm_DAmp = ff(1);
DLS_R_Gyr_norm_DFreq = ff(2);
DLS_R_Gyr_norm_PSD_mean = ff(3);
DLS_R_Gyr_norm_PSD_std = ff(4);
DLS_R_Gyr_norm_PSD_skew = ff(5);
DLS_R_Gyr_norm_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_L_Gyr(:,1), Hz);
DLS_L_Gyr_x_DAmp = ff(1);
DLS_L_Gyr_x_DFreq = ff(2);
DLS_L_Gyr_x_PSD_mean = ff(3);
DLS_L_Gyr_x_PSD_std = ff(4);
DLS_L_Gyr_x_PSD_skew = ff(5);
DLS_L_Gyr_x_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_L_Gyr(:,2), Hz);
DLS_L_Gyr_y_DAmp = ff(1);
DLS_L_Gyr_y_DFreq = ff(2);
DLS_L_Gyr_y_PSD_mean = ff(3);
DLS_L_Gyr_y_PSD_std = ff(4);
DLS_L_Gyr_y_PSD_skew = ff(5);
DLS_L_Gyr_y_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_L_Gyr(:,3), Hz);
DLS_L_Gyr_z_DAmp = ff(1);
DLS_L_Gyr_z_DFreq = ff(2);
DLS_L_Gyr_z_PSD_mean = ff(3);
DLS_L_Gyr_z_PSD_std = ff(4);
DLS_L_Gyr_z_PSD_skew = ff(5);
DLS_L_Gyr_z_PSD_kurtosis = ff(6);

ff = FFeatures(DLS_L_Gyr_norm, Hz);
DLS_L_Gyr_norm_DAmp = ff(1);
DLS_L_Gyr_norm_DFreq = ff(2);
DLS_L_Gyr_norm_PSD_mean = ff(3);
DLS_L_Gyr_norm_PSD_std = ff(4);
DLS_L_Gyr_norm_PSD_skew = ff(5);
DLS_L_Gyr_norm_PSD_kurtosis = ff(6);




