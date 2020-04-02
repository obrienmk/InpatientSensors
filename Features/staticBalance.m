function features=staticBalance(transformData,t)

% Funtion staticBalance = computes static balance features form the Berg
% Balance Scale Clinical Test.

% 1. Standing Unsupported
% 2. Standing Eyes Closed
% 3. Standing Feet Together
% 4. Tamdem Standing
% 5. Standing One Leg

% Some Citations
% [1]. S. Del Din, A. Godfrey, S. Coleman, B. Galna, S. Lord, and L. Rochester,
% “Time-dependent changes in postural control in early Parkinson’s disease:
% what are we missing?,” Med. Biol. Eng. Comput., vol. 54, no. 2–3,
% pp. 401–410, 2016.

% Inputs
% t = time series vector;
% transformData = structure of the transformed accelerations by the coordinate transformation
% method

% Outputs
% features = time and ffrequency domain features

%Acceleration (Remove Linear Trend)
aAP=detrend(transformData.AP,'linear');
aML=detrend(transformData.ML,'linear');

%% Butterworth Filter 4th order

fs = 1/(t(2)-t(1));       % sampling frequency (Hz)
fc=3.5;                         % Low pass at 3.5Hz
[b,a] = butter(4,fc/(fs/2));
aAPFilt = filtfilt(b,a,aAP);    % Filter Signal AP
aMLFilt = filtfilt(b,a,aML);    % Filter Signal AP

%% FFT (AP & ML Acceleration signals)

n = length(t);                  % number of samples
f = (1:n/2+1)*(fs/n);           % frequency range
delta = 1/fs;                   % sampling period (s)

% Power spectrum
accML_fft = fft(aMLFilt);  % FFT
accML_pow =  1/(n*fs)*abs(accML_fft(1:floor(n/2)+1)).^2;   % power of the FFT
accML_pow=accML_pow';       %[:,1]
accAP_fft = fft(aAPFilt);  % FFT
accAP_pow =  abs(accAP_fft(1:floor(n/2)+1)).^2;  % power of the FFT
accAP_pow=accAP_pow';       %[:,1]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Frequency Domain Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_ML=cumtrapz(f,accML_pow);
area_AP=cumtrapz(f,accAP_pow);

Ptotal_ML=area_ML(end);
Ptotal_AP=area_AP(end);

Area_ML = area_ML / Ptotal_ML;
Area_AP = area_AP / Ptotal_AP;

f95=0.95;  f50=0.50;
[~, ind95ML] = min(abs(Area_ML-f95));
[~, ind50ML] = min(abs(Area_ML-f50));

[~, ind95AP] = min(abs(Area_AP-f95));
[~, ind50AP] = min(abs(Area_AP-f50));

%A. Frequency up to 95% Power(F95)
f95_ML=f(ind95ML);
f95_AP=f(ind95AP);

%B.Frequency up tospectral_centroid_AP 50% Power(F50)
f50_ML=f(ind50ML);
f50_AP=f(ind50AP);

%C.Centroidal Frequency (CF)
spectral_centroid_AP = sum(f*(accAP_pow))/sum((accAP_pow));
spectral_centroid_ML = sum(f*(accML_pow))/sum((accML_pow));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Time Domain Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A.Root Mean Square ACCELERATION (RMS)
rms_AP = rms(aAPFilt);
rms_ML = rms(aMLFilt);

%B.Ellipsoid 95% ML & AP
data=[aML' aAP']; %2 columns
[center,angles,R,area,axis]=ellipsoid(data);
axis=axis';
%D.Mean Sway Velocity (magnitude)
mean_velAP=mean(abs(cumtrapz(t,aAPFilt)));
mean_velML=mean(abs(cumtrapz(t,aMLFilt)));

%F.Path Length of Acceleration
length_swayAP=sum(abs(diff(aAPFilt',1,1)));  % path length of acceleration (m/s^2)
length_swayML=sum(abs(diff(aMLFilt',1,1)));

%G.Jerk (magnitude)
jerk_AP=mean(abs(diff(aAPFilt')./diff(t)));
jerk_ML=mean(abs(diff(aMLFilt')./diff(t)));

%E.Mean Acceleration (magnitude)
mean_accAP=mean(abs(aAPFilt));
mean_accML=mean(abs(aMLFilt));

%C.Maximum Acceleration (magnitude)
max_accAP=max(abs(aAPFilt));
max_accML=max(abs(aMLFilt));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Store Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Features Frequency Domain
features.data.f50_ML=f50_ML;
features.data.f50_AP=f50_AP;
features.data.f95_ML=f95_ML;
features.data.f95_AP=f95_AP;
features.data.spectral_centroid_AP=spectral_centroid_AP;
features.data.spectral_centroid_ML=spectral_centroid_ML;

% Features Time Domain
features.data.max_accAP = max_accAP;
features.data.max_accML = max_accML;
features.data.mean_accAP = mean_accAP;
features.data.mean_accML = mean_accML;

features.data.rms_AP=rms_AP;
features.data.rms_ML=rms_ML;

features.data.ellipAngleY=angles(1);
features.data.ellipAngleX=angles(2);
features.data.ellipArea=area;
features.data.ellipAxisY=axis(1);
features.data.ellipAxisX=axis(2);

features.data.jerk_AP=jerk_AP;
features.data.jerk_ML=jerk_ML;

features.data.mean_velAP=mean_velAP;
features.data.mean_velML=mean_velML;

features.data.length_swayAPAcc=length_swayAP;
features.data.length_swayMLAcc=length_swayML;



% Plot Info
features.plot.ML.fftmag=accML_pow;
features.plot.ML.f50=f50_ML;
features.plot.ML.f95=f95_ML;
features.plot.ML.spectral_centroid=spectral_centroid_ML;
features.plot.ML.f=f;
features.plot.ML.signal=aML';

features.plot.AP.fftmag=accAP_pow;
features.plot.AP.f50=f50_AP;
features.plot.AP.f95=f95_AP;
features.plot.AP.spectral_centroid=spectral_centroid_AP;
features.plot.AP.f=f;
features.plot.t=t;
features.plot.AP.signal=aAP';

features.plot.center=center;
features.plot.angles=angles;
features.plot.R=R;
features.plot.axis=axis;
end