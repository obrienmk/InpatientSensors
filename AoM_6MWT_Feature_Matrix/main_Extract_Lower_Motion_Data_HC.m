clc
clear all
close all

% Healthy Controls (HC)
% ID: 26, 38 --> no data
% ID: 4 --> no sensor at DLS_L

TS = [5 10 20 30 60 90 120 180 240 300 360]

for k = 1:1:length(TS)
    ID = [1:3 5:25 27:37 39:51]
    
    Type_of_Subject = 'CONTROLS'

    for n = 1:1:length(ID)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID(n)) '.mat']
        load(file_input);

        SN = length(data.Session);     % Session #
        TN = 1;     % Trial #

        Hz = 31.25;
        dt = 1/Hz;

        % Gyroscope data
        SC_Gyr = data.Session{SN}.Motion.SC.Gyr{TN};
        DLS_R_Gyr = data.Session{SN}.Motion.DLS_R.Gyr{TN};
        DLS_L_Gyr = data.Session{SN}.Motion.DLS_L.Gyr{TN};


        % To match final time
%         final = min([length(SC_Gyr) length(SC_Acc) length(DLS_R_Gyr) length(DLS_R_Acc) length(DLS_L_Gyr) length(DLS_L_Acc)]);
        final = floor(TS(k)*Hz);
        Time = data.Session{SN}.Motion.Time{TN}(1:final,:);
        SC_Gyr = SC_Gyr(1:final,:);
        DLS_R_Gyr = DLS_R_Gyr(1:final,:);
        DLS_L_Gyr = DLS_L_Gyr(1:final,:);

        % Norm of Gyro and Acc Data
        for i = 1:1:length(Time)
            SC_Gyr_norm(i,:) = norm(SC_Gyr(i,:));
            DLS_R_Gyr_norm(i,:) = norm(DLS_R_Gyr(i,:));
            DLS_L_Gyr_norm(i,:) = norm(DLS_L_Gyr(i,:));
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


        % AoM
        for i = 1:1:length(Time)
            % Gyro
            if SC_Gyr(i,1) >= 0
                Pel_Tilt_BW(i,:) = SC_Gyr(i,1);
                Pel_Tilt_FW(i,:) = 0;
            elseif SC_Gyr(i,1) < 0
                Pel_Tilt_FW(i,:) = SC_Gyr(i,1);
                Pel_Tilt_BW(i,:) = 0;
            end

            if SC_Gyr(i,2) >= 0
                Pel_Ro_L(i,:) = SC_Gyr(i,2);
                Pel_Ro_R(i,:) = 0;
            elseif SC_Gyr(i,2) < 0
                Pel_Ro_R(i,:) = SC_Gyr(i,2);
                Pel_Ro_L(i,:) = 0;
            end

            if SC_Gyr(i,3) >= 0
                Pel_Oblq_L(i,:) = SC_Gyr(i,3);
                Pel_Oblq_R(i,:) = 0;
            elseif SC_Gyr(i,3) < 0
                Pel_Oblq_R(i,:) = SC_Gyr(i,3);
                Pel_Oblq_L(i,:) = 0;
            end

        end

        % Amount of Activity
        % Pelvis
        Pel_tilt_BW = sum(abs(Pel_Tilt_BW)) * dt;  
        Pel_tilt_FW = sum(abs(Pel_Tilt_FW)) * dt; 
        Pel_Ro_L = sum(abs(Pel_Ro_L)) * dt;
        Pel_Ro_R = sum(abs(Pel_Ro_R)) * dt; 
        Pel_Oblq_L = sum(abs(Pel_Oblq_L)) * dt; 
        Pel_Oblq_R = sum(abs(Pel_Oblq_R)) * dt;
        Pel_Gyr_norm = sum(abs(SC_Gyr_norm)) * dt;  


        % Ankle sensor
        Gyr_Ankle_L_x = sum(abs(DLS_L_Gyr(:,1))) * dt;  
        Gyr_Ankle_R_x = sum(abs(DLS_R_Gyr(:,1))) * dt;
        Gyr_Ankle_L_y = sum(abs(DLS_L_Gyr(:,2))) * dt;  
        Gyr_Ankle_R_y = sum(abs(DLS_R_Gyr(:,2))) * dt;
        Gyr_Ankle_L_z = sum(abs(DLS_L_Gyr(:,3))) * dt;  
        Gyr_Ankle_R_z = sum(abs(DLS_R_Gyr(:,3))) * dt;
        Gyr_Ankle_L_norm = sum(abs(DLS_L_Gyr_norm)) * dt; 
        Gyr_Ankle_R_norm = sum(abs(DLS_R_Gyr_norm)) * dt;


        % Data out
        TSteps(n,:) = Step_R + Step_L;

        Gyr.AoM.Pel_Tilt(n,:) = [Pel_tilt_BW Pel_tilt_FW];
        Gyr.AoM.Pel_Ro(n,:) = [Pel_Ro_L Pel_Ro_R];
        Gyr.AoM.Pel_Oblq(n,:) = [Pel_Oblq_L Pel_Oblq_R];
        Gyr.AoM.Ankle_L(n,:) = [Gyr_Ankle_L_x Gyr_Ankle_L_y Gyr_Ankle_L_z];
        Gyr.AoM.Ankle_R(n,:) = [Gyr_Ankle_R_x Gyr_Ankle_R_y Gyr_Ankle_R_z];
        Gyr.AoM.norm(n,:) = [Pel_Gyr_norm Gyr_Ankle_L_norm Gyr_Ankle_R_norm];

        Ltime(n,:) = length(Time);

    end
    
    clear data

    data.NoS = length(ID);
    data.ID = ID;
    data.Gyr = Gyr;
    data.Steps= TSteps;

    fileout = ['./HC_6MWT_lower_motion' num2str(TS(k)) '.mat']
    save(fileout,'data')
    
    
    clear data Gyr Time NoSteps ...
         Pel_Tilt_BW Pel_Tilt_FW Pel_Ro_L Pel_Ro_R Pel_Oblq_L Pel_Oblq_R ...
         SC_Gyr DLS_L_Gyr DLS_R_Gyr SC_Gyr_norm DLS_L_Gyr_norm DLS_R_Gyr_norm ...
         TSteps step_R step_L DLS_L_Gyr_filt DLS_R_Gyr_filt pks_G_L locs_G_L pks_G_R locs_G_R Ltime 

     
end
