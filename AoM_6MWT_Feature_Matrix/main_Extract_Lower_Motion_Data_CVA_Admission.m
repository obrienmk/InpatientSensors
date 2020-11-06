clc
clear all
close all

% CVA
% ID: 1, 25  --> no data
% ID: 7, 40 --> sensor at DLS_L missing
% ID: 8 --> sensor at DLS_R missing
% ID: 32 --> long data 3 times*6MWT? 
% ID: 3, 22, 30 --> No admission sensor data 
% ID: 35, 43 --> Short data 5595 (179sec, 43), 10167 (325sec, 35)  frames 
% ID: 26, 54 --> Outlier. Sensor data broken


TS = [5 10 20 30 60 90 120 180 240 300 360]

for k = 1:1:length(TS)
    ID = [2 4:6 9:21 23:24 27:29 31 33:34 36:39 41:42 44:53 55]

    Type_of_Subject = 'CVA'

    for n = 1:1:length(ID)
        file_input = ['../Sensor_Data/' Type_of_Subject '_MWT6_ID' sprintf('%.2d',ID(n)) '.mat']
        load(file_input);

        Side = get_side('./CVA_Paretic_Side.csv');
        AS = Side.Side(ID(n));

        SN = 1;     % Session #
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
        Pel_Gyr_norm = sum(abs(SC_Gyr_norm)) * dt;  

        if strcmp(AS,'L') == 1
            Pel_Ro_AS = sum(abs(Pel_Ro_L)) * dt;
            Pel_Ro_US = sum(abs(Pel_Ro_R)) * dt; 
            Pel_Oblq_AS = sum(abs(Pel_Oblq_L)) * dt; 
            Pel_Oblq_US = sum(abs(Pel_Oblq_R)) * dt;

        elseif strcmp(AS,'R') == 1
            Pel_Ro_AS = sum(abs(Pel_Ro_R)) * dt;
            Pel_Ro_US = sum(abs(Pel_Ro_L)) * dt; 
            Pel_Oblq_AS = sum(abs(Pel_Oblq_R)) * dt; 
            Pel_Oblq_US = sum(abs(Pel_Oblq_L)) * dt;
        end
        

        % Ankle sensor
        if strcmp(AS,'L') == 1
            Gyr_Ankle_AS_x = sum(abs(DLS_L_Gyr(:,1))) * dt;  
            Gyr_Ankle_US_x = sum(abs(DLS_R_Gyr(:,1))) * dt;
            Gyr_Ankle_AS_y = sum(abs(DLS_L_Gyr(:,2))) * dt;  
            Gyr_Ankle_US_y = sum(abs(DLS_R_Gyr(:,2))) * dt;
            Gyr_Ankle_AS_z = sum(abs(DLS_L_Gyr(:,3))) * dt;  
            Gyr_Ankle_US_z = sum(abs(DLS_R_Gyr(:,3))) * dt;
            Gyr_Ankle_AS_norm = sum(abs(DLS_L_Gyr_norm)) * dt; 
            Gyr_Ankle_US_norm = sum(abs(DLS_R_Gyr_norm)) * dt;

        elseif strcmp(AS,'R') == 1
            Gyr_Ankle_AS_x = sum(abs(DLS_R_Gyr(:,1))) * dt;  
            Gyr_Ankle_US_x = sum(abs(DLS_L_Gyr(:,1))) * dt;
            Gyr_Ankle_AS_y = sum(abs(DLS_R_Gyr(:,2))) * dt;  
            Gyr_Ankle_US_y = sum(abs(DLS_L_Gyr(:,2))) * dt;
            Gyr_Ankle_AS_z = sum(abs(DLS_R_Gyr(:,3))) * dt;  
            Gyr_Ankle_US_z = sum(abs(DLS_L_Gyr(:,3))) * dt;
            Gyr_Ankle_AS_norm = sum(abs(DLS_R_Gyr_norm)) * dt; 
            Gyr_Ankle_US_norm = sum(abs(DLS_L_Gyr_norm)) * dt;
        end

        
        % Data out
        TSteps(n,:) = Step_R + Step_L;
        
        Gyr.AoM.Pel_Tilt(n,:) = [Pel_tilt_BW Pel_tilt_FW];
        Gyr.AoM.Pel_Ro(n,:) = [Pel_Ro_US Pel_Ro_AS];
        Gyr.AoM.Pel_Oblq(n,:) = [Pel_Oblq_US Pel_Oblq_AS];
        Gyr.AoM.Ankle_US(n,:) = [Gyr_Ankle_US_x Gyr_Ankle_US_y Gyr_Ankle_US_z];
        Gyr.AoM.Ankle_AS(n,:) = [Gyr_Ankle_AS_x Gyr_Ankle_AS_y Gyr_Ankle_AS_z];
        Gyr.AoM.norm(n,:) = [Pel_Gyr_norm Gyr_Ankle_US_norm Gyr_Ankle_AS_norm];

        Ltime(n,:) = length(Time);

    end

    clear data

    data.NoS = length(ID);
    data.ID = ID;
    data.Gyr = Gyr;
    data.Steps= TSteps;

    fileout = ['./CVA_Admission_6MWT_lower_motion' num2str(TS(k)) '.mat']
    save(fileout,'data')
    
    
    clear data Gyr Time NoSteps ...
         Pel_Tilt_BW Pel_Tilt_FW Pel_Ro_L Pel_Ro_R Pel_Oblq_L Pel_Oblq_R ...
         SC_Gyr DLS_L_Gyr DLS_R_Gyr SC_Gyr_norm DLS_L_Gyr_norm DLS_R_Gyr_norm ...
         TSteps step_R step_L DLS_L_Gyr_filt DLS_R_Gyr_filt pks_G_L locs_G_L pks_G_R locs_G_R Ltime 

end