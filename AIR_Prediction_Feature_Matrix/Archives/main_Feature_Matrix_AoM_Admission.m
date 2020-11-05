clc
clear all
close all

% CVA
% ID: 1, 3, 22, 25, 30, 32, 42, 44, 49  --> no data
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

    if n == 1 || n == 3 || n == 7 || n == 8 || n == 22 || n == 25 || n == 30 || n == 32 || n == 40 || n == 42 || n == 44 || n == 49
         SN = nan;     % Session #
         TN = nan;     % Trial #
        
        % Amount of Activity
        % Pelvis
        Gyr_Pel_tilt_BW(n,:) = nan;
        Gyr_Pel_tilt_FW(n,:) = nan;
        Acc_Pel_y_P(n,:) = nan;
        Acc_Pel_y_N(n,:) = nan;
        Acc_Pel_z_P(n,:) = nan;
        Acc_Pel_z_N(n,:) = nan;
        Gyr_Pel_norm(n,:) = nan;
        Acc_Pel_norm(n,:) = nan;

        Gyr_Pel_Ro_AS(n,:) = nan;
        Gyr_Pel_Ro_US(n,:) = nan; 
        Gyr_Pel_Oblq_AS(n,:) = nan;
        Gyr_Pel_Oblq_US(n,:) = nan;
        Acc_Pel_x_AS(n,:) = nan;
        Acc_Pel_x_US(n,:) = nan;

        Gyr_SI_Pel_Tilt(n,:) = nan;
        Gyr_SI_Pel_Ro(n,:) = nan;
        Gyr_SI_Pel_Oblq(n,:) = nan;
        Acc_SI_Pel_x(n,:) = nan;
        Acc_SI_Pel_y(n,:) = nan;
        Acc_SI_Pel_z(n,:) = nan;


        % Ankle sensor
        Gyr_Ankle_AS_x(n,:) = nan;
        Gyr_Ankle_US_x(n,:) = nan;
        Gyr_Ankle_AS_y(n,:) = nan;
        Gyr_Ankle_US_y(n,:) = nan;
        Gyr_Ankle_AS_z(n,:) = nan;
        Gyr_Ankle_US_z(n,:) = nan;
        Gyr_Ankle_AS_norm(n,:) = nan;
        Gyr_Ankle_US_norm(n,:) = nan;
        Acc_Ankle_AS_x(n,:) = nan;
        Acc_Ankle_US_x(n,:) = nan;
        Acc_Ankle_AS_y(n,:) = nan;
        Acc_Ankle_US_y(n,:) = nan;
        Acc_Ankle_AS_z(n,:) = nan;
        Acc_Ankle_US_z(n,:) = nan;
        Acc_Ankle_AS_norm(n,:) = nan;
        Acc_Ankle_US_norm(n,:) = nan;

        Gyr_SI_Ankle_x(n,:) = nan;
        Gyr_SI_Ankle_y(n,:) = nan;
        Gyr_SI_Ankle_z(n,:) = nan;
        Gyr_SI_Ankle_norm(n,:) = nan;
        Acc_SI_Ankle_x(n,:) = nan;
        Acc_SI_Ankle_y(n,:) = nan;
        Acc_SI_Ankle_z(n,:) = nan;
        Acc_SI_Ankle_norm(n,:) = nan;

        Ltime(n,:) = nan;
        
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

            % Acc
            if SC_Acc(i,1) >= 0
                Pel_x_R(i,:) = SC_Acc(i,1);
                Pel_x_L(i,:) = 0;
            elseif SC_Acc(i,1) < 0
                Pel_x_L(i,:) = SC_Acc(i,1);
                Pel_x_R(i,:) = 0;
            end

            if SC_Acc(i,2) >= 0
                Pel_y_P(i,:) = SC_Acc(i,2);
                Pel_y_N(i,:) = 0;
            elseif SC_Acc(i,2) < 0
                Pel_y_N(i,:) = SC_Acc(i,2);
                Pel_y_P(i,:) = 0;
            end

            if SC_Acc(i,3) >= 0
                Pel_z_P(i,:) = SC_Acc(i,3);
                Pel_z_N(i,:) = 0;
            elseif SC_Acc(i,3) < 0
                Pel_z_N(i,:) = SC_Acc(i,3);
                Pel_z_P(i,:) = 0;
            end

        end

        % Amount of Activity
        % Pelvis
        Gyr_Pel_tilt_BW(n,:) = sum(abs(Pel_Tilt_BW)) * dt;  
        Gyr_Pel_tilt_FW(n,:) = sum(abs(Pel_Tilt_FW)) * dt;
        Acc_Pel_y_P(n,:) = sum(abs(Pel_y_P)) * dt;
        Acc_Pel_y_N(n,:) = sum(abs(Pel_y_N)) * dt;
        Acc_Pel_z_P(n,:) = sum(abs(Pel_z_P)) * dt;
        Acc_Pel_z_N(n,:) = sum(abs(Pel_z_N)) * dt;
        Gyr_Pel_norm(n,:) = sum(abs(SC_Gyr_norm)) * dt;  
        Acc_Pel_norm(n,:) = sum(abs(SC_Acc_norm)) * dt; 

    %     if AS == 'L' 
        if strcmp(AS,'L') == 1
            Gyr_Pel_Ro_AS(n,:) = sum(abs(Pel_Ro_L)) * dt;
            Gyr_Pel_Ro_US(n,:) = sum(abs(Pel_Ro_R)) * dt; 
            Gyr_Pel_Oblq_AS(n,:) = sum(abs(Pel_Oblq_L)) * dt; 
            Gyr_Pel_Oblq_US(n,:) = sum(abs(Pel_Oblq_R)) * dt;
            Acc_Pel_x_AS(n,:) = sum(abs(Pel_x_L)) * dt;
            Acc_Pel_x_US(n,:) = sum(abs(Pel_x_R)) * dt;

        elseif strcmp(AS,'R') == 1
            Gyr_Pel_Ro_AS(n,:) = sum(abs(Pel_Ro_R)) * dt;
            Gyr_Pel_Ro_US(n,:) = sum(abs(Pel_Ro_L)) * dt; 
            Gyr_Pel_Oblq_AS(n,:) = sum(abs(Pel_Oblq_R)) * dt; 
            Gyr_Pel_Oblq_US(n,:) = sum(abs(Pel_Oblq_L)) * dt;
            Acc_Pel_x_AS(n,:) = sum(abs(Pel_x_R)) * dt;
            Acc_Pel_x_US(n,:) = sum(abs(Pel_x_L)) * dt;
        end

        Gyr_SI_Pel_Tilt(n,:) = (Gyr_Pel_tilt_BW(n,:) - Gyr_Pel_tilt_FW(n,:)) / (0.5 * (Gyr_Pel_tilt_BW(n,:) + Gyr_Pel_tilt_FW(n,:)));
        Gyr_SI_Pel_Ro(n,:) = (Gyr_Pel_Ro_US(n,:) - Gyr_Pel_Ro_AS(n,:)) / (0.5 * (Gyr_Pel_Ro_US(n,:) + Gyr_Pel_Ro_AS(n,:)));
        Gyr_SI_Pel_Oblq(n,:) = (Gyr_Pel_Oblq_US(n,:) - Gyr_Pel_Oblq_AS(n,:)) / (0.5 * (Gyr_Pel_Oblq_US(n,:) + Gyr_Pel_Oblq_AS(n,:)));
        Acc_SI_Pel_x(n,:) = (Acc_Pel_x_US(n,:) - Acc_Pel_x_AS(n,:)) / (0.5 * (Acc_Pel_x_US(n,:) + Acc_Pel_x_AS(n,:)));
        Acc_SI_Pel_y(n,:) = (Acc_Pel_y_P(n,:) - Acc_Pel_y_N(n,:)) / (0.5 * (Acc_Pel_y_P(n,:) + Acc_Pel_y_N(n,:)));
        Acc_SI_Pel_z(n,:) = (Acc_Pel_z_P(n,:) - Acc_Pel_z_N(n,:)) / (0.5 * (Acc_Pel_z_P(n,:) + Acc_Pel_z_N(n,:)));


        % Ankle sensor
        if strcmp(AS,'L') == 1
            Gyr_Ankle_AS_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;  
            Gyr_Ankle_US_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;
            Gyr_Ankle_AS_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;  
            Gyr_Ankle_US_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;
            Gyr_Ankle_AS_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;  
            Gyr_Ankle_US_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;
            Gyr_Ankle_AS_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt; 
            Gyr_Ankle_US_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt;
            Acc_Ankle_AS_x(n,:) = sum(abs(DLS_L_Acc(:,1))) * dt;  
            Acc_Ankle_US_x(n,:) = sum(abs(DLS_R_Acc(:,1))) * dt;
            Acc_Ankle_AS_y(n,:) = sum(abs(DLS_L_Acc(:,2))) * dt;  
            Acc_Ankle_US_y(n,:) = sum(abs(DLS_R_Acc(:,2))) * dt;
            Acc_Ankle_AS_z(n,:) = sum(abs(DLS_L_Acc(:,3))) * dt;  
            Acc_Ankle_US_z(n,:) = sum(abs(DLS_R_Acc(:,3))) * dt;
            Acc_Ankle_AS_norm(n,:) = sum(abs(DLS_L_Acc_norm)) * dt; 
            Acc_Ankle_US_norm(n,:) = sum(abs(DLS_R_Acc_norm)) * dt;

        elseif strcmp(AS,'R') == 1
            Gyr_Ankle_AS_x(n,:) = sum(abs(DLS_R_Gyr(:,1))) * dt;  
            Gyr_Ankle_US_x(n,:) = sum(abs(DLS_L_Gyr(:,1))) * dt;
            Gyr_Ankle_AS_y(n,:) = sum(abs(DLS_R_Gyr(:,2))) * dt;  
            Gyr_Ankle_US_y(n,:) = sum(abs(DLS_L_Gyr(:,2))) * dt;
            Gyr_Ankle_AS_z(n,:) = sum(abs(DLS_R_Gyr(:,3))) * dt;  
            Gyr_Ankle_US_z(n,:) = sum(abs(DLS_L_Gyr(:,3))) * dt;
            Gyr_Ankle_AS_norm(n,:) = sum(abs(DLS_R_Gyr_norm)) * dt; 
            Gyr_Ankle_US_norm(n,:) = sum(abs(DLS_L_Gyr_norm)) * dt;
            Acc_Ankle_AS_x(n,:) = sum(abs(DLS_R_Acc(:,1))) * dt;  
            Acc_Ankle_US_x(n,:) = sum(abs(DLS_L_Acc(:,1))) * dt;
            Acc_Ankle_AS_y(n,:) = sum(abs(DLS_R_Acc(:,2))) * dt;  
            Acc_Ankle_US_y(n,:) = sum(abs(DLS_L_Acc(:,2))) * dt;
            Acc_Ankle_AS_z(n,:) = sum(abs(DLS_R_Acc(:,3))) * dt;  
            Acc_Ankle_US_z(n,:) = sum(abs(DLS_L_Acc(:,3))) * dt;
            Acc_Ankle_AS_norm(n,:) = sum(abs(DLS_R_Acc_norm)) * dt; 
            Acc_Ankle_US_norm(n,:) = sum(abs(DLS_L_Acc_norm)) * dt;
        end

        Gyr_SI_Ankle_x(n,:) = (Gyr_Ankle_US_x(n,:) - Gyr_Ankle_AS_x(n,:)) / (0.5 * (Gyr_Ankle_US_x(n,:) + Gyr_Ankle_AS_x(n,:)));
        Gyr_SI_Ankle_y(n,:) = (Gyr_Ankle_US_y(n,:) - Gyr_Ankle_AS_y(n,:)) / (0.5 * (Gyr_Ankle_US_y(n,:) + Gyr_Ankle_AS_y(n,:)));
        Gyr_SI_Ankle_z(n,:) = (Gyr_Ankle_US_z(n,:) - Gyr_Ankle_AS_z(n,:)) / (0.5 * (Gyr_Ankle_US_z(n,:) + Gyr_Ankle_AS_z(n,:)));
        Gyr_SI_Ankle_norm(n,:) = (Gyr_Ankle_US_norm(n,:) - Gyr_Ankle_AS_norm(n,:)) / (0.5 * (Gyr_Ankle_US_norm(n,:) + Gyr_Ankle_AS_norm(n,:)));
        Acc_SI_Ankle_x(n,:) = (Acc_Ankle_US_x(n,:) - Acc_Ankle_AS_x(n,:)) / (0.5 * (Acc_Ankle_US_x(n,:) + Acc_Ankle_AS_x(n,:)));
        Acc_SI_Ankle_y(n,:) = (Acc_Ankle_US_y(n,:) - Acc_Ankle_AS_y(n,:)) / (0.5 * (Acc_Ankle_US_y(n,:) + Acc_Ankle_AS_y(n,:)));
        Acc_SI_Ankle_z(n,:) = (Acc_Ankle_US_z(n,:) - Acc_Ankle_AS_z(n,:)) / (0.5 * (Acc_Ankle_US_z(n,:) + Acc_Ankle_AS_z(n,:)));
        Acc_SI_Ankle_norm(n,:) = (Acc_Ankle_US_norm(n,:) - Acc_Ankle_AS_norm(n,:)) / (0.5 * (Acc_Ankle_US_norm(n,:) + Acc_Ankle_AS_norm(n,:)));

        Ltime(n,:) = length(Time)/Hz;
    end
    
end
ID = ID';
AS = Side.Side;
Time = Ltime;

CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_Discharge.csv');
FIM = CO(:,2);
BBS = CO(:,3);
MWT10_SSV = CO(:,4);
MWT10_FV = CO(:,5);
MWT6 = CO(:,6);
TUG = CO(:,7);


tbl_AoM = table(ID, AS, Time, ... % subject info
                Gyr_Pel_tilt_BW, Gyr_Pel_tilt_FW, Gyr_Pel_Oblq_US, Gyr_Pel_Oblq_AS, Gyr_Pel_Ro_US, Gyr_Pel_Ro_AS, Gyr_Pel_norm, ...  gyro pelvis
                Gyr_Ankle_US_x, Gyr_Ankle_AS_x, Gyr_Ankle_US_y, Gyr_Ankle_AS_y, Gyr_Ankle_US_z, Gyr_Ankle_AS_z, Gyr_Ankle_US_norm, Gyr_Ankle_AS_norm, ...    gyro foot
                Acc_Pel_x_US, Acc_Pel_x_AS, Acc_Pel_y_P, Acc_Pel_y_N, Acc_Pel_z_P, Acc_Pel_z_N, Acc_Pel_norm, ...    acc pelvis
                Acc_Ankle_US_x, Acc_Ankle_AS_x, Acc_Ankle_US_y, Acc_Ankle_AS_y, Acc_Ankle_US_z, Acc_Ankle_AS_z, Acc_Ankle_US_norm, Acc_Ankle_AS_norm, ...     acc foot
                Gyr_SI_Pel_Tilt, Gyr_SI_Pel_Oblq, Gyr_SI_Pel_Ro, Acc_SI_Pel_x, Acc_SI_Pel_y, Acc_SI_Pel_z, ... SI Pel
                Gyr_SI_Ankle_x, Gyr_SI_Ankle_y, Gyr_SI_Ankle_z, Gyr_SI_Ankle_norm, Acc_SI_Ankle_x, Acc_SI_Ankle_y, Acc_SI_Ankle_z, Acc_SI_Ankle_norm, ...   SI foot
                FIM, BBS, MWT10_SSV, MWT10_FV, MWT6, TUG ... Clincal outcomes
                )   

writetable(tbl_AoM,'Feature_Matrix_AoM_Admission_6MWT.csv','Delimiter',',','QuoteStrings',true)



