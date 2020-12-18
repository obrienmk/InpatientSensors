clc
clear all
close all

Sub_ID = [];
Sub_Type = [];
FIM = [];
BBS = [];
MWT10_SSV = [];
MWT10_FV = [];
MWT6 = [];
TUG = [];
AoM_Pel_tilt = [];
AoM_Pel_ro = [];
AoM_Pel_oblq = [];
AoM_Ankle_US_x = [];
AoM_Ankle_US_y = [];
AoM_Ankle_US_z = [];
AoM_Ankle_AS_x = [];
AoM_Ankle_AS_y = [];
AoM_Ankle_AS_z = [];
AoM_Pel_Norm = [];
AoM_Ankle_US_Norm = [];
AoM_Ankle_AS_Norm = [];
Steps = [];
Time = [];

TS = [5 10 20 30 60 90 120 180 240 300 360]
for k = 1:1:length(TS)
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

    CO = import_Clincal_Outcomes_Discharge('./CVA_Clinical_Outcome_First.csv');
    % Baseline Clinical Outcomes
    Out_ID = [54 43 40 35 32 30 26 25 22 8 7 3 1];

    for i = 1:1:length(Out_ID)
        CO(Out_ID(i),:) = [];
    end

    FIM_AD = CO(:,2);
    BBS_AD = CO(:,3);
    MWT10_SSV_AD = CO(:,4);
    MWT10_FV_AD = CO(:,5);
    MWT6_AD = CO(:,6);
    TUG_AD = CO(:,7);

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


    FIM_ = [FIM_AD; FIM_DC; FIM_HC];
    BBS_ = [BBS_AD; BBS_DC; BBS_HC];
    MWT10_SSV_ = [MWT10_SSV_AD; MWT10_SSV_DC; MWT10_SSV_HC];
    MWT10_FV_ = [MWT10_FV_AD; MWT10_FV_DC; MWT10_FV_HC];
    MWT6_ = [MWT6_AD; MWT6_DC; MWT6_HC];
    TUG_ = [TUG_AD; TUG_DC; TUG_HC];


    file_AD = ['./CVA_Admission_6MWT_lower_motion' num2str(TS(k)) '.mat']
    load(file_AD);
    CVA_AD = data;

    file_DC = ['./CVA_Discharge_6MWT_lower_motion' num2str(TS(k)) '.mat']
    load(file_DC);
    CVA = data

    file_HC = ['./HC_6MWT_lower_motion' num2str(TS(k)) '.mat']
    load(file_HC);
    HC = data;



    % CVA Admission: AoM
    Steps_AD = CVA_AD.Steps;

    Pel_tilt_CVA_AD = sum(CVA_AD.Gyr.AoM.Pel_Tilt')';
    Pel_ro_CVA_AD = sum(CVA_AD.Gyr.AoM.Pel_Ro')';
    Pel_oblq_CVA_AD = sum(CVA_AD.Gyr.AoM.Pel_Oblq')';
    Ankle_US_AD = CVA_AD.Gyr.AoM.Ankle_US;
    Ankle_AS_AD = CVA_AD.Gyr.AoM.Ankle_AS;

    Pel_Norm_CVA_AD = CVA_AD.Gyr.AoM.norm(:,1);
    AoM_Ankle_US_Norm_AD = CVA_AD.Gyr.AoM.norm(:,2);
    AoM_Ankle_AS_Norm_AD = CVA_AD.Gyr.AoM.norm(:,3);

    AoM_Ankle_US_x_AD = Ankle_US_AD(:,1);
    AoM_Ankle_US_y_AD = Ankle_US_AD(:,2);
    AoM_Ankle_US_z_AD = Ankle_US_AD(:,3);
    AoM_Ankle_AS_x_AD = Ankle_AS_AD(:,1);
    AoM_Ankle_AS_y_AD = Ankle_AS_AD(:,2);
    AoM_Ankle_AS_z_AD = Ankle_AS_AD(:,3);

    % CVA: AoM
    Steps_CVA = CVA.Steps;

    Pel_tilt_CVA = sum(CVA.Gyr.AoM.Pel_Tilt')';
    Pel_ro_CVA = sum(CVA.Gyr.AoM.Pel_Ro')';
    Pel_oblq_CVA = sum(CVA.Gyr.AoM.Pel_Oblq')';
    Ankle_US = CVA.Gyr.AoM.Ankle_US;
    Ankle_AS = CVA.Gyr.AoM.Ankle_AS;

    Pel_Norm_CVA = CVA.Gyr.AoM.norm(:,1);
    AoM_Ankle_US_Norm_CVA = CVA.Gyr.AoM.norm(:,2);
    AoM_Ankle_AS_Norm_CVA = CVA.Gyr.AoM.norm(:,3);

    AoM_Ankle_US_x_CVA = Ankle_US(:,1);
    AoM_Ankle_US_y_CVA = Ankle_US(:,2);
    AoM_Ankle_US_z_CVA = Ankle_US(:,3);
    AoM_Ankle_AS_x_CVA = Ankle_AS(:,1);
    AoM_Ankle_AS_y_CVA = Ankle_AS(:,2);
    AoM_Ankle_AS_z_CVA = Ankle_AS(:,3);

    % HC: AoM
    Steps_HC = HC.Steps;

    Pel_tilt_HC = sum(HC.Gyr.AoM.Pel_Tilt')';
    Pel_ro_HC = sum(HC.Gyr.AoM.Pel_Ro')';
    Pel_oblq_HC = sum(HC.Gyr.AoM.Pel_Oblq')';
    Ankle_L = HC.Gyr.AoM.Ankle_L;
    Ankle_R = HC.Gyr.AoM.Ankle_R;

    Pel_Norm_HC = HC.Gyr.AoM.norm(:,1);
    AoM_Ankle_L_Norm = HC.Gyr.AoM.norm(:,2);
    AoM_Ankle_R_Norm = HC.Gyr.AoM.norm(:,3);

    AoM_Ankle_L_x = Ankle_L(:,1);
    AoM_Ankle_L_y = Ankle_L(:,2);
    AoM_Ankle_L_z = Ankle_L(:,3);
    AoM_Ankle_R_x = Ankle_R(:,1);
    AoM_Ankle_R_y = Ankle_R(:,2);
    AoM_Ankle_R_z = Ankle_R(:,3);

    % CVA+HC AoM
    Steps_ = [Steps_AD; Steps_CVA; Steps_HC];
    AoM_Pel_tilt_ = [Pel_tilt_CVA_AD; Pel_tilt_CVA; Pel_tilt_HC];
    AoM_Pel_oblq_ = [Pel_oblq_CVA_AD; Pel_oblq_CVA; Pel_oblq_HC];
    AoM_Pel_ro_ = [Pel_ro_CVA_AD; Pel_ro_CVA; Pel_ro_HC];
    AoM_Pel_Norm_ = [Pel_Norm_CVA_AD; Pel_Norm_CVA; Pel_Norm_HC];

    AoM_Ankle_US_x_ = [AoM_Ankle_US_x_AD; AoM_Ankle_US_x_CVA; AoM_Ankle_R_x];
    AoM_Ankle_US_y_ = [AoM_Ankle_US_y_AD; AoM_Ankle_US_y_CVA; AoM_Ankle_R_y];
    AoM_Ankle_US_z_ = [AoM_Ankle_US_z_AD; AoM_Ankle_US_z_CVA; AoM_Ankle_R_z];
    AoM_Ankle_AS_x_ = [AoM_Ankle_AS_x_AD; AoM_Ankle_AS_x_CVA; AoM_Ankle_L_x];
    AoM_Ankle_AS_y_ = [AoM_Ankle_AS_y_AD; AoM_Ankle_AS_y_CVA; AoM_Ankle_L_y];
    AoM_Ankle_AS_z_ = [AoM_Ankle_AS_z_AD; AoM_Ankle_AS_z_CVA; AoM_Ankle_L_z];
    AoM_Ankle_US_Norm_ = [AoM_Ankle_US_Norm_AD; AoM_Ankle_US_Norm_CVA; AoM_Ankle_R_Norm];
    AoM_Ankle_AS_Norm_ = [AoM_Ankle_AS_Norm_AD; AoM_Ankle_AS_Norm_CVA; AoM_Ankle_L_Norm];

    ID = [CVA_AD.ID CVA.ID HC.ID];

    for i = 1:1:length(FIM_)
        if i >= 1 && i < 43
            Sub_Type_{i} = 'CVA_AD'
            Sub_ID_{i} = ['CVA' num2str(ID(i))];
        elseif i >= 43 && i < 93
            Sub_Type_{i} = 'CVA_DC'
            Sub_ID_{i} = ['CVA' num2str(ID(i))];
        elseif i >= 93 && i <= 140
            Sub_Type_{i} = 'HC'
            Sub_ID_{i} = ['HC' num2str(ID(i))];
        end
    end

    Sub_ID_ = Sub_ID_';
    Sub_Type_ = Sub_Type_';
    Time_ = TS(k)*ones(140,1);
    
    Sub_ID = [Sub_ID; Sub_ID_];
    Sub_Type = [Sub_Type; Sub_Type_];
    Time = [Time; Time_];
    FIM = [FIM; FIM_];
    BBS = [BBS; BBS_];
    MWT10_SSV = [MWT10_SSV; MWT10_SSV_];
    MWT10_FV = [MWT10_FV; MWT10_FV_];
    MWT6 = [MWT6; MWT6_];
    TUG = [TUG; TUG_];
    AoM_Pel_tilt = [AoM_Pel_tilt; AoM_Pel_tilt_];
    AoM_Pel_ro = [AoM_Pel_ro; AoM_Pel_ro_];
    AoM_Pel_oblq = [AoM_Pel_oblq; AoM_Pel_oblq_];
    AoM_Ankle_US_x = [AoM_Ankle_US_x; AoM_Ankle_US_x_];
    AoM_Ankle_US_y = [AoM_Ankle_US_y; AoM_Ankle_US_y_];
    AoM_Ankle_US_z = [AoM_Ankle_US_z; AoM_Ankle_US_z_];
    AoM_Ankle_AS_x = [AoM_Ankle_AS_x; AoM_Ankle_AS_x_];
    AoM_Ankle_AS_y = [AoM_Ankle_AS_y; AoM_Ankle_AS_y_];
    AoM_Ankle_AS_z = [AoM_Ankle_AS_z; AoM_Ankle_AS_z_];
    AoM_Pel_Norm = [AoM_Pel_Norm; AoM_Pel_Norm_];
    AoM_Ankle_US_Norm = [AoM_Ankle_US_Norm; AoM_Ankle_US_Norm_];
    AoM_Ankle_AS_Norm = [AoM_Ankle_AS_Norm; AoM_Ankle_AS_Norm_];
    Steps = [Steps; Steps_];

    
    clear Sub_Type_ Sub_ID_ Time_
end

tbl_AoM = table(Sub_ID, Sub_Type, Time, FIM, BBS, MWT10_SSV, MWT10_FV, MWT6, TUG, ...
    AoM_Pel_tilt, AoM_Pel_ro, AoM_Pel_oblq, ...
    AoM_Ankle_US_x, AoM_Ankle_US_y, AoM_Ankle_US_z, AoM_Ankle_AS_x, AoM_Ankle_AS_y, AoM_Ankle_AS_z, ...
    AoM_Pel_Norm, AoM_Ankle_US_Norm, AoM_Ankle_AS_Norm, ...
    Steps);

fileout = ['./MWT6_AoM_Feature_Matrix.csv']

writetable(tbl_AoM,fileout)   