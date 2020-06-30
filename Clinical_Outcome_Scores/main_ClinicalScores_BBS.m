%-------------------------------------------------------------------
%
% Extract BBS Total Scores and Subscores for Inpatient Sensor Data
% Adam P. Horin
% June 29, 2020
%
%------------------------------------------------------------------

clc
clear all
close all

%% Process BBS Clinical Scores
% Define Variables
Type_of_Subject = {'CVA', 'CONTROLS'}; % enter CONTROLS and/or CVA
BBS = 'BBS';

TrialNames = {'1. SIT TO STAND', ...
    '2. STAND UNSUPPORTED', ...
    '3. SIT W/ BACK', ...
    '4. STAND TO SIT', ...
    '5. TRANSFERS', ...
    '6. STAND W/ EYES CLOSED', ...
    '7. STAND W/ FEET TOGETHER', ...
    '8. REACH FORWARD', ...
    '9. PICK UP OBJECT', ...
    '10. TURN', ...
    '11. TURN 360 DEGS', ...
    '12. FOOT ON STEP/STOOL', ...
    '13. STAND W/ ONE FOOT IN FRONT', ...
    '14. STAND ON ONE LEG'};

TrialNumbers = {'N01', ...
    'N02', ...
    'N03', ...
    'N04', ...
    'N05', ...
    'N06', ...
    'N07', ...
    'N08', ...
    'N09', ...
    'N10', ...
    'N11', ...
    'N12', ...
    'N13', ...
    'N14'};

BBS_Scores = [];
List_WrongDate = [];

%% Main Loop
% Identify group to loop through
for h = 1:1:length(Type_of_Subject)
    if strcmp(Type_of_Subject{h}, 'CVA') == 1
        ID = [1]; %[1:55] %53 has an extra row of sensor measurements
        Type_of_Subject_Group = 'CVA';
        BBS_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_BBS.xlsx';
    %test date order with mock data BBS_path = "C:\Users\ahorin\Documents\Inpatient Sensor Code AH\General_Feature_Matrix_Extraction_Code\mock_BBS_score.xlsx";
    elseif strcmp(Type_of_Subject{h}, 'CONTROLS') == 1
        ID = []; %[1:51]
        Type_of_Subject_Group = 'CONTROLS';
        BBS_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_BBS.xlsx';
    end
    
    % loop for subject ID
    for n = 1:1:length(ID)
        
        %Get Subject ID
        if isnumeric(ID(n))
            str = num2str(ID(n));
            if length(str) == 1
                if strcmp(Type_of_Subject_Group,'CVA')
                    SubjectID=['CVA0' str];
                else strcmp(Type_of_Subject_Group,'CONTROLS')
                    SubjectID=['HC0' str];
                end
            else
                if strcmp(Type_of_Subject_Group,'CVA')
                    SubjectID=['CVA' str];
                else strcmp(Type_of_Subject_Group,'CONTROLS')
                    SubjectID=['HC' str];
                end
            end
        end
        
        % Load BBS scores from excel file and sheet
        [Sig, TStr, Raw] = xlsread(BBS_path, SubjectID)
        Raw([1],:) = []; % remove the header row from the dates
        Raw(strcmp(Raw, 'DNA')) = {NaN}; % replace 'DNA' with nan
        RawScore = cell2mat(Raw(:,3)); % get the raw scores from the Raw cells to make indexing easier
        %dimBBS = size(Raw);
        
        % Get the index of each group of dates for each session
        SN1_start=1;
        for k = 1:1:length(Raw)
            if datenum(Raw(SN1_start)) == datenum(Raw(k))
                SN1_end = SN1_start+k-1;
                SN2_start = SN1_end+1;
                SN1 = SN1_start:SN1_end;
                SN2 = [];
                SN3 = [];
                SN4 = [];
            elseif datenum(Raw(SN2_start)) == datenum(Raw(k))
                SN2_end = SN1_start+k-1;
                SN3_start = SN2_end+1;
                SN2 = [SN2_start:SN2_end];
                SN3 = [];
                SN4 = [];
            elseif datenum(Raw(SN3_start)) == datenum(Raw(k))
                SN3_end = SN1_start+k-1;
                SN4_start = SN3_end+1;
                SN3 = [SN3_start:SN3_end];
                SN4 = [];
            elseif datenum(Raw(SN4_start)) == datenum(Raw(k))
                SN4_end = SN1_start+k-1;
                SN4 = [SN4_start:SN4_end];
            end
            
        end
        SN = {SN1, SN2, SN3, SN4};
        
        
        % check to make sure the the dates are in the correct order
        % make the last session SN4
        
        if isempty(SN2) == 1 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(Raw(SN1_start))];
            order = 1;
        elseif isempty(SN2) == 0 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(Raw(SN1_start)), datenum(Raw(SN2_start))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                SN4 = SN2; %makes the discharge session SN4
                SN2 = [];
                order = 1;
            else
                %dates are not in correct order
                %date_index = find(contains(dates(),dates_sorted(1))))
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                SN1 = SN{date_index1};
                SN2 = [];
                SN3 = [];
                SN4 = SN{date_index2};
                order = 0;
            end
        elseif isempty(SN2) == 0 && isempty(SN3) == 0 && isempty(SN4) == 1
            
            dates = [datenum(Raw(SN1_start)), datenum(Raw(SN2_start)), datenum(Raw(SN3_start))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                SN4 = SN3;
                SN3 = [];
                order = 1;
            else
                %dates are not in correct order
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                date_index3 = find(ismember(dates, dates_sorted(3)));
                SN1 = SN{date_index1};
                SN2 = SN{date_index2};
                SN3 = [];
                SN4 = SN{date_index3};
                order = 0;
                
            end
        elseif isempty(SN2) == 0 && isempty(SN3) == 0 && isempty(SN4) == 0
            dates = [datenum(Raw(SN1_start)), datenum(Raw(SN2_start)), datenum(Raw(SN3_start)), datenum(Raw(SN4_start))];
            dates_sorted = sort(dates);
            if dates == dates_sorted
                %dates in the correct order
                order = 1;
            else
                %dates are not in correct order
                date_index1 = find(ismember(dates, dates_sorted(1)));
                date_index2 = find(ismember(dates, dates_sorted(2)));
                date_index3 = find(ismember(dates, dates_sorted(3)));
                date_index4 = find(ismember(dates, dates_sorted(4)));
                SN1 = SN{date_index1};
                SN2 = SN{date_index2};
                SN3 =  SN{date_index3};
                SN4 =  SN{date_index4};
                order = 0;
            end
        end
        
        
        % Redefine SN in case dates were out of order
        SN = {SN1, SN2, SN3, SN4};
        
        % Loop Through Sessions
        % SN1 = Admission, SN4 = Discharge, SN2 and SN3 are any sessions additional
        
        
        for j = 1:1:length(SN)
            
            for g = 1:1:length(TrialNames)
                if isempty(SN{j}) == 1
                    subject(g,:) = {SubjectID};
                    group(g,:) = Type_of_Subject(h);
                    session(g,:) = j;
                    activity(g,:) = BBS;
                    trial_no(g,:) = TrialNumbers(g);
                    BBS_TotalScore(g,:) = nan;
                    BBS_Subscore(g,:) = nan;
                else
                    BBS_index = find(contains(Raw(SN{j},2), TrialNames(g))) %use contains() in case trials were repeated; will select the first attempt
                    
                    subject(g,:) = {SubjectID};
                    group(g,:) = Type_of_Subject(h);
                    session(g,:) = j;
                    activity(g,:) = BBS;
                    trial_no(g,:) = TrialNumbers(g);
                    BBS_TotalScore(g,:) = nansum(RawScore(SN{j}));
                    BBS_Subscore(g,:) = RawScore(SN{j}(BBS_index(1))); %use the first attempt of any repeated trials
                end
            end
            
            
            
            
            
            % BBS Clinical Scores Matrix
            BBS_Scores_temp = table(subject, group, session, activity, trial_no, BBS_TotalScore, BBS_Subscore);
            BBS_Scores = [BBS_Scores; BBS_Scores_temp]
            
            % Check if any date orders needed to be changed
            if order(1) == 0
                DateOrderID = {SubjectID};
                WrongSession = j;
                List_WrongDate_temp = table(DateOrderID, WrongSession);
                List_WrongDate = [List_WrongDate; List_WrongDate_temp];
            end
            
            
        end
    end
end

writetable(BBS_Scores,'BBS_ClinicalScores.csv','Delimiter',',','QuoteStrings',true)



