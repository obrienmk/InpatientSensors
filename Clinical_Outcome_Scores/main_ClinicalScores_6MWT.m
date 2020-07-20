%-------------------------------------------------------------------
%
% Extract 6MWT Scores for Inpatient Sensor Data
% Adam P. Horin
% June 29, 2020
%
%------------------------------------------------------------------

clc
clear all
close all

%% Process 6MWT Clinical Scores
% Define Variables
Type_of_Subject = {'CVA', 'CONTROLS'}; % enter CONTROLS and/or CVA
Activity = '6MWT';

Minute = [1,2,3,4,5,6];

MWT6_Scores = [];
List_WrongDate = [];

%% Main Loop
% Identify group to loop through
for h = 1:1:length(Type_of_Subject)
    if strcmp(Type_of_Subject{h}, 'CVA') == 1
        ID = [1:55]; %[1:55]
        Type_of_Subject_Group = 'CVA';
        MWT6_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_6MWT.xlsx';
    elseif strcmp(Type_of_Subject{h}, 'CONTROLS') == 1
        ID = [1:51]; %[1:51]
        Type_of_Subject_Group = 'CONTROLS';
        MWT6_path = '\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_6MWT.xlsx';
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
        [Sig, TStr, Raw] = xlsread(MWT6_path, SubjectID)
        Raw([1],:) = []; % remove the header row from the dates
        TStr([1],:) = []; % remove the header row from the dates
        Raw(strcmp(Raw, 'DNA')) = {NaN}; % replace 'DNA' with nan
        %RawScore = cell2mat(Raw(:,3)); % get the raw scores from the Raw cells to make indexing easier
        %dimBBS = size(Raw);
        
        % Get the index of each group of dates for each session
        SN1_start=1;
        for k = 1:1:length(TStr(:,1))
            if datenum(TStr(SN1_start)) == datenum(TStr(k))
                SN1_end = SN1_start+k-1;
                SN2_start = SN1_end+1;
                SN1 = SN1_start:SN1_end;
                SN2 = [];
                SN3 = [];
                SN4 = [];
            elseif datenum(TStr(SN2_start)) == datenum(TStr(k))
                SN2_end = SN1_start+k-1;
                SN3_start = SN2_end+1;
                SN2 = [SN2_start:SN2_end];
                SN3 = [];
                SN4 = [];
            elseif datenum(TStr(SN3_start)) == datenum(TStr(k))
                SN3_end = SN1_start+k-1;
                SN4_start = SN3_end+1;
                SN3 = [SN3_start:SN3_end];
                SN4 = [];
            elseif datenum(TStr(SN4_start)) == datenum(TStr(k))
                SN4_end = SN1_start+k-1;
                SN4 = [SN4_start:SN4_end];
            end
            
        end
        SN = {SN1, SN2, SN3, SN4};
        
        
        % check to make sure the the dates are in the correct order
        % make the last session SN4
        
        if isempty(SN2) == 1 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(TStr(SN1_start))];
            order = 1;
        elseif isempty(SN2) == 0 && isempty(SN3) == 1 && isempty(SN4) == 1
            dates = [datenum(TStr(SN1_start)), datenum(TStr(SN2_start))];
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
            
            dates = [datenum(TStr(SN1_start)), datenum(TStr(SN2_start)), datenum(TStr(SN3_start))];
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
            dates = [datenum(TStr(SN1_start)), datenum(TStr(SN2_start)), datenum(TStr(SN3_start)), datenum(TStr(SN4_start))];
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
            
            for g = 1:1:length(Minute)
                if strcmp(Type_of_Subject_Group, 'CVA')
                    if isempty(SN{j}) == 1 || isempty(Sig) == 1 || ID(n) == 32 %CVA32 has no data
                    subject(g,:) = {SubjectID};
                    group(g,:) = Type_of_Subject(h);
                    session(g,:) = j;
                    activity(g,:) = Activity;
                    minute(g,:) = Minute(g);
                    stepcount(g,:) = nan;
                    distance_m(g,:) = nan;
                    else
                        if g > length(SN{j}(:)) || isempty(find(ismember(Sig(SN{j},2), Minute(g)))) == 1
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        session(g,:) = j;
                        activity(g,:) = Activity;
                        minute(g,:) = Minute(g);
                        stepcount(g,:) = nan;
                        distance_m(g,:) = nan;
                        
                        else
                        Minute_index = find(ismember(Sig(SN{j},2), Minute(g)));
                        
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        session(g,:) = j;
                        activity(g,:) = Activity;
                        minute(g,:) = Minute(g);
                        stepcount(g,:) = Sig(SN{j}(Minute_index),3);
                        distance_m(g,:) = Sig(SN{j}(Minute_index),5);
                        end
                    end
                elseif strcmp(Type_of_Subject_Group, 'CONTROLS')
                    if isempty(SN{j}) == 1 || isempty(Sig) == 1
                    subject(g,:) = {SubjectID};
                    group(g,:) = Type_of_Subject(h);
                    session(g,:) = j;
                    activity(g,:) = Activity;
                    minute(g,:) = Minute(g);
                    stepcount(g,:) = nan;
                    distance_m(g,:) = nan;
                    else
                        if g > length(SN{j}(:)) || isempty(find(ismember(Sig(SN{j},1), Minute(g)))) == 1
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        session(g,:) = j;
                        activity(g,:) = Activity;
                        minute(g,:) = Minute(g);
                        stepcount(g,:) = nan;
                        distance_m(g,:) = nan;
                        
                        else
                        Minute_index = find(ismember(Sig(SN{j},1), Minute(g)));
                        
                        subject(g,:) = {SubjectID};
                        group(g,:) = Type_of_Subject(h);
                        session(g,:) = j;
                        activity(g,:) = Activity;
                        minute(g,:) = Minute(g);
                        stepcount(g,:) = nan;
                        distance_m(g,:) = Sig(SN{j}(Minute_index),3);
                        end
                    end
                    
                
                end
            end
            
            % 6MWT Clinical Scores Matrix
            MWT6_Scores_temp = table(subject, group, session, activity, minute, stepcount, distance_m);
            MWT6_Scores = [MWT6_Scores; MWT6_Scores_temp]
            
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

writetable(MWT6_Scores,'MWT6_ClinicalScores.csv','Delimiter',',','QuoteStrings',true)






