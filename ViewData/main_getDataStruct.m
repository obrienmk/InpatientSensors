clc
close all
clear all
%**************************************************************************
% main_getDataStruct.m
%
% Create primary data structure of raw MC10 data
% All subjects, sessions, tests, locations, modalities.
%
% Loads all session dates for chosen subjects, clinical tests, and sensor
% locations.
%
% Structure is:
%   Subject.SessionDate.ClinicalTest.SensorLocation.TestItem.SensorModality
%
% FIX TestItem
% ADD LABEL from Excel
% ADD PT Session
%
% Written by Megan O'Brien
% version 09.18.2018
%**************************************************************************

Subject_type = 'cva'; % cva or controls

Subj2Load = 4;
Act2Load = 3:7;
Loc2Load = [5 6 14];

% Directories
dirMfile = pwd;
dirData = ['Z:\Inpatient Sensors -Stroke\MC10 Study\Data\biostamp_data\' Subject_type filesep];
dirDem = ['Z:\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\'];

% Activities
ActivitiesAll = {'Clinical - MAS','Clinical - MMT','Clinical - BBS','Clinical - 10MWT SSV','Clinical - 10MWT FV',...
    'Clinical - 6MWT','Clinical - TUG','Activity Recognition','Physical Therapy'}; 
ActivitiesAll_nickname = {'MAS','MMT','BBS','MWT10_SSV','MWT10_FV',...
    'MWT6','TUG','AR','PT'}; %9

% Locations
LocationsAll = {'bicep_left','bicep_right','biceps_femoris_left','biceps_femoris_right',...
    'distal_lateral_shank_left','distal_lateral_shank_right','gastrocnemius_left','gastrocnemius_right',...
    'medial_chest','posterior_forearm_left','posterior_forearm_right','rectus_femoris_left','rectus_femoris_right',...
    'sacrum','tibialis_anterior_left','tibialis_anterior_right'}; %16

% Subjects
SubjectsStruct = dir(dirData);
SubjectsAll = {SubjectsStruct.name};
SubjectsAll(1:2) = [];

% Demographics
Demographics=get_demographics(Subject_type,dirDem);

tic
for indSub = Subj2Load
    DataSensor = [];
    
    Subject = SubjectsAll{indSub};
    dirSubj = [dirData Subject filesep];
    
    fprintf('%s\n',Subject)
    %------------------------- DEMOGRAPHICS -----------------------
    
    
    %------------------------- OPEN LABELS / ANNOTATIONS -----------------------
    Labels=table2cell(readtable([dirSubj 'annotations.csv'],'ReadVariableNames',false));
    Labels=Labels(2:end,:);  %remove header
    
    DataTime = get_timestamps(Labels);

    %--------------------- GET TASKS AND DAYS OF EACH ACTIVITY -----------------
    for indAct = Act2Load
        Activity = ActivitiesAll{indAct};
        Activity_nickname = ActivitiesAll_nickname{indAct};
        
        fprintf('    %s\n',Activity)

        %[mode,tasks,sensors] = get_sensors({Activity});
        tasks = get_tasks({Activity});
        [taskInfo] = get_annotationind(tasks,{Activity},Labels);

        % Final tasks we actually have
        Activity_tasks = unique([taskInfo{:,2}]);
        nTasksTotal = length(Activity_tasks);

        % Annotation details
        Ann_idx_AllDays = taskInfo(:,1);
        Ann_act_AllDays = taskInfo(:,2);
        
        nTasks_AllDays = size(taskInfo,1);
        whichact_All = NaN*ones(1,nTasks_AllDays);
        for i=1:nTasks_AllDays
            % Get annotation index of desired data label (first index in annotation block)
            whichact_All(i) = Ann_idx_AllDays{i}(1);
        end

        LabelsAll = Labels(:,3);                               % names of all Labels
        LabelsAll_startTsMs = str2num(cell2mat(Labels(:,5)));  % starttime of all Labels
        LabelsAll_endTsMs = str2num(cell2mat(Labels(:,6)));    % stoptime of all Labels

        % Get date of each desired activity and store in annotation block
        whichdays = DataTime(whichact_All);
        for i=1:nTasks_AllDays
            taskInfo{i,3} = cellstr(char(whichdays(i)));
        end
        whichdays_unq = unique(cellstr(char(whichdays)));
        
        for indDay = 1:size(whichdays_unq,1)
            Day = ['d' regexprep(whichdays_unq{indDay},'[-]','')];
            
            % Separate taskInfo structure into tasks on just this day
            findDay = strcmp([taskInfo{:,3}]',whichdays_unq(indDay));
            dayInfo=taskInfo(findDay,:);
            nTasks_day = size(dayInfo,1);
                        
            % Number of tasks on this day
            Ann_idx = dayInfo(:,1);
            whichact = NaN*ones(1,nTasks_day);
            for i=1:nTasks_day
                % Get annotation index of desired data label (first index in annotation block)
                whichact(i) = Ann_idx{i}(1);
            end
            
            fprintf('        %s  ( /%i)\n',Day,nTasks_day)

            %----------------------- LOAD DATA FROM EACH SENSOR ------------------------
            % Data is stored as: Subj \ sensor location \ sensor ID \ test ID (timestamp) \ sensor mode
            LocsStruct = dir(dirSubj);
            LocsStruct(1:2) = [];  
            findAnn = strfind({LocsStruct.name},'annotations');
            LocsStruct(~cellfun(@isempty,findAnn)) = []; % removes annotation file(s) too
            LocsPresent = {LocsStruct.name};
                        
            indEcg = 0;  % increment for medial_chest col
            for indLoc = Loc2Load  % Sensor location on body
                Loc = LocationsAll{indLoc};
                
                dirLoc = [dirSubj Loc filesep];
                
                % If this location does not exist (was not recorded for this patient), move on to next
                if ~exist(dirLoc)
                    fprintf('Location not found\n')
                    continue
                end

                fprintf('            *%s\n',Loc)

                % Look for desired data in each sensor and recording for this subject and location
                flag_datafound = 0;  % flag for finding data (0 for not found, 1 for found)

                % Each sensor
                SensorsStruct = dir(dirLoc);
                SensorsStruct(1:2) = [];
                Sensors = {SensorsStruct.name};
                
                % Initialize data structure for this day and loc
                DataSensor = [];

                % Each recording from sensor
                for indSens = 1:length(Sensors)
                    Sensor = Sensors{indSens};
                    dirSensor = [dirLoc Sensor filesep];

                    RecStruct = dir(dirSensor);
                    RecStruct(1:2) = [];
                    RecIDs = {RecStruct.name};

                    for indRec = 1:length(RecIDs)  % check each recording within sensor
                        RecID = RecIDs{indRec};

                        % Does this recording date match with date of activity of interest?
                        RecID_char = char(RecID);

                        if ~strcmp(RecID_char(1:10),char(whichdays_unq(indDay)))  % dates do not match                        
                            continue  % if no data found, move on to next recording (indRec)
                        end

                        dirFile = [dirSensor RecID filesep];

                        ModesStruct = dir(dirFile);
                        ModesStruct(1:2) = [];
                        Modes = {ModesStruct.name};
                        
                        for indMode = 1:length(Modes)
                            Mode = Modes{indMode};
                            filename = [dirFile Mode];
                            
                            fid=fopen(filename,'rt');
                            SensorData=textscan(fid,'%f%f%f%f','Delimiter',',','HeaderLines',1);
                            fclose(fid);

                            if strcmp(Mode,'elec.csv')
                                Mode_nickname = 'Elec';
                                fs = 1000; datacol = 2;
                            elseif strcmp(Mode,'accel.csv')
                                Mode_nickname = 'Acc';
                                fs = 62.5; datacol = [2 3 4];
                            elseif strcmp(Mode,'gyro.csv')
                                Mode_nickname = 'Gyr';
                                fs=62.5; datacol = [2 3 4];
                            elseif strfind(Mode,'-errors')  % Skip error files
                                continue
                            end

                            tMs = cell2mat(SensorData(1)); % unpack
                            y = cell2mat(SensorData(datacol));
                            
                            fprintf('                 %s  ',Mode_nickname)
                            
                            for indTask = 1:nTasks_day
                                
                                [c1, ind1] = min(abs(tMs-LabelsAll_startTsMs( whichact(indTask) )));  %find closest data time to task start time
                                [c2, ind2] = min(abs(tMs-LabelsAll_endTsMs( whichact(indTask) )));  %find closest data time to task end time

                                % If task is not in this recording, if will find closest start and end time as the last index
                                % Move on to next loop
                                if ind1==ind2
                                    fprintf('[X]' )
                                    continue
                                else
                                    %indData(indPlot) = indData(indPlot) + 1;  % increment plots
                                    fprintf('%i ',indTask) % print to show found data
                                end
                                
                                if strcmp(Mode_nickname,'Elec') && strcmp(Loc,'medial_chest')
                                    indEcg = indEcg+1;
                                    Mode_nickname = [Mode_nickname '_' num2str(indEcg)];
                                end


                                indStart = ind1-0*1/fs;  % -2000 Add time before to see activation onset
                                indEnd = ind2+0*1/fs;    % +3000 Add time after to see activation cessation

                                DataSensor = y(indStart:indEnd,:);
                                dayInfo{indTask,4}.(Mode_nickname) = DataSensor;
                                
                                
                            end %indTask
                            fprintf('\n')
                        end  %indMode
                    end %indRec
                end %indSens

                DataStruct.(Subject).(Day).(Activity_nickname).(Loc)=dayInfo;
                if size(dayInfo,2)>3
                    dayInfo(:,4)=[];
                end

            end %indLoc
        end %indDay
    end %indAct
end %indSub
toc                
                
        
% Save final output
clearvars -except DataStruct dirMfile dirSave

dirSave = [dirMfile '\data\'];
save([dirSave,'DataStruct'])
fprintf('***\nSaved!\n%sDataStruct.mat\n***\n',dirSave)


%% % -------------------------------- PLOT ---------------------------------

FinalSubj = fieldnames(DataStruct); 
for s=1:numel(FinalSubj)
    Subject = FinalSubj{s};
    
    FinalDate = fieldnames(DataStruct.(Subject));   
    for d=1:numel(FinalDate)
        Date = FinalDate{d};
        
        FinalAct = fieldnames(DataStruct.(Subject).(Date));
        for a=1:numel(FinalAct)
            Act = FinalAct{a};
            
            figure(100*s + 10*d + a)
            set(gcf,'name',[Subject '-' Act '-' Date])
            
            FinalLoc = fieldnames(DataStruct.(Subject).(Date).(Act));
            for l=1:numel(FinalLoc)
                Loc = FinalLoc{l};
                
                % Choose which task to plot
                Task2Load = 1;
                TaskName = DataStruct.(Subject).(Date).(Act).(Loc){Task2Load,2};
                
                if isempty(TaskName)
                    TaskName = ['Trial ' num2str(Task2Load)];
                else
                    TaskName = TaskName{1};
                end
                
                % If there is a column with sensor data for this task, plot it
                % Otherwise, missing sensor data for this task
                if size(DataStruct.(Subject).(Date).(Act).(Loc),2)>=4
                    
                    FinalMod = fieldnames(DataStruct.(Subject).(Date).(Act).(Loc){Task2Load,4});
                    for m=1:numel(FinalMod)
                        Mod = FinalMod{m};

                        subplot(numel(FinalLoc),numel(FinalMod),2*(l-1)+m)
                        plot(DataStruct.(Subject).(Date).(Act).(Loc){Task2Load,4}.(Mod))
                        ylabel(Mod)
                        title([TaskName ' (' Loc ')'])
                    end %m
                end %if sensor data
            end %l
        end %a
    end %d
end %s

% % PT
% counter_PT = 0;
% nDays = length(fieldnames((DataStruct.CVA08)));
% for indDay=1:nDays
%     Day = days{indDay};
%     nLabels = size(DataStruct.CVA08.(Day).PT.sacrum,1);
%     for indLabel=1:nLabels
%         counter_PT = counter_PT+1;
%         subplot(3,2,counter_PT)
%         plot(DataStruct.CVA08.(Day).PT.sacrum{indLabel,4}.Acc)
%         title(strcat(DataStruct.CVA08.(Day).PT.sacrum{indLabel,2},': ',Day))
%     end
% end
% fprintf('# of PT sessions: %i\n',counter_PT)
