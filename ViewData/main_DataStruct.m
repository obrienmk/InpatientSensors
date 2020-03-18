clc
clear all
close all
% *************************************************************************
% main_DataStruct.m
%
% Load and clean all sensor data from selected subjects for a desired
% activity (e.g. 10MWT or TUG). Save data structure to folder in current directory.
%
%   ACTIVITY NAMES
% -------------------
%  - 'MAS     - 'MWT10_SSV'    - 'TUG'           - AR (controls only)
%  - 'MMT'    - 'MWT10_FV'     - 'Resting_ECG'
%  - 'BBS'    - 'MWT6'         - 'PT'
% -------------------
% (names taken from folders in Z:\Inpatient Sensors -Stroke\MC10 Study\Data analysis\Segmented_Data\
%
% Run after segmenting data using Lars' code.
%
% Megan O'Brien
% Version date 02-19-2020
% *************************************************************************

%% ------------------------------- INPUT ----------------------------------
% Select data to load

Subject_type = 'cva';    % cva or controls
Activity = 'MWT10_SSV';  % Select one from list of nicknames below
Subj2Load = [1:4];   % Subject numbers to load (1:55)

saveon = 0;              % flag to save data structure

%% ---------------------------- INITIALIZE --------------------------------
% Directories
dirMfile = pwd;
dirData = ['Z:\Inpatient Sensors -Stroke\MC10 Study\Data analysis\1_Segmented_Data\' Subject_type filesep Activity filesep];  % Segmented sensor data
dirMeta = 'Z:\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\';              % Metadata for participants
addpath([dirMfile filesep 'clean'])

% Initialize structure
DataStruct=struct('SubjID', '', 'Timepoint','','Activity', '', 'Trial', '', 'Location', '', 'SensorData', []);

%% ---------------------------- MAIN LOOP ---------------------------------

% SUBJECTS
struct_counter = 1;
for indSub = 1:length(Subj2Load)
    Subject_num = Subj2Load(indSub);
    Subject = get_subjectId(Subject_num,Subject_type);
    
    fprintf('* %s \n',Subject)
    
    % 1. --------------------- LOAD SEGMENTED DATA ------------------------
    TempStruct = load([dirData Subject '.mat']);
    
    % DATES
    Dates = fieldnames(TempStruct.(Subject));
    
    % Sort by date
    DateTimes = datetime(erase(Dates,'ymd'),'InputFormat','yyyy_MM_dd');
    [DateTimes_sorted, indsort] = sort(DateTimes);
    Dates = Dates(indsort);
    
%     Dates(:,1) = fieldnames(TempStruct.(Subject));
%     % Identify Admission & Discharge dates (first and last dates)
%     DateTimes = datetime(erase(Dates(:,1),'ymd'),'InputFormat','yyyy_MM_dd');
%     [DateTimes_sorted, indsort] = sort(DateTimes);
%     Dates(:,1) = Dates(indsort,1);
%     numDates = numel(Dates(:,1));
%     
%     Dates{1,2} = 'Adm';
%     if numDates>1
%         Dates{end,2} = 'Dis';
%         if numDates>2
%             for indMid=1:numDates-2
%                 Dates{indMid+1,2} = ['Mid_' num2str(indMid)]
%             end
%         end
%     else
%         fprintf('Only one timepoint!')
%     end
    
    for indDate = 1:numel(Dates)
        Date = Dates{indDate};
%         Date = Dates{indDate,1};
        DateTime = datetime(erase(Date,'ymd'),'InputFormat','yyyy_MM_dd');

        
        % ACTIVITIES
        Acts = fieldnames(TempStruct.(Subject).(Date));
        for indAct = 1:numel(Acts)
            Act = Acts{indAct};

            % LOCATIONS
            Locs = fieldnames(TempStruct.(Subject).(Date).(Act));
            for indLoc = 1:numel(Locs)
                Loc = Locs{indLoc};
                
                % TRIALS / TASKS
                Trials = fieldnames(TempStruct.(Subject).(Date).(Act).(Loc));
                for indTrial = 1:numel(Trials)
                    Trial = Trials{indTrial};
                    
                    % Will concatenate all ecg data from this trial into single Ecg struct (time, Sensor 1, Sensor 2, Sensor 3), unordered
                    if strcmp('medial_chest',Loc)
                        ecg_counter = 1;
                    end
                    
                    % SENSOR MODALITIES
                    Mods = fieldnames(TempStruct.(Subject).(Date).(Act).(Loc).(Trial));
                    AccData = []; GyrData = []; ElecData = [];
                    for indMod = 1:numel(Mods)
                        Mod = Mods{indMod};
                        
                        % 2. ----------------------------- CLEAN DATA -----------------------------
                        if contains(Mod,'acc')
                            Fs = 31.25; % Sampling freq (Hz)
                            ncol = 3;   % Data col (X,Y,Z)
                        elseif contains(Mod,'gyr')
                            Fs = 31.25; % Sampling freq (Hz)
                            ncol = 3;   % Data col (X,Y,Z)
                        elseif contains(Mod,'elec')
                            Fs = 1000;  % Sampling freq (Hz)
                            ncol = 1;   % Data col (V)
                        end
                        
                        SensorData.time = TempStruct.(Subject).(Date).(Act).(Loc).(Trial).(Mod).data(:,1)/1000; %time in sec
                        SensorData.time = SensorData.time-SensorData.time(1);
                        SensorData.dataRaw = TempStruct.(Subject).(Date).(Act).(Loc).(Trial).(Mod).data(:,2:(ncol+1));
                        
                        % Remove duplicate timestamps
                        SensorData = removeduplicateTime( SensorData );
                        
                        % Resampling
                        t=SensorData.time(1):1/Fs:max(SensorData.time);
                        SensorData.dataClean=spline((SensorData.time.'),SensorData.dataRaw(:,1:ncol).',t);  % spline interpolation on all columns of data
                        SensorData.dataClean=SensorData.dataClean';
                        SensorData.time = t';
                        
                        if contains(Mod,'acc')
                            AccData = [SensorData.time SensorData.dataClean];
                        elseif contains(Mod,'gyr')
                            GyrData = [SensorData.time SensorData.dataClean];
                        elseif contains(Mod,'elec')
                            if strcmp('medial_chest',Loc)
                                if ecg_counter == 1
                                    ElecData = [SensorData.time SensorData.dataClean];
                                else
                                    ElecData = [ElecData SensorData.dataClean];
                                end
                                ecg_counter = ecg_counter+1;
                            else
                                ElecData = [SensorData.time SensorData.dataClean];    
                            end
                        end
                        
                    end %Mods
                    
                    SensorDataStruct=struct('Acc',AccData,'Gyr',GyrData,'Elec',ElecData);
                    DataStruct(struct_counter)=struct('SubjID', Subject, 'Timepoint', erase(Date,'ymd'),'Activity', Act, 'Trial', Trial, 'Location', Loc, 'SensorData', SensorDataStruct);
                    struct_counter = struct_counter+1;
                end %Trials
            end %Locs
        end %Acts
    end %Dates
 end %Subj

%% ------------------------------- OUTPUT ---------------------------------
% Save structure
if saveon
    dirSave = [dirMfile filesep 'Saved' filesep];
    if ~exist(dirSave, 'dir'); mkdir(dirSave); end
    
    filenameSave = ['DataStruct' Subject_type '_' Activity '.mat'];
    
    save([dirSave filenameSave],'DataStruct')
end

% % Plot (final loaded subject)
% figure; hold on;
% plot(SensorData.time,SensorData.dataClean);
% xlabel('Time (s)');