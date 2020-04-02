function main_Data_Extract(Subject_type,Activity,Subj2Load)
    % *************************************************************************
    % main_Data_Extract.m
    %
    % Load and clean all sensor data from selected subjects for a desired
    % activity (e.g. 10MWT or TUG). Save data structure to a specified folder.
    %
    % Run after segmenting data using Lars' code.
    %
    % Subjec_type: 'cva' or 'controls'
    %
    %   ACTIVITY NAMES
    % -------------------
    %  - 'MAS     - 'MWT10_SSV'    - 'TUG'           - AR (controls only)
    %  - 'MMT'    - 'MWT10_FV'     - 'Resting_ECG'
    %  - 'BBS'    - 'MWT6'         - 'PT'
    % -------------------
    % (names taken from folders in Z:\Inpatient Sensors -Stroke\MC10 Study\Data analysis\Segmented_Data\
    %
%   % Subj2Load: Subject numbers to load (1:55)
    %
    %   OUTPUT: structure with
    % -------------------
    %   Subject            % Session_trials (cell of trials from each Activity)
    %   Admission_date
    %   Final_date
    % -------------------
    %
    % Version date 04-02-2020
    % *************************************************************************
    
    %% ---------------------------- INITIALIZE --------------------------------
    % Directories
    if ~exist('RTOdir.mat','file') % Have user select path to RTO drive and save it for later use
        RTOdir = uigetdir(matlabroot,'Choose Local Path to RTO Drive'); 
        save('RTOdir.mat','RTOdir')
    else
        load('RTOdir.mat') % Load stored loca path to drive
    end

    dirMfile = pwd;
    dirData = fullfile(RTOdir,'Inpatient Sensors -Stroke','MC10 Study','Data analysis','1_Segmented_Data',Subject_type,Activity,filesep);  % Segmented sensor data
    dirMeta = fullfile(RTOdir,'Inpatient Sensors -Stroke','MC10 Study','Outcome Measures',filesep);              % Metadata for participants

    addpath([dirMfile filesep 'clean'])
    if ~exist('Data_out','dir')
        mkdir('Data_out')
    end

    % Initialize structure
    DataStruct=struct('SubjID', '', 'Timepoint','','Activity', '', 'Trial', '', 'Location', '', 'SensorData', []);
    NoS = 1;    % Number of Sessions (Date)
    %% ---------------------------- MAIN LOOP ---------------------------------

    % SUBJECTS
    struct_counter = 1;
    for indSub = 1:length(Subj2Load)

        % Data out folder, filename (Sung)
        file_out = ['C:\Users\mobrien\Desktop\Data_out\' upper(Subject_type) '_' Activity '_ID' sprintf('%02d',Subj2Load(indSub)) '.mat']


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

        % Subject & Date info (Sung)
        data.Subject = Subject;
        data.Admission_date = DateTimes_sorted(1);
        data.Final_date = DateTimes_sorted(end);

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


                if isfield(TempStruct.(Subject).(Date).(Act),'sacrum') == 0 || isempty(fieldnames(TempStruct.(Subject).(Date).(Act).sacrum)) == 1 
                    % skip
                else
                    % LOCATIONS
                    Locs = fieldnames(TempStruct.(Subject).(Date).(Act));
                    for indLoc = 1:numel(Locs)
                        Loc = Locs{indLoc};

                        % TRIALS / TASKS
                        Trials = fieldnames(TempStruct.(Subject).(Date).(Act).(Loc));
                        data.Session_trials(NoS) = {Trials};
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
                                if contains(Mod,'error')
                                    % skip
                                else
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
                                                if length(SensorData.dataClean) > length(ElecData(:,1))
                                                    ElecData = [ElecData SensorData.dataClean(1:length(ElecData(:,end)))];
                                                elseif length(SensorData.dataClean) < length(ElecData(:,1))
                                                    ElecData = [ElecData(1:length(SensorData.dataClean),:) SensorData.dataClean];
                                                else
                                                    ElecData = [ElecData SensorData.dataClean];
                                                end

                                            end
                                            ecg_counter = ecg_counter+1;
                                        else
                                            ElecData = [SensorData.time SensorData.dataClean];    
                                        end
                                    end
                                end

    %                             if contains(Mod,'acc') 
    %                                 sig_acc = TempStruct.(Subject).(Date).(Act).(Loc).(Trial).(Mod).data(:,2:4)
    %                             elseif contains(Mod,'elec')
    %                                 sig_EMG = TempStruct.(Subject).(Date).(Act).(Loc).(Trial).(Mod).data(:,2)
    %                             end

                            end %Mods

                            SensorDataStruct=struct('Acc',AccData,'Gyr',GyrData,'Elec',ElecData);
                            DataStruct(struct_counter)=struct('SubjID', Subject, 'Timepoint', erase(Date,'ymd'),'Activity', Act, 'Trial', Trial, 'Location', Loc, 'SensorData', SensorDataStruct);
                            struct_counter = struct_counter+1;

                            % Add here to save data (Sung)
                            data.Session_date(NoS) = DateTime;

                            % Location
                            % Motion Ssensors
                            if strcmp('sacrum',Loc)
                                data.Session{NoS}.Motion.Time{indTrial} = SensorDataStruct.Acc(:,1);

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.SC.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Gyr) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.SC.Gyr{indTrial} = SensorDataStruct.Gyr(:,2:4);
                                end

                            elseif strcmp('distal_lateral_shank_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.DLS_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Gyr) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.DLS_L.Gyr{indTrial} = SensorDataStruct.Gyr(:,2:4);
                                end

                            elseif strcmp('distal_lateral_shank_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.DLS_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Gyr) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.DLS_R.Gyr{indTrial} = SensorDataStruct.Gyr(:,2:4);
                                end

                            elseif strcmp('posterior_forearm_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.FA_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Gyr) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.FA_L.Gyr{indTrial} = SensorDataStruct.Gyr(:,2:4);
                                end

                            elseif strcmp('posterior_forearm_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.FA_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Gyr) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Motion.FA_R.Gyr{indTrial} = SensorDataStruct.Gyr(:,2:4);
                                end

                            % Physiological (EMG, ECG) Sensors
                            elseif strcmp('medial_chest',Loc)
                                data.Session{NoS}.Physio.Time_Elec{indTrial} = SensorDataStruct.Elec(:,1);
                                if isempty(SensorDataStruct.Acc) == 1
                                    data.Session{NoS}.Physio.Time_Acc{indTrial} = data.Session{NoS}.Motion.Time{indTrial};
                                else
                                    data.Session{NoS}.Physio.Time_Acc{indTrial} = SensorDataStruct.Acc(:,1);
                                end

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.MC.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.MC.ECG{indTrial} = SensorDataStruct.Elec(:,2:end);
                                end

                            elseif strcmp('rectus_femoris_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.RF_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.RF_L.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('rectus_femoris_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.RF_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.RF_R.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('tibialis_anterior_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.TA_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.TA_L.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('tibialis_anterior_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.TA_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.TA_R.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('biceps_femoris_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BF_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BF_L.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('biceps_femoris_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BF_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BF_R.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('gastrocnemius_left',Loc)  

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.GA_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.GA_L.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('gastrocnemius_right',Loc) 

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.GA_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.GA_R.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('bicep_left',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BC_L.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BC_L.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            elseif strcmp('bicep_right',Loc)

                                if isempty(SensorDataStruct.Acc) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BC_R.Acc{indTrial} = SensorDataStruct.Acc(:,2:4);
                                end
                                if isempty(SensorDataStruct.Elec) == 1
                                    % Skip
                                else
                                    data.Session{NoS}.Physio.BC_R.EMG{indTrial} = SensorDataStruct.Elec(:,2);
                                end

                            end  

                        end %Trials
                    end %Locs
                    NoS = NoS + 1;  % Remove empty session
                end % is empty?
            end %Acts  
        end %Dates 


        save(file_out,'data')
        NoS = 1;
        clear data 


    end %Subj






    % figure
    % plot(sig_acc)
    % hold on
    % plot(data.Session{2}.Physio.RF_R.Acc{1},'r--')
    %  
    % 
    % figure
    % plot(sig_EMG)
    % hold on
    % plot(data.Session{2}.Physio.RF_R.EMG{1},'r--')





