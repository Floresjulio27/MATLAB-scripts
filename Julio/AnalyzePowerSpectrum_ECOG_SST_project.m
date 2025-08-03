function [Results_PowerSpectrum_ECOG] = AnalyzePowerSpectrum_ECOG_SST_project(animalID,group,setName,rootFolder,delim,Results_PowerSpectrum_ECOG,days)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted and modified by Julio Flores-Cuadra
% Department of Biology
%________________________________________________________________________________________________________________________
%
%  Purpose: Analyze the spectral power of ECOG signal
%________________________________________________________________________________________________________________________

%% Consider the animals that have good ECOG for plotting
%% function parameters
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% Set folder
dataLocation = [rootFolder delim 'Data' delim group delim setName delim days delim animalID delim 'CombinedImaging']; % here is the key how to organize it
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
if isempty (procDataFileStruct)
    return 
end 
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat strut
SleepDataFileStruct = dir('*_SleepData.mat');
SleepDataFile = {SleepDataFileStruct.name}';
SleepDataFileID = char(SleepDataFile);
load(SleepDataFileID,'-mat')
% scoring results
% find and load manual baseline event information
scoringResultsFileStruct = dir('*Manual_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% filter characteristics & resting criteria
samplingRate = RestData.CBV.P_ACh.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {10};%5

% go through each valid data type for behavior-based power spectrum analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze power spectra during periods of rest
    % pull data from RestData.mat structure
    if strcmp(dataType,'deltaBandPower') == true
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        restData = RestData.cortical_LH.(dataType).data(combRestLogical,:);
    elseif strcmp(dataType,'thetaBandPower') == true
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        restData = RestData.cortical_LH.(dataType).data(combRestLogical,:);
    elseif strcmp(dataType,'alphaBandPower') == true
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        restData = RestData.cortical_LH.(dataType).data(combRestLogical,:);
    elseif strcmp(dataType,'betaBandPower') == true
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        restData = RestData.cortical_LH.(dataType).data(combRestLogical,:);
    elseif strcmp(dataType,'gammaBandPower') == true
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        restData = RestData.cortical_LH.(dataType).data(combRestLogical,:);
    end

    % keep only the data that occurs within the manually-approved awake regions
    [finalRestData,~,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    clear procRestData
    if isempty(finalRestData) == false
        % detrend and truncate data to minimum length to match events
        for bb = 1:length(finalRestData)
            if length(finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{bb,1});
                restPad = (ones(1,restChunkSampleDiff))*finalRestData{bb,1}(end);
                procRestData{bb,1} = horzcat(finalRestData{bb,1},restPad); %#ok<*AGROW>
                procRestData{bb,1} = detrend(procRestData{bb,1},'constant');
            else
                procRestData{bb,1} = detrend(finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end

        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for cc = 1:length(procRestData)
            if sum(isnan(procRestData{cc,1})) == 0
                restDataMat(:,zz) = procRestData{cc,1};
                zz = zz + 1;
            end
        end

        % parameters for mtspectrumc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];

        % calculate the power spectra of the desired signals
        [rest_S,rest_f,rest_sErr] = mtspectrumc(restDataMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).S = rest_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).f = rest_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).sErr = rest_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Rest.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of alert
    zz = 1;
    clear awakeData procAwakeData
    awakeData = [];
    for cc = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(cc,:);
        [~,~,awakeDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
        scoringLabels = [];
        for dd = 1:length(ScoringResults.fileIDs)
            if strcmp(awakeDataFileID,ScoringResults.fileIDs{dd,1}) == true
                scoringLabels = ScoringResults.labels{dd,1};
            end
        end
        % cut the data into three 10 minutes chunk
        ScoreSize = 10*12;% 10 minutes X 12 5 sec periods per minutes

        load(procDataFileID,'-mat')
        DataSize = 10*60*ProcData.notes.dsFs; % 15 minutes X 60 secs X 30 Hz;
        % only run on files with good pupil measurement
        if strcmp(ProcData.data.Pupil,'Diameter') == true
            try
                puffs = ProcData.data.stimulations.LPadSol;
            catch
                puffs = ProcData.data.solenoids.LPadSol;
            end
            % don't include trials with stimulation
            if isempty(puffs) == true
                if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                    for SSL =  1:1:5
                        ScoreLabels = scoringLabels(((SSL-1)*ScoreSize)+1:(SSL*ScoreSize));
                        % check labels to match arousal state
                        if sum(strcmp(ScoreLabels,'Not Sleep')) > 108 % 90% of the time awake
                            % pull data based on data type
                            if strcmp(dataType,'deltaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.deltaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'thetaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.thetaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'alphaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.alphaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'betaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.betaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'gammaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.gammaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            end
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
    end
    % calculate the power spectra
    if isempty(awakeData) == false
        % detrend data
        for bb = 1:length(awakeData)
            procAwakeData{bb,1} = detrend(awakeData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        awakeDataMat = zeros(length(procAwakeData{1,1}),length(procAwakeData));
        for cc = 1:length(procAwakeData)
            awakeDataMat(:,cc) = procAwakeData{cc,1}(1:length(procAwakeData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [5,9];   % Tapers [n, 2n - 1] [5,9]; Kevin never changed the tapers but Shak does. To reduce noise?
        [awake_S,awake_f,awake_sErr] = mtspectrumc(awakeDataMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).S = awake_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).f = awake_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).sErr = awake_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Awake.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of Asleep
    zz = 1;
    clear asleepData procAsleepData

    asleepData = [];
    for cc = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(cc,:);
        [~,~,asleepDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
        scoringLabels = [];
        for dd = 1:length(ScoringResults.fileIDs)
            if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                scoringLabels = ScoringResults.labels{dd,1};
            end
        end
        % cut the data into three 10 minutes chunk
        ScoreSize = 10*12;% 10 minutes X 12 5 sec periods per minutes

        load(procDataFileID,'-mat')
        DataSize = 10*60*ProcData.notes.dsFs; % 10 minutes X 60 secs X 30 Hz;
        % only run on files with good pupil measurement
        if strcmp(ProcData.data.Pupil,'Diameter') == true
            try
                puffs = ProcData.data.stimulations.LPadSol;
            catch
                puffs = ProcData.data.solenoids.LPadSol;
            end
            % don't include trials with stimulation
            if isempty(puffs) == true
                if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                    for SSL =  1:1:5
                        ScoreLabels = scoringLabels(((SSL-1)*ScoreSize)+1:(SSL*ScoreSize));
                        % check labels to match arousal state
                        if sum(strcmp(ScoreLabels,'Not Sleep')) < 36 % 70% of the time asleep
                            % pull data based on data type
                            if strcmp(dataType,'deltaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.deltaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'thetaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.thetaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'alphaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.alphaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'betaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.betaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                            elseif strcmp(dataType,'gammaBandPower') == true
                                awakeData{zz,1} = ProcData.data.cortical_LH.gammaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSizeResults_PowerSpectrum));
                            end
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
    end

    % calculate the power spectra
    if isempty(asleepData) == false
        % detrend data
        for bb = 1:length(asleepData)
            procAsleepData{bb,1} = detrend(asleepData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        asleepDataMat = zeros(length(procAsleepData{1,1}),length(procAsleepData));
        for cc = 1:length(procAsleepData)
            asleepDataMat(:,cc) = procAsleepData{cc,1}(1:length(procAsleepData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
        [asleep_S,asleep_f,asleep_sErr] = mtspectrumc(asleepDataMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).S = asleep_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).f = asleep_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).sErr = asleep_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.Asleep.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of all data
    zz = 1;
    clear allData procAllData
    allData = [];

    for cc = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(cc,:);
        [~,~,~] = GetFileInfo_JNeurosci2022(procDataFileID);
        load(procDataFileID,'-mat')

        DataSize = 10*60*ProcData.notes.dsFs; % 15 minutes X 60 secs X 30 Hz;
        % only run on files with good pupil measurement
        if strcmp(ProcData.data.Pupil,'Diameter') == true
            try
                puffs = ProcData.data.stimulations.LPadSol;
            catch
                puffs = ProcData.data.solenoids.LPadSol;
            end
            % don't include trials with stimulation
            if isempty(puffs) == true
                if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                    % cut the data into three 10 minutes chunk
                    for SSL =  1:1:5
                        % pull data based on data type
                        if strcmp(dataType,'deltaBandPower') == true
                            awakeData{zz,1} = ProcData.data.cortical_LH.deltaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                        elseif strcmp(dataType,'thetaBandPower') == true
                            awakeData{zz,1} = ProcData.data.cortical_LH.thetaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                        elseif strcmp(dataType,'alphaBandPower') == true
                            awakeData{zz,1} = ProcData.data.cortical_LH.alphaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                        elseif strcmp(dataType,'betaBandPower') == true
                            awakeData{zz,1} = ProcData.data.cortical_LH.betaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSize));
                        elseif strcmp(dataType,'gammaBandPower') == true
                            awakeData{zz,1} = ProcData.data.cortical_LH.gammaBandPower(((SSL-1)*DataSize)+1:(SSL*DataSizeResults_PowerSpectrum));
                        end
                        zz = zz + 1;
                    end
                end
            end
        end
    end

    % calculate the power spectra
    if isempty(allData) == false
        % detrend data
        for bb = 1:length(allData)
            procAllData{bb,1} = detrend(allData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        allDataMat = zeros(length(procAllData{1,1}),length(procAllData));
        for cc = 1:length(procAllData)
            allDataMat(:,cc) = procAllData{cc,1}(1:length(procAllData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
        [all_S,all_f,all_sErr] = mtspectrumc(allDataMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).S = all_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).f = all_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).sErr = all_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.All.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of NREM
    % pull data from SleepData.mat structure
    if strcmp(dataType,'deltaBandPower') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.deltaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'thetaBandPower') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.thetaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'alphaBandPower') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.alphaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'betaBandPower') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.betaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'gammaBandPower') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.gammaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end

    if isempty(nremData) == false
        % detrend and truncate data to minimum length to match events
        for dd = 1:length(nremData)
            nremData{dd,1} = detrend(nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for ee = 1:length(nremData)
            if sum(isnan(nremData{ee,1})) == 0
                nremMat(:,zz) = nremData{ee,1};
                zz = zz + 1;
            end
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        [nrem_S,nrem_f,nrem_sErr] = mtspectrumc(nremMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).S = nrem_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).f = nrem_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).sErr = nrem_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.NREM.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of REM
    % pull data from SleepData.mat structure
     if strcmp(SleepData.(modelType),'REM') == true
        if strcmp(dataType,'deltaBandPower') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.deltaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'thetaBandPower') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.thetaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'alphaBandPower') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.alphaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'betaBandPower') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.betaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'gammaBandPower') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.gammaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        end
    else
        remData = [];
    end
    if isempty(remData) == false
        % detrend and truncate data to minimum length to match events
        for dd = 1:length(remData)
            remData{dd,1} = detrend(remData{dd,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for ee = 1:length(remData)
            if sum(isnan(remData{ee,1})) == 0
                remMat(:,zz) = remData{ee,1};
                zz = zz + 1;
            end
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        [rem_S,rem_f,rem_sErr] = mtspectrumc(remMat,params);
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).S = rem_S;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).f = rem_f;
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).sErr = rem_sErr;
    else
        % save results
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).S = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).f = [];
        Results_PowerSpectrum_ECOG.(group).(days).(animalID).PowerSpectrum.REM.(dataType).sErr = [];
    end
end
%% save data
cd([rootFolder delim 'Results_SST_project'])
save('Results_PowerSpectrum_ECOG.mat','Results_PowerSpectrum_ECOG','-v7.3')
cd([rootFolder delim 'Data'])

end
