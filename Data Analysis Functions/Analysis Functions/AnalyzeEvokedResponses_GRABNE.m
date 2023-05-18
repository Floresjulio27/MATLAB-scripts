function [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
% extInd = strfind(group,delim);
% setName = group(extInd + 1:end);
% if strcmp(setName,'IOS Set A') == true
%     dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
% elseif strcmp(setName,'IOS Set B') == true
%     dataLocation = [rootFolder delim group delim animalID delim 'Combined Imaging'];
% end
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end
% dataTypes = {'LH','RH'};
%% only run analysis for valid animal IDs
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% criteria for whisking
WhiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaA.Comparison = {'gt','lt','gt'};
WhiskCriteriaA.Value = {0.5,2,5};
WhiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaB.Comparison = {'gt','lt','gt'};
WhiskCriteriaB.Value = {2,5,5};
WhiskCriteriaC.Fieldname = {'duration','puffDistance'};
WhiskCriteriaC.Comparison = {'gt','gt'};
WhiskCriteriaC.Value = {5,5};
WhiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% criteria for stimulation
StimCriteriaA.Value = {'LPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'RPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
StimCriteriaC.Value = {'AudSol'};
StimCriteriaC.Fieldname = {'solenoidName'};
StimCriteriaC.Comparison = {'equal'};
stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
%% analyze whisking-evoked responses
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.Rhodamine.Z_NE.whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.Rhodamine.Z_NE.whisk.trialDuration_sec;
    timeVector = (0:(EventData.Rhodamine.Z_NE.whisk.epoch.duration*samplingRate))/samplingRate - EventData.Rhodamine.Z_NE.whisk.epoch.offset;
    offset = EventData.Rhodamine.Z_NE.whisk.epoch.offset;
    for bb = 1:length(WhiskCriteriaNames)
        whiskCriteriaName = WhiskCriteriaNames{1,bb};
        if strcmp(whiskCriteriaName,'ShortWhisks') == true
            WhiskCriteria = WhiskCriteriaA;
        elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
            WhiskCriteria = WhiskCriteriaB;
        elseif strcmp(whiskCriteriaName,'LongWhisks') == true
            WhiskCriteria = WhiskCriteriaC;
        end
        % pull data from EventData.mat structure
        [whiskLogical] = FilterEvents_IOS(EventData.Rhodamine.Z_NE.whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [AchWhiskRhodamineData] = EventData.Rhodamine.Z_Ach.whisk.NormData(combWhiskLogical,:);
        [NEWhiskRhodamineData] = EventData.Rhodamine.Z_NE.whisk.NormData(combWhiskLogical,:);
        [AchWhiskGFPData] = EventData.GFP.Z_Ach.whisk.NormData(combWhiskLogical,:);
        [NEWhiskGFPData] = EventData.GFP.Z_NE.whisk.NormData(combWhiskLogical,:);

        [allWhiskCorticalMUAData] = EventData.cortical_LH.corticalPower.whisk.NormData(combWhiskLogical,:);
%         [allWhiskHippocampalMUAData] = EventData.hippocampus.corticalPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskCorticalGamData] = EventData.cortical_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
%         [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskFileIDs] = EventData.Rhodamine.Z_NE.whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.Rhodamine.Z_NE.whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.Rhodamine.Z_NE.whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [AchfinalWhiskRhodamineData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS(AchWhiskRhodamineData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [NEfinalWhiskRhodamineData,~,~,~] = RemoveInvalidData_IOS(NEWhiskRhodamineData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
               
        [AchfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(AchWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [NEfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(NEWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);

        [finalWhiskCorticalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
%         [finalWhiskHippocampalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskCorticalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
%         [finalWhiskHippocampalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskGFPData procWhiskRhodamineData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
        dd = 1;
        for cc = 1:size(NEfinalWhiskRhodamineData,1)
            whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 0;%2;
            whiskEndTime = whiskStartTime + 15;
            finalWhiskFileID = finalWhiskFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                NEwhiskRhodaminearray = NEfinalWhiskRhodamineData(cc,:);
                AchwhiskRhodaminearray = AchfinalWhiskRhodamineData(cc,:);

                NEwhiskGFParray = NEfinalWhiskGFPData(cc,:);
                AchwhiskGFParray = AchfinalWhiskGFPData(cc,:);

                whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
%                 whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
%                 whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                
                AchfiltWhiskRhodaminearray = sgolayfilt(AchwhiskRhodaminearray,3,9) - mean(AchwhiskRhodaminearray(1:(offset*samplingRate)));
                NEfiltWhiskRhodaminearray = sgolayfilt(NEwhiskRhodaminearray,3,9) - mean(NEwhiskRhodaminearray(1:(offset*samplingRate)));
                AchfiltWhiskGFParray = sgolayfilt(AchwhiskGFParray,3,9) - mean(AchwhiskGFParray(1:(offset*samplingRate)));
                NEfiltWhiskGFParray = sgolayfilt(NEwhiskGFParray,3,9) - mean(NEwhiskGFParray(1:(offset*samplingRate)));

                AchprocWhiskRhodamineData(dd,:) = AchfiltWhiskRhodaminearray;
                NEprocWhiskRhodamineData(dd,:) = NEfiltWhiskRhodaminearray;
                AchprocWhiskGFPData(dd,:) = AchfiltWhiskGFParray;
                NEprocWhiskGFPData(dd,:) = NEfiltWhiskGFParray;

                procWhiskCorticalMUAData(dd,:) = whiskCorticalMUAarray - mean(whiskCorticalMUAarray(1:(offset*samplingRate)));
%                 procWhiskHippocampalMUAData(dd,:) = whiskHippocampalMUAarray - mean(whiskHippocampalMUAarray(1:(offset*samplingRate)));
                procWhiskCorticalGamData(dd,:) = whiskCorticalGamArray - mean(whiskCorticalGamArray(1:(offset*samplingRate)));
%                 procWhiskHippocampalGamData(dd,:) = whiskHippocampalGamArray - mean(whiskHippocampalGamArray(1:(offset*samplingRate)));
                finalWhiskStartTimes(dd,1) = whiskStartTime;
                finalWhiskEndTimes(dd,1) = whiskEndTime;
                finalWhiskFiles{dd,1} = finalWhiskFileID;
                dd = dd + 1;
            end
        end
        AchmeanWhiskRhodamineData = mean(AchprocWhiskRhodamineData,1);%*100;
        AchstdWhiskRhodamineData = std(AchprocWhiskRhodamineData,0,1);%*100;
        NEmeanWhiskRhodamineData = mean(NEprocWhiskRhodamineData,1);%*100;
        NEstdWhiskRhodamineData = std(NEprocWhiskRhodamineData,0,1);%*100;

        AchmeanWhiskGFPData = mean(AchprocWhiskGFPData,1);%*100;
        AchstdWhiskGFPData = std(AchprocWhiskGFPData,0,1);%*100;
        NEmeanWhiskGFPData = mean(NEprocWhiskGFPData,1);%*100;
        NEstdWhiskGFPData = std(NEprocWhiskGFPData,0,1);%*100;

        meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
        stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
%         meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
%         stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
        meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
        stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
%         meanWhiskHippocampalGamData = mean(procWhiskHippocampalGamData,1)*100;
%         stdWhiskHippocampalGamData = std(procWhiskHippocampalGamData,0,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
%         whiskHippocampalZhold = [];
        for ee = 1:length(finalWhiskFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = 'cortical_LH';
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
%                     whiskHippocampalS_Data = AllSpecData.hippocampus.normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
%             whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskCorticalZhold = cat(3,whiskCorticalZhold,whiskCorticalS_Vals);
%             whiskHippocampalZhold = cat(3,whiskHippocampalZhold,whiskHippocampalS_Vals);
        end
        % cortical mean-subtract by first 2 seconds prior to stimulus
        meanWhiskCortS = mean(whiskCorticalZhold,3);
        baseWhiskCortS_Vals = mean(meanWhiskCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixWhiskCortS_Vals = baseWhiskCortS_Vals.*ones(size(meanWhiskCortS));
        msStimWhiskS_Vals = (meanWhiskCortS - baseMatrixWhiskCortS_Vals);
        % hippocampal mean-subtract by first 2 seconds prior to stimulus
%         meanWhiskHipS = mean(whiskHippocampalZhold,3);
%         baseWhiskHipS_Vals = mean(meanWhiskHipS(:,1:1.5*specSamplingRate),2);
%         baseMatrixWhiskHipS_Vals = baseWhiskHipS_Vals.*ones(size(meanWhiskHipS));
%         msWhiskHipS_Vals = (meanWhiskHipS - baseMatrixWhiskHipS_Vals);
        T2 = -20:(1/specSamplingRate):20;
        % save results
        AnalysisResults.(animalID).Whisk.Z_Ach.(whiskCriteriaName).Rhodamine.Rhodamine = AchmeanWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.Z_Ach.(whiskCriteriaName).Rhodamine.RhodamineStD = AchstdWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).Rhodamine.Rhodamine = NEmeanWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).Rhodamine.RhodamineStD = NEstdWhiskRhodamineData;

        AnalysisResults.(animalID).Whisk.Z_Ach.(whiskCriteriaName).GFP.GFP = AchmeanWhiskGFPData;
        AnalysisResults.(animalID).Whisk.Z_Ach.(whiskCriteriaName).GFP.GFPStD = AchstdWhiskGFPData;
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).GFP.GFP = NEmeanWhiskGFPData;
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).GFP.GFPStD = NEstdWhiskGFPData;

        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.hippocampalData = meanWhiskHippocampalGamData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.hippocampalStD = stdWhiskHippocampalGamData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.corticalS = msStimWhiskS_Vals;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.hippocampalS = msWhiskHipS_Vals;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.T = T2;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.F = F;
    end
    %% analyze stimulus-evoked responses
    if firstHrs == "true"
        for gg = 1:length(stimCriteriaNames)
            stimCriteriaName = stimCriteriaNames{1,gg};
            if strcmp(stimCriteriaName,'stimCriteriaA') == true
                StimCriteria = StimCriteriaA;
                solenoid = 'LPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
                StimCriteria = StimCriteriaB;
                solenoid = 'RPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
                StimCriteria = StimCriteriaC;
                solenoid = 'AudSol';
            end
            % pull data from EventData.mat structure
            allStimFilter = FilterEvents_IOS(EventData.Rhodamine.Z_NE.stim,StimCriteria);
            [AchallStimRhodamineData] = EventData.Rhodamine.Z_Ach.stim.NormData(allStimFilter,:);
            [NEallStimRhodamineData] = EventData.Rhodamine.Z_NE.stim.NormData(allStimFilter,:);
            [AchallStimGFPData] = EventData.GFP.Z_Ach.stim.NormData(allStimFilter,:);
            [NEallStimGFPData] = EventData.GFP.Z_NE.stim.NormData(allStimFilter,:); 
    
            [allStimCortMUAData] = EventData.cortical_LH.corticalPower.stim.NormData(allStimFilter,:);
%             [allStimHipMUAData] = EventData.hippocampus.corticalPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData] = EventData.cortical_LH.gammaBandPower.stim.NormData(allStimFilter,:);
%             [allStimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.Rhodamine.Z_NE.stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.Rhodamine.Z_NE.stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalStimRhodamineData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS(NEallStimRhodamineData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [AchfinalStimRhodamineData,~,~,~] = RemoveInvalidData_IOS(AchallStimRhodamineData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [NEfinalStimGFPData,~,~,~] = RemoveInvalidData_IOS(NEallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);        
            [AchfinalStimGFPData,~,~,~] = RemoveInvalidData_IOS(AchallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
%             [finalStimHipMUAData,~,~,~] = RemoveInvalidData_IOS(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
%             [finalStimHipGamData,~,~,~] = RemoveInvalidData_IOS(allStimHipGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procStimGFPData procStimRhodamineData procStimRhodamineData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
            ii = 1;
            for hh = 1:size(NEfinalStimRhodamineData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 0;%2;
                stimEndTime = stimStartTime + 15;%15;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    AchstimRhodaminearray = AchfinalStimRhodamineData(hh,:);
                    NEstimRhodaminearray = NEfinalStimRhodamineData(hh,:);
                    AchstimGFParray = AchfinalStimGFPData(hh,:);
                    NEstimGFParray = NEfinalStimGFPData(hh,:);
    
                    stimCortMUAarray = finalStimCortMUAData(hh,:);
%                     stimHipMUAarray = finalStimHipMUAData(hh,:);
                    stimCortGamArray = finalStimCortGamData(hh,:);
%                     stimHipGamArray = finalStimHipGamData(hh,:);
    
                    AchfiltStimRhodaminearray = sgolayfilt(AchstimRhodaminearray,3,9);
                    NEfiltStimRhodaminearray = sgolayfilt(NEstimRhodaminearray,3,9);
                    AchfiltStimGFParray = sgolayfilt(AchstimGFParray,3,9);
                    NEfiltStimGFParray = sgolayfilt(NEstimGFParray,3,9);
    
                    AchprocStimRhodamineData(hh,:) = AchfiltStimRhodaminearray - mean(AchfiltStimRhodaminearray(1:(offset*samplingRate)));
                    NEprocStimRhodamineData(hh,:) = NEfiltStimRhodaminearray - mean(NEfiltStimRhodaminearray(1:(offset*samplingRate)));
                    AchprocStimGFPData(hh,:) = AchfiltStimGFParray - mean(AchfiltStimGFParray(1:(offset*samplingRate)));
                    NEprocStimGFPData(hh,:) = NEfiltStimGFParray - mean(NEfiltStimGFParray(1:(offset*samplingRate)));
    
                    procStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
%                     procStimHipMUAData(hh,:) = stimHipMUAarray - mean(stimHipMUAarray(1:(offset*samplingRate)));
                    procStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
%                     procStimHipGamData(hh,:) = stimHipGamArray - mean(stimHipGamArray(1:(offset*samplingRate)));
                    finalStimStartTimes(ii,1) = stimStartTime;
                    finalStimEndTimes(ii,1) = stimEndTime;
                    finalStimFiles{ii,1} = finalStimFileID;
                    ii = ii + 1;
                end
            end
            AchmeanStimRhodamineData = mean(AchprocStimRhodamineData,1);%*100;
            AchstdStimRhodamineData = std(AchprocStimRhodamineData,0,1);%*100;
            NEmeanStimRhodamineData = mean(NEprocStimRhodamineData,1);%*100;
            NEstdStimRhodamineData = std(NEprocStimRhodamineData,0,1);%*100;
            AchmeanStimGFPData = mean(AchprocStimGFPData,1);%*100;
            AchstdStimGFPData = std(AchprocStimGFPData,0,1);%*100;
            NEmeanStimGFPData = mean(NEprocStimGFPData,1);%*100;
            NEstdStimGFPData = std(NEprocStimGFPData,0,1);%*100;
    
            meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
            stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
%             meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
%             stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
            meanStimCortGamData = mean(procStimCortGamData,1)*100;
            stdStimCortGamData = std(procStimCortGamData,0,1)*100;
%             meanStimHipGamData = mean(procStimHipGamData,1)*100;
%             stdStimHipGamData = std(procStimHipGamData,0,1)*100;
            % extract LFP from spectrograms associated with the stimuli indecies
            stimCortZhold = [];
%             stimHipZhold = [];
            for jj = 1:length(finalStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_LH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
%                         stimHippocampalS_Data = AllSpecData.hippocampus.normS{kk,1};
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
%                 stimHipS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
%                 stimHipZhold = cat(3,stimHipZhold,stimHipS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanStimCortS = mean(stimCortZhold,3);
            baseStimCortS_Vals = mean(meanStimCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixStimCortS_Vals = baseStimCortS_Vals.*ones(size(meanStimCortS));
            msStimCortS_Vals = (meanStimCortS - baseMatrixStimCortS_Vals);
            % hippocampal mean-subtract by first 2 seconds prior to stimulus
%             meanStimHipS = mean(stimHipZhold,3);
%             baseStimHipS_Vals = mean(meanStimHipS(:,1:1.5*specSamplingRate),2);
%             baseMatrixStimHipS_Vals = baseStimHipS_Vals.*ones(size(meanStimHipS));
%             msStimHipS_Vals = (meanStimHipS - baseMatrixStimHipS_Vals);
            % save results
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).count = size(procStimCortMUAData,1);
    
            AnalysisResults.(animalID).Stim.Z_Ach.(solenoid).Rhodamine.Rhodamine = AchmeanStimRhodamineData;
            AnalysisResults.(animalID).Stim.Z_Ach.(solenoid).Rhodamine.RhodamineStD = AchstdStimRhodamineData;
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).Rhodamine.Rhodamine = NEmeanStimRhodamineData;
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).Rhodamine.RhodamineStD = NEstdStimRhodamineData;
            AnalysisResults.(animalID).Stim.Z_Ach.(solenoid).GFP.GFP= AchmeanStimGFPData;
            AnalysisResults.(animalID).Stim.Z_Ach.(solenoid).GFP.GFPStD = AchstdStimGFPData;
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).GFP.GFP= NEmeanStimGFPData;
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).GFP.GFPStD = NEstdStimGFPData;
    
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalStD = stdStimCortMUAData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.hippocampalData = meanStimHipMUAData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalData = meanStimCortGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalStD = stdStimCortGamData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.hippocampalData = meanStimHipGamData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.hippocampalStD = stdStimHipGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.corticalS = msStimCortS_Vals;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.hippocampalS = msStimHipS_Vals;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.F = F;
        end
    end
 % save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end

    %% 

%     %% [1-S2o] Rhodamine auditory stim
%     figure;
% ax15 = subplot(1,2,1);
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.Rhodamine.Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.Rhodamine.Rhodamine   + AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.Rhodamine.RhodamineStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.Rhodamine.Rhodamine   - AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.Rhodamine.RhodamineStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% title('Auditory stim Blood Volume')
% ylabel('Z \DeltaF/F (Ach)')
% ax15.YLim = [-3 6];
% 
% yyaxis right
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.Rhodamine.Rhodamine  ,'-','color',colors('army green'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.Rhodamine.Rhodamine   + AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.Rhodamine.RhodamineStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.Rhodamine.Rhodamine   - AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.Rhodamine.RhodamineStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% ylabel('Z \DeltaF/F (NE)')
% ax15.YAxis(1).Color = colors('indian red');
% ax15.YAxis(2).Color = colors('army green');
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax15.TickLength = [0.03,0.03];
% ax15.YLim = [-3 6];
% 
% %% [1-S2r] GFP auditory stim
% ax18 = subplot(1,2,2);
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.GFP.GFP  ,'-','color',colors('indian red'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.GFP.GFP   + AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.GFP.GFPStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.GFP.GFP   - AnalysisResults_firstHrs.NEACh002.Stim.Z_Ach.AudSol.GFP.GFPStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% title('Auditory stim GFP')
% ylabel('Z \DeltaF/F GRAB Ach')
% ax18.YLim = [-3 6];
% 
% yyaxis right
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.GFP.GFP  ,'-','color',colors('army green'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.GFP.GFP   + AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.GFP.GFPStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.GFP.GFP   - AnalysisResults_firstHrs.NEACh002.Stim.Z_NE.AudSol.GFP.GFPStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% title('Auditory stim GFP')
% ylabel('Z \DeltaF/F GRAB NE')
% ax18.YAxis(1).Color = colors('indian red');
% ax18.YAxis(2).Color = colors('army green');
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax18.TickLength = [0.03,0.03];
% ax18.YLim = [-3 6];

end

