function [AnalysisResults] = AnalyzeOptoStimEvokedResponses_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs,saveFigs)
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

% criteria for optostimulation
StimCriteriaO.Value = {'OptoStim'};
StimCriteriaO.Fieldname = {'solenoidName'};
StimCriteriaO.Comparison = {'equal'};
OptostimCriteriaNames = {'stimCriteriaO'};
%%
samplingRate = EventData.Rhodamine.Z_NE.whisk.samplingRate;
specSamplingRate = 10;
trialDuration_sec = EventData.Rhodamine.Z_NE.whisk.trialDuration_sec;
timeVector = (0:(EventData.Rhodamine.Z_NE.whisk.epoch.duration*samplingRate))/samplingRate - EventData.Rhodamine.Z_NE.whisk.epoch.offset;
offset = EventData.Rhodamine.Z_NE.whisk.epoch.offset;
%% analyze opto stimulus-evoked responses
    if firstHrs == "false"
        for gg = 1:length(OptostimCriteriaNames)
            stimCriteriaName = OptostimCriteriaNames{1,gg};
            if strcmp(stimCriteriaName,'stimCriteriaO') == true
                OptoStimCriteria = StimCriteriaO;
                solenoid = 'OptoStim';
            end
            % pull data from EventData.mat structure
            allOptoStimFilter = FilterEvents_IOS(EventData.Rhodamine.Z_NE.stim,OptoStimCriteria);
            [AchallOptoStimRhodamineData] = EventData.Rhodamine.Z_Ach.stim.NormData(allOptoStimFilter,:);
            [NEallOptoStimRhodamineData] = EventData.Rhodamine.Z_NE.stim.NormData(allOptoStimFilter,:);
            [AchallOptoStimGFPData] = EventData.GFP.Z_Ach.stim.NormData(allOptoStimFilter,:);
            [NEallOptoStimGFPData] = EventData.GFP.Z_NE.stim.NormData(allOptoStimFilter,:); 
    
            [allOptoStimCortMUAData] = EventData.cortical_LH.corticalPower.stim.NormData(allOptoStimFilter,:);
            [allOptoStimCortGamData] = EventData.cortical_LH.gammaBandPower.stim.NormData(allOptoStimFilter,:);
            [allOptoStimFileIDs] = EventData.Rhodamine.Z_NE.stim.fileIDs(allOptoStimFilter,:);
            [allOptoStimEventTimes] = EventData.Rhodamine.Z_NE.stim.eventTime(allOptoStimFilter,:);
            allOptoStimDurations = zeros(length(allOptoStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalOptoStimRhodamineData,finalOptoStimFileIDs,~,finalOptoStimFileEventTimes] = RemoveInvalidData_IOS(NEallOptoStimRhodamineData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [AchfinalOptoStimRhodamineData,~,~,~] = RemoveInvalidData_IOS(AchallOptoStimRhodamineData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [NEfinalOptoStimGFPData,~,~,~] = RemoveInvalidData_IOS(NEallOptoStimGFPData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);        
            [AchfinalOptoStimGFPData,~,~,~] = RemoveInvalidData_IOS(AchallOptoStimGFPData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [finalOptoStimCortMUAData,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortMUAData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [finalOptoStimCortGamData,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortGamData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procOptoStimGFPData procOptoStimRhodamineData procOptoStimRhodamineData procOptoStimCortMUAData procOptoStimCortGamData finalOptoStimStartTimes finalOptoStimEndTimes finalOptoStimFiles
            ii = 1;
            for hh = 1:size(NEfinalOptoStimRhodamineData,1)
                stimStartTime = round(finalOptoStimFileEventTimes(hh,1),1) - 0;%2;
                stimEndTime = stimStartTime + 15;%15;
                finalOptoStimFileID = finalOptoStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    AchstimRhodaminearray = AchfinalOptoStimRhodamineData(hh,:);
                    NEstimRhodaminearray = NEfinalOptoStimRhodamineData(hh,:);
                    AchstimGFParray = AchfinalOptoStimGFPData(hh,:);
                    NEstimGFParray = NEfinalOptoStimGFPData(hh,:);
    
                    stimCortMUAarray = finalOptoStimCortMUAData(hh,:);
                    stimCortGamArray = finalOptoStimCortGamData(hh,:);
    
                    AchfiltOptoStimRhodaminearray = sgolayfilt(AchstimRhodaminearray,3,9);
                    NEfiltOptoStimRhodaminearray = sgolayfilt(NEstimRhodaminearray,3,9);
                    AchfiltOptoStimGFParray = sgolayfilt(AchstimGFParray,3,9);
                    NEfiltOptoStimGFParray = sgolayfilt(NEstimGFParray,3,9);
    
                    AchprocOptoStimRhodamineData(hh,:) = AchfiltOptoStimRhodaminearray - mean(AchfiltOptoStimRhodaminearray(1:(offset*samplingRate)));
                    NEprocOptoStimRhodamineData(hh,:) = NEfiltOptoStimRhodaminearray - mean(NEfiltOptoStimRhodaminearray(1:(offset*samplingRate)));
                    AchprocOptoStimGFPData(hh,:) = AchfiltOptoStimGFParray - mean(AchfiltOptoStimGFParray(1:(offset*samplingRate)));
                    NEprocOptoStimGFPData(hh,:) = NEfiltOptoStimGFParray - mean(NEfiltOptoStimGFParray(1:(offset*samplingRate)));
    
                    procOptoStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
                    procOptoStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
                    finalOptoStimStartTimes(ii,1) = stimStartTime;
                    finalOptoStimEndTimes(ii,1) = stimEndTime;
                    finalOptoStimFiles{ii,1} = finalOptoStimFileID;
                    ii = ii + 1;
                end
            end
            AchmeanOptoStimRhodamineData = mean(AchprocOptoStimRhodamineData,1);%*100;
            AchstdOptoStimRhodamineData = std(AchprocOptoStimRhodamineData,0,1);%*100;
            NEmeanOptoStimRhodamineData = mean(NEprocOptoStimRhodamineData,1);%*100;
            NEstdOptoStimRhodamineData = std(NEprocOptoStimRhodamineData,0,1);%*100;
            AchmeanOptoStimGFPData = mean(AchprocOptoStimGFPData,1);%*100;
            AchstdOptoStimGFPData = std(AchprocOptoStimGFPData,0,1);%*100;
            NEmeanOptoStimGFPData = mean(NEprocOptoStimGFPData,1);%*100;
            NEstdOptoStimGFPData = std(NEprocOptoStimGFPData,0,1);%*100;
    
            meanOptoStimCortMUAData = mean(procOptoStimCortMUAData,1)*100;
            stdOptoStimCortMUAData = std(procOptoStimCortMUAData,0,1)*100;
            meanOptoStimCortGamData = mean(procOptoStimCortGamData,1)*100;
            stdOptoStimCortGamData = std(procOptoStimCortGamData,0,1)*100;
            % extract LFP from spectrograms associated with the stimuli indecies
            stimCortZhold = [];
            for jj = 1:length(finalOptoStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalOptoStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_LH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                        F = AllSpecData.(stimSpecField).F{kk,1};
                        T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                    end
                end
                stimStartTimeIndex = find(T == round(finalOptoStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalOptoStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanOptoStimCortS = mean(stimCortZhold,3);
            baseOptoStimCortS_Vals = mean(meanOptoStimCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixOptoStimCortS_Vals = baseOptoStimCortS_Vals.*ones(size(meanOptoStimCortS));
            msOptoStimCortS_Vals = (meanOptoStimCortS - baseMatrixOptoStimCortS_Vals);   
            % save results
            T2 = -20:(1/specSamplingRate):20;
            AnalysisResults.(animalID).OptoStim.Z_NE.(solenoid).count = size(procOptoStimCortMUAData,1);
    
            AnalysisResults.(animalID).OptoStim.Z_Ach.(solenoid).Rhodamine.Rhodamine = AchmeanOptoStimRhodamineData;
            AnalysisResults.(animalID).OptoStim.Z_Ach.(solenoid).Rhodamine.RhodamineStD = AchstdOptoStimRhodamineData;
            AnalysisResults.(animalID).OptoStim.Z_NE.(solenoid).Rhodamine.Rhodamine = NEmeanOptoStimRhodamineData;
            AnalysisResults.(animalID).OptoStim.Z_NE.(solenoid).Rhodamine.RhodamineStD = NEstdOptoStimRhodamineData;
            AnalysisResults.(animalID).OptoStim.Z_Ach.(solenoid).GFP.GFP= AchmeanOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.Z_Ach.(solenoid).GFP.GFPStD = AchstdOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.Z_NE.(solenoid).GFP.GFP= NEmeanOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.Z_NE.(solenoid).GFP.GFPStD = NEstdOptoStimGFPData;
    
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).MUA.corticalData = meanOptoStimCortMUAData;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).MUA.corticalStD = stdOptoStimCortMUAData;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).Gam.corticalData = meanOptoStimCortGamData;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).Gam.corticalStD = stdOptoStimCortGamData;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).LFP.corticalS = msOptoStimCortS_Vals;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).OptoStim.cortical.(solenoid).LFP.F = F;
        end
    end
 %% save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end

%% Opto stim
summaryFigureN = figure;
plot(AnalysisResults.(animalID).OptoStim.cortical.OptoStim.timeVector, AnalysisResults.(animalID).OptoStim.Z_Ach.OptoStim.Rhodamine.Rhodamine,'color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(AnalysisResults.(animalID).OptoStim.cortical.OptoStim.timeVector, AnalysisResults.(animalID).OptoStim.Z_NE.OptoStim.Rhodamine.Rhodamine,'color',[0.6350 0.0780 0.1840],'LineWidth',2)

plot(AnalysisResults.(animalID).OptoStim.cortical.OptoStim.timeVector, AnalysisResults.(animalID).OptoStim.Z_Ach.OptoStim.GFP.GFP,'color',[0 0.4470 0.7410],'LineWidth',2)

plot(AnalysisResults.(animalID).OptoStim.cortical.OptoStim.timeVector, AnalysisResults.(animalID).OptoStim.Z_NE.OptoStim.GFP.GFP,'color',[0.4660 0.6740 0.1880],'LineWidth',2)

legend('CBV LH','CBV RH','RH NE','NE')
ylabel('Z \DeltaF/F ')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
title('Optogenetic Stimulus Evoked Response')

% save figure(s)
delim = '\';
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'OptoStim-FiberSignals-consolidated']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'OptoStim-FiberSignals-consolidated'])
    close 
end
end

