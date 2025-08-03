function [] = MainScript_SST
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. eLife2020
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________
% addpath(genpath('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Documents\Research_Codes\FiberPhotometry\Data-Analysis-master'))
clear; clc; close all;
%% Animal IDs
%FP_animalIDs = {'JF048'};
%FP_animalIDs = {'JF037','JF038','JF039'};
%FP_animalIDs = {'JF037','JF039','JF040'};
%FP_animalIDs = {'JF048','SST002'};
%FP_animalIDs = {'JF037','JF038','JF039','JF040'};
%FP_animalIDs = {'JF037','JF038','JF039','JF040','JF048','JF050'};
FP_animalIDs = {'JF037','JF038','JF039','JF050','JF048'};
%FP_animalIDs = {'JF037','JF038','JF039','JF048','JF050'};
%FP_animalIDs = {'JF037','JF038','JF039'};
%% make sure the code repository and data are present in the current directory
firstHrs = "false";
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'y';
saveFigs = 'y';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
%     multiWaitbar('Analyzing sleep probability',0,'Color','B'); pause(0.25);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults] = AnalyzeData(rootFolder,FP_animalIDs);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end
%% main figure panels
%OldCodeShakVersion: Ones below this message
%[AnalysisResults] = Fig4_FP_Transition_SST_SingleMouse_consolidated_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = Fig1_S2_Stim_SST_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%% New functions
%[AnalysisResults] = Fig1_S2_Stim_SST_ResponseIndividual_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = GCaMPvsGRAB_test(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = Fig4_FP_Transition_SST_SingleMouse_consolidated_CI_S1BF(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = Fig4_FP_Transition_SST_SingleMouse_consolidated_CI_PFC(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = Fig4_FP_Transition_SST_SingleMouse_consolidated_CI_PFC_FP(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%[AnalysisResults] = Fig1_S5_Movement_SST_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
 [AnalysisResults] = Transition_SST_SingleMouse_consolidated_CI_PFC_FP(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData(rootFolder,FP_animalIDs)

saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end
% these data are not from first hours
firstHrs = "false";
%% Block [1] Analyze the arousal-state probability of trial duration and resting events (IOS)
% runFromStart = 'y';
% for aa = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,aa})) == false || isfield(AnalysisResults.(FP_animalIDs{1,aa}),'SleepProbability') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeAwakeProbability_GRBANE(FP_animalIDs{1,aa},saveFigs,rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing sleep probability','Value',aa/length(FP_animalIDs));
% end
%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% runFromStart = 'y';
% for pp = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Movement') == false ||isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
% end
%% Block [3] Analyze the transitions between different arousal-states (IOS)
runFromStart = 'n';
for dd = 1:length(FP_animalIDs)
    if isfield(AnalysisResults,(FP_animalIDs{1,dd})) == false || isfield(AnalysisResults.(FP_animalIDs{1,dd}),'Transitions') == false || strcmp(runFromStart,'y') == true
          [AnalysisResults] = AnalyzeTransitionalAverages_FP_Movements_GRABNE(FP_animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults,firstHrs);
    end
    multiWaitbar('Analyzing behavioral transitions triggered changes','Value',dd/length(FP_animalIDs));
end
%% Block [4] Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
% runFromStart = 'n';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeMeanRhodamine_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%                 [AnalysisResults] = AnalyzeMeanCBV_FP_GRABNE_Sleep(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing behavior evoked CBV','Value',ff/length(FP_animalIDs));
% end
%% Block [5] Analyze the hemodynamic signal GRAB during different arousal states (IOS)
% runFromStart = 'n';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanGFP') == false || strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%                 [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE_Sleep(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing behavior evoked GRAB Sensors','Value',ff/length(FP_animalIDs));
% end
end