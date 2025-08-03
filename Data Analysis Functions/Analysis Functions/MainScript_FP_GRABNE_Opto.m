function [] = MainScript_FP_GRABNE_Opto()
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
FP_animalIDs =    {'NEChR2001'};
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
% [AnalysisResults] = Fig1_S4_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
% [AnalysisResults] = Fig1_S3_Whisk_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
% [AnalysisResults] = Fig1_S2_Stim_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig1_S5_Movement_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);

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

%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'y';
for pp = 1:length(FP_animalIDs)
    if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeOptoStimEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs,saveFigs);
    end
    multiWaitbar('Analyzing optostim evoked responses','Value',pp/length(FP_animalIDs));
end

%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
