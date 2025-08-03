function [] = MainScript_SST_project()
%%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner

%Modified by: Julio Flores-Cuadra
%----------------------------------------------------------------------------------------------------------
currentFolder = pwd;
addpath(genpath(currentFolder));
zap;
saveState = true;
runAnalysis = true;
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
addpath(genpath(rootFolder))
if runAnalysis == true
    RunAnalysis_SST_project(rootFolder,delim)
end
% analysis subfunctions
disp('Loading analysis results and generating figures...'); disp(' ')

%_______MAIN FIGURE PANELS_______%

%% Alcohol consumption Data 
Fig_Alcohol_MeanGFP_relationship(rootFolder,saveState,delim)  %DONE I think will have to do this per inidividual day
%% Sleep Data
%Fig_Sleep_Histograms2_SST_project(rootFolder,saveState,delim);
%Fig_Total_Sleep_Duration_SST_project(rootFolder,saveState,delim); %Graph
%with days 
%Fig_Total_Sleep_Bouts_SST_project(rootFolder,saveState,delim)
%Fig_Total_Sleep_Bouts_Histfit_SST_project(rootFolder,saveState,delim);

%% Evoked Reponses Pooled
%Fig_Stim_SST_ResponseIndividual_Alcohol_ANIMALS(rootFolder,saveState,delim)
%Fig_Stim_SST_ResponseIndividual_WATER_ANIMALS(rootFolder,saveState,delim)

%Fig_Stim_SST_ResponseIndividual_Water(rootFolder,saveState,delim);
%Fig_Stim_SST_ResponseIndividual_Alcohol(rootFolder,saveState,delim)
%Fig_Stim_SST_ResponseIndividual_Alcohol_Stats(rootFolder,saveState,delim) %Stats are always calculated per animal 
%% Arousal State Transitions_Animal mode
%Fig2_FP_Transition_SST_PFC_GFP_Water(rootFolder,saveState,delim)
%Fig2_FP_Transition_SST_PFC_CBV_Water(rootFolder,saveState,delim)
%Fig2_FP_Transition_SST_PFC_GFP_Alcohol(rootFolder,saveState,delim) 
%Fig2_FP_Transition_SST_PFC_CBV_Alcohol(rootFolder,saveState,delim)

%Fig2_FP_Transition_SST_S1BF_GFP_Water(rootFolder,saveState,delim)
%Fig2_FP_Transition_SST_S1BF_CBV_Water(rootFolder,saveState,delim)
%Fig2_FP_Transition_SST_S1BF_GFP_Alcohol(rootFolder,saveState,delim)
%Fig2_FP_Transition_SST_S1BF_CBV_Alcohol(rootFolder,saveState,delim) 

%% Arousal State Transitions ANIMALS
%Fig_FP_Transition_SST_PFC_GFP_Alcohol_ANIMALS(rootFolder,saveState,delim)

%% MeanGFP data
%Fig_MeanGFP_NORM_AcrossDays_WaterVsAlcohol_SST_project(rootFolder,saveState,delim)%xnot
%all animals have baselines :(
%Fig_MeanGFP_RAW_AcrossDays_WaterVsAlcohol_SST_project(rootFolder,saveState,delim) %This. 
%Fig_MeanGFP_POOLED_NORM_AcrossDays_WaterVsAlcohol_SST_project(rootFolder,saveState,delim) %Pooled data is for visualization
%Fig_MeanGFP_POOLED_RAW_AcrossDays_WaterVsAlcohol_SST_project(rootFolder,saveState,delim)%x

%One measurement per animal do this
%% Mean CBV data 
%Fig_MeanCBV_RAW_AcrossDays_WaterVsAlcohol_SST_project(rootFolder,saveState,delim)

%%

%% WhiskerProbability
%FigS01_WhiskerProbabilityBefandAftPuff(rootFolder,saveState,delim)

% go back to original directory
cd(currentFolder)

function [] = RunAnalysis_SST_project(rootFolder,delim)

%%% Get Alcohol data from excel file 
%GetAlcoholData_NewVersion

%%% Sleep Data
%AnalyzeSleepData_SST_project_handler(rootFolder,delim,true) %Missing to
%re-run this one

%%% Evoked Responses
%AnalyzeEvokedResponses_SST_project_handler(rootFolder,delim,true)

%%% Arousal state Transitions
%AnalyzeArousalTransitions_SST_project_handler(rootFolder,delim,true)

%%% Analyze Whisker Probability
%AnalyzeWhiskerProbability_SST_project_handler(rootFolder,delim,true)

%%% Analyze Mean GFP across arousal states
%AnalyzeMeanGFP_SST_project_handler(rootFolder,delim,true)

%%% Analyze Mean CBV across arousal states 
%AnalyzeMeanCBV_SST_project_handler(rootFolder,delim,true)

%%% Power Spectrum
%AnalyzePowerSpectrum_CBVandGFP_SST_project_handler(rootFolder,delim,false) % In progress...
%AnalyzePowerSpectrum_ECOG_project_handler(rootFolder,delim,false) %DONE

%%% Coherence 

%%% Cross-Correlation
