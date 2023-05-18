
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Create ProcData structure using threshholds for the observed data
%            2) Extract the average pixel reflectance changes within those ROIs and save to RawData/ProcData files
%            3) Use spectral analysis of the reflectance data to pull out the animal's heart rate
%            4) Remove pixel-drift trends from each file, if necessary
%________________________________________________________________________________________________________________________   

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
% zap;
function StageTwoProcessing_FP_GRABNE_SingleFiber()
clc; clear all; close all
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
[animalID,~,~] = GetFileInfo_FP(rawDataFileIDs(1,:));
% Identify whether this analysis is for bilateral or single hemisphere imaging
curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
%% BLOCK PURPOSE: [1] Fiber Photometry data
disp('Analyzing Block [1] Processing TDT fiber photometry data.'); disp(' ')
% generate locations for fiberData
filelocation = pwd;
% Plot_Raw_TDTFiberPhotometry_bilateralGCaMP7s(filelocation,rawDataFiles)

TDTFiberPhotometry_bilateralGRABNE(filelocation,rawDataFiles)

% TDTFiberPhotometry_bilateralGCaMP7s_CheckFilter(filelocation,rawDataFiles)
%% BLOCK PURPOSE: [2] Correct the offset between the MScan and LabVIEW acquisiton.
disp('Analyzing Block [2] Correcting LabVIEW time offset.'); disp(' ')
fiberDataFileStruct = dir('*_FiberData.mat');
fiberDataFiles = {fiberDataFileStruct.name}';
fiberDataFileID = char(fiberDataFiles);
TemplateMatchFiberData_FP_GRABNE_SingleFiber(fiberDataFileID,rawDataFileIDs)
%% BLOCK PURPOSE: [3] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [3] Creating ProcData files and processing analog data.'); disp(' ')
ProcessRawDataFiles_FP_GRABNE_SingleFiber(rawDataFileIDs)
%% fin.
disp('Fiber Photometry Stage Two Processing - Complete.'); disp(' ')
