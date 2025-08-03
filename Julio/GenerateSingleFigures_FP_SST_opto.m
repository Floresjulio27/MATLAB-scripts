function [figHandle,ax1,ax2,ax3,ax4,ax5] = GenerateSingleFigures_FP_SST_opto(procDataFileID,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain (modified)
% The Pennsylvania State University, Dept. of Biomedical Engineering
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure for a single n minute trial with rearranged subplots.
%          This version unifies the EEG and OptoStim events (originally ax4 and ax5) into a single subplot.
%________________________________________________________________________________________________________________________

%% Load file and gather information
load(procDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_FP(procDataFileID);
strDay = ConvertDate_FP(fileDate);

%% Setup Butterworth filter coefficients
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);

%% Process Signals
% Whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;

% Pupil area
filteredpupilarea = ProcData.data.Pupil.Diameter;

% Force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;

% EMG signals
EMG = ProcData.data.EMG.emg;
EMG_Signal = ProcData.data.EMG.emgSignal;

% Stimulations (we will only use OptoStim here)
OptoStim = ProcData.data.stimulations.OptoStim;

% (CBV and GFP processing is retained if needed later but will not be plotted)
ACh_CBV = ProcData.data.CBV.P_ACh;
filtACh_CBV = filtfilt(sos2,g2,ACh_CBV);
NE_CBV = ProcData.data.CBV.P_NE;
filtNE_CBV = filtfilt(sos2,g2,NE_CBV);
ACh_GFP = ProcData.data.GFP.P_ACh;
filtACh_GFP = filtfilt(sos2,g2,ACh_GFP);
NE_GFP = ProcData.data.GFP.P_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);

% EEG
EEG_LH = ProcData.data.cortical_LH.corticalSignal;

%% Load Spectrogram Data
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;

%% Compute Y-values for behavioral indices and stimulation markers
indecesMax = max([filtACh_CBV; filtACh_CBV; filtACh_GFP; filtACh_GFP]);
OptoStim_Yvals = 1.30 * max(indecesMax) * ones(size(OptoStim));  % For plotting OptoStim

whisking_Yvals = 1.10 * max(indecesMax) * ones(size(binWhiskers));
force_Yvals    = 1.20 * max(indecesMax) * ones(size(binForce));
forceInds = binForce .* force_Yvals;
whiskInds = binWhiskers .* whisking_Yvals;

% Replace zero values with NaN (for cleaner scatter plots)
forceInds(forceInds==0) = NaN;
whiskInds(whiskInds==0) = NaN;

%% Create Figure
figHandle = figure;
if strcmp(saveFigs,'y')
    figHandle.WindowState = 'minimized';
end

%% Define 5 Subplots in the following order:
%
% 1. Force Sensor and EMG
% 2. Whisker Angle and Raw EMG
% 3. Pupil Area and Behavioral Indices
% 4. Combined EEG and OptoStim Events
% 5. Spectrogram

% 1. Force Sensor and EMG
ax1 = subplot(5,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs, filtForceSensor, 'color', colors('north texas green'), 'LineWidth', 1);
title([animalID ' Behavioral characterization for ' fileID2])
ylabel('Force Sensor (Volts)')
xlim([0, ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs, EMG, 'color', 'k', 'LineWidth', 1);
ylabel('EMG (Volts^2)')
xlim([0, ProcData.notes.trialDuration_sec])
legend([p1,p2],'Force Sensor','EMG Power')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% 2. Whisker Angle and Raw EMG
ax2 = subplot(5,1,2);
p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs, filteredWhiskerAngle, 'color', colors('dark pink'), 'LineWidth', 1);
ylabel('Angle (deg)')
xlim([0, ProcData.notes.trialDuration_sec])
ylim([-20,60])
yyaxis right
p2 = plot((1:length(EMG_Signal))/ProcData.notes.dsFs, EMG_Signal, 'color', 'k', 'LineWidth', 1);
ylabel('EMG (pixels)')
xlim([0, ProcData.notes.trialDuration_sec])
legend([p1,p2],'Whisker Angle','EMG')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% 3. Pupil Area and Behavioral Indices
ax3 = subplot(5,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs, forceInds, '.', 'MarkerEdgeColor', colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs, whiskInds, '.', 'MarkerEdgeColor', colors('dark pink'));
p2 = plot((1:length(filteredpupilarea))/ProcData.notes.dsFs, filteredpupilarea, 'color', colors('deep carrot orange'), 'LineWidth', 1);
ylabel('Pupil Area (pixels)')
legend([s1,s2,p2],'Movement','Whisking','Pupil')
xlim([0, ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% 4. Combined EEG and OptoStim Events
ax4 = subplot(5,1,4);
scatter(OptoStim, OptoStim_Yvals,'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
hold on
% Plot EEG on left y-axis
yyaxis left
plot((1:length(EEG_LH))/ProcData.notes.dsFs, EEG_LH, 'color', colors('deep carrot orange'), 'LineWidth', 1);
ylabel('ECoG (V)')
xlim([0, ProcData.notes.trialDuration_sec])
% Plot OptoStim events on right y-axis
ylabel('ECoG')
title('EEG and OptoStim Events')
legend('OptoStim','ECoG','Location','best')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight


% 5. Spectrogram
ax5 = subplot(5,1,5);
Semilog_ImageSC(T, F, cortical_LHnormS, 'y')
axis xy
c5 = colorbar;
ylabel(c5, '\DeltaP/P (%)')
clim([-100,100])
ylabel('Frequency (Hz)')
xlim([0, ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('ECOG')
set(gca,'Yticklabel', [])

%% Link x-axes for synchronized zooming/panning
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

%% (Optional) Adjust Axes Positions
% The following block is optional and can be customized as needed.
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');

% Example adjustments:
ax5Pos(3:4) = ax1Pos(3:4);

ax1Pos(4) = ax1Pos(4) + 0.03;
ax2Pos(4) = ax2Pos(4) + 0.03;
ax3Pos(4) = ax3Pos(4) + 0.03;
ax4Pos(4) = ax4Pos(4) + 0.03;
ax5Pos(4) = ax5Pos(4) + 0.03;

ax1Pos(3) = ax1Pos(3) + 0.05;
ax2Pos(3) = ax2Pos(3) + 0.05;
ax3Pos(3) = ax3Pos(3) + 0.05;
ax4Pos(3) = ax4Pos(3) + 0.05;
ax5Pos(3) = ax5Pos(3) + 0.05;

ax1Pos(1) = ax1Pos(1) - 0.05;
ax2Pos(1) = ax2Pos(1) - 0.05;
ax3Pos(1) = ax3Pos(1) - 0.05;
ax4Pos(1) = ax4Pos(1) - 0.05;
ax5Pos(1) = ax5Pos(1) - 0.05;

ax1Pos(2) = ax1Pos(2) - 0.05;
ax2Pos(2) = ax2Pos(2) - 0.05;
ax3Pos(2) = ax3Pos(2) - 0.05;
ax4Pos(2) = ax4Pos(2) - 0.05;
ax5Pos(2) = ax5Pos(2) - 0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);

%% Save the Figure (if requested)
if strcmp(saveFigs,'y')
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle, [dirpath animalID '_' fileID '_SingleTrialFig']);
    saveas(figHandle, [dirpath animalID '_' fileID '_SingleTrialFig'], 'tiff')
end

end
