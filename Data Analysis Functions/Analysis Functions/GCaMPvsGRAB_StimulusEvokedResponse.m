saveFigs = 'y';
rootFolder = 'F:\GCaMP_SST.vs.GRAB'; % Replace with your actual root folder path

% Whicn animals do you want?
FPanimalIDs_mouse1 = {'JF048'};
FPanimalIDs_mouse2 = {'SST002'};

% Initialize structures to hold data for each mouse
data_mouse1 = struct();
data_mouse2 = struct();

%Copied from Kevins'code
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};

%You need to store the data of each mouse individually. In the other code
%it extracts and avergage them all together
% Extract data for mouse 1
for aa = 1:length(FPanimalIDs_mouse1)
    animalID = FPanimalIDs_mouse1{aa};
    for dd = 1:length(solenoidNames)
        solenoidName = solenoidNames{dd};
        data_mouse1.P_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.P_NE.(solenoidName).GFP.GFP;
        data_mouse1.P_ACh.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).GFP.GFP;
        data_mouse1.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
    end
end

% Extract data for mouse 2
for aa = 1:length(FPanimalIDs_mouse2)
    animalID = FPanimalIDs_mouse2{aa};
    for dd = 1:length(solenoidNames)
        solenoidName = solenoidNames{dd};
        data_mouse2.P_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.P_NE.(solenoidName).GFP.GFP;
        data_mouse2.P_ACh.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).GFP.GFP;
        data_mouse2.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
    end
end

% Compute means and 95% CI for GFP data
for dd = 1:length(solenoidNames)
    solenoidName = solenoidNames{dd};
    
    % Mouse 1 data
    data_mouse1.P_NE.(solenoidName).GFP_mean = mean(data_mouse1.P_NE.(solenoidName).GFP, 2);
    data_mouse1.P_ACh.(solenoidName).GFP_mean = mean(data_mouse1.P_ACh.(solenoidName).GFP, 2);
    
    n_mouse1 = size(data_mouse1.P_NE.(solenoidName).GFP, 2);
    ci95_factor_mouse1 = tinv(0.975, n_mouse1-1);
    
    data_mouse1.P_NE.(solenoidName).GFP_ci95 = ci95_factor_mouse1 * (std(data_mouse1.P_NE.(solenoidName).GFP, 0, 2) / sqrt(n_mouse1));
    data_mouse1.P_ACh.(solenoidName).GFP_ci95 = ci95_factor_mouse1 * (std(data_mouse1.P_ACh.(solenoidName).GFP, 0, 2) / sqrt(n_mouse1));
    
    % Mouse 2 data
    data_mouse2.P_NE.(solenoidName).GFP_mean = mean(data_mouse2.P_NE.(solenoidName).GFP, 2);
    data_mouse2.P_ACh.(solenoidName).GFP_mean = mean(data_mouse2.P_ACh.(solenoidName).GFP, 2);
    
    n_mouse2 = size(data_mouse2.P_NE.(solenoidName).GFP, 2);
    ci95_factor_mouse2 = tinv(0.975, n_mouse2-1);
    
    data_mouse2.P_NE.(solenoidName).GFP_ci95 = ci95_factor_mouse2 * (std(data_mouse2.P_NE.(solenoidName).GFP, 0, 2) / sqrt(n_mouse2));
    data_mouse2.P_ACh.(solenoidName).GFP_ci95 = ci95_factor_mouse2 * (std(data_mouse2.P_ACh.(solenoidName).GFP, 0, 2) / sqrt(n_mouse2));
end
%%
summaryFigureN = figure;
sgtitle('Stimulus evoked responses in fiber photometry signals - Individual Mice GFP Data')

color_mouse1 = [0 0.4470 0.7410]; % Blue for mouse 1
color_mouse2 = [0.4660 0.6740 0.1880]; % Dark green for mouse 2

% Plot for Contralateral Stimulation in S1BF Fiber
ax1 = subplot(2,3,1);
% Mouse 1 data
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_ACh.RPadSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_ACh.RPadSol.GFP_mean + data_mouse1.P_ACh.RPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_ACh.RPadSol.GFP_mean - data_mouse1.P_ACh.RPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
title('Contra stim S1BF fiber')
ylabel('\DeltaF/F SST GCaMP (%)')
ax1.YLim = [-1.5 1.5];

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_ACh.RPadSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
hold on
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_ACh.RPadSol.GFP_mean + data_mouse2.P_ACh.RPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_ACh.RPadSol.GFP_mean - data_mouse2.P_ACh.RPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
ylabel('\DeltaF/F GRAB-SST (%)')
ax1.YAxis(1).Color = color_mouse1;
ax1.YAxis(2).Color = color_mouse2;
xlabel('Peri-stimulus time (s)')
axis square
set(gca, 'box', 'off')
ax1.TickLength = [0.03, 0.03];
ax1.YLim = [-1.5 1.5];
xlim([-5 15])
axis square

% Plot for Ipsilateral Stimulation in S1BF Fiber
ax2 = subplot(2,3,2);
yyaxis left
% Mouse 1 data
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_ACh.LPadSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_ACh.LPadSol.GFP_mean + data_mouse1.P_ACh.LPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_ACh.LPadSol.GFP_mean - data_mouse1.P_ACh.LPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
ylim([-1.5 1.5])
ylabel('\DeltaF/F SST GCaMP (%)')

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_ACh.LPadSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
hold on
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_ACh.LPadSol.GFP_mean + data_mouse2.P_ACh.LPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_ACh.LPadSol.GFP_mean - data_mouse2.P_ACh.LPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
title('Ipsi stim S1BF fiber')
ylabel('\DeltaF/F GRAB SST (%)')
xlabel('Peri-stimulus time (s)')
ax2.YAxis(1).Color = color_mouse1;
ax2.YAxis(2).Color = color_mouse2;
xlim([-5 15])
ylim([-1.5 1.5])
axis square
set(gca, 'box', 'off')
ax2.TickLength = [0.03, 0.03];

% Plot for Auditory Stimulation in S1BF Fiber
ax3 = subplot(2,3,3);
yyaxis left
% Mouse 1 data
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_ACh.AudSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_ACh.AudSol.GFP_mean + data_mouse1.P_ACh.AudSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_ACh.AudSol.GFP_mean - data_mouse1.P_ACh.AudSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
ylim([-1.5 1.5])
ylabel('\DeltaF/F SST GCaMP (%)')

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_ACh.AudSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_ACh.AudSol.GFP_mean + data_mouse2.P_ACh.AudSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_ACh.AudSol.GFP_mean - data_mouse2.P_ACh.AudSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
title('Auditory stim S1BF fiber')
ylabel('\DeltaF/F GRAB SST (%)')
xlabel('Peri-stimulus time (s)')
ax3.YAxis(1).Color = color_mouse1;
ax3.YAxis(2).Color = color_mouse2;
xlim([-5 15])
ylim([-1.5 1.5])
axis square
set(gca, 'box', 'off')
ax3.TickLength = [0.03, 0.03];

% Plot for Contralateral Stimulation in PFC Fiber
ax4 = subplot(2,3,4);
yyaxis left
% Mouse 1 data
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_NE.LPadSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_NE.LPadSol.GFP_mean + data_mouse1.P_NE.LPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.LPadSol.timeVector, 2), data_mouse1.P_NE.LPadSol.GFP_mean - data_mouse1.P_NE.LPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
ylim([-1.5 1.5])
ylabel('\DeltaF/F SST GCaMP (%)')

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_NE.LPadSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_NE.LPadSol.GFP_mean + data_mouse2.P_NE.LPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.LPadSol.timeVector, 2), data_mouse2.P_NE.LPadSol.GFP_mean - data_mouse2.P_NE.LPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
title('Contra stim PFC fiber')
ylabel('\DeltaF/F GRAB SST (%)')
xlabel('Peri-stimulus time (s)')
ax4.YAxis(1).Color = color_mouse1;
ax4.YAxis(2).Color = color_mouse2;
xlim([-5 15])
ylim([-1.5 1.5])
axis square
set(gca, 'box', 'off')
ax4.TickLength = [0.03, 0.03];

% Plot for Ipsilateral Stimulation in PFC Fiber
ax5 = subplot(2,3,5);
yyaxis left
% Mouse 1 data
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_NE.RPadSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_NE.RPadSol.GFP_mean + data_mouse1.P_NE.RPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.RPadSol.timeVector, 2), data_mouse1.P_NE.RPadSol.GFP_mean - data_mouse1.P_NE.RPadSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
ylim([-1.5 1.5])
ylabel('\DeltaF/F SST GCaMP (%)')

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_NE.RPadSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_NE.RPadSol.GFP_mean + data_mouse2.P_NE.RPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.RPadSol.timeVector, 2), data_mouse2.P_NE.RPadSol.GFP_mean - data_mouse2.P_NE.RPadSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
title('Ipsi stim PFC fiber')
ylabel('\DeltaF/F GRAB SST (%)')
xlabel('Peri-stimulus time (s)')
ax5.YAxis(1).Color = color_mouse1;
ax5.YAxis(2).Color = color_mouse2;
xlim([-5 15])
ylim([-1.5 1.5])
axis square
set(gca, 'box', 'off')
ax5.TickLength = [0.03, 0.03];

% Plot for Auditory Stimulation in PFC Fiber
ax6 = subplot(2,3,6);
yyaxis left
% Mouse 1 data
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_NE.AudSol.GFP_mean, '-', 'color', color_mouse1, 'LineWidth', 2)
hold on
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_NE.AudSol.GFP_mean + data_mouse1.P_NE.AudSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
plot(mean(data_mouse1.cortical.AudSol.timeVector, 2), data_mouse1.P_NE.AudSol.GFP_mean - data_mouse1.P_NE.AudSol.GFP_ci95, '-', 'color', color_mouse1, 'LineWidth', 0.10)
ylim([-1.5 1.5])
ylabel('\DeltaF/F SST GCaMP (%)')

% Mouse 2 data
yyaxis right
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_NE.AudSol.GFP_mean, '-', 'color', color_mouse2, 'LineWidth', 2)
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_NE.AudSol.GFP_mean + data_mouse2.P_NE.AudSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
plot(mean(data_mouse2.cortical.AudSol.timeVector, 2), data_mouse2.P_NE.AudSol.GFP_mean - data_mouse2.P_NE.AudSol.GFP_ci95, '--', 'color', color_mouse2, 'LineWidth', 0.10)
title('Auditory stim PFC fiber')
ylabel('\DeltaF/F GRAB SST (%)')
xlabel('Peri-stimulus time (s)')
ax6.YAxis(1).Color = color_mouse1;
ax6.YAxis(2).Color = color_mouse2;
xlim([-5 15])
ylim([-1.5 1.5])
axis square
set(gca, 'box', 'off')
ax6.TickLength = [0.03, 0.03];

%% Save figure(s)
if strcmp(saveFigs, 'y') == true
    % Save the figure in the rootFolder
    savefig(summaryFigureN, fullfile(rootFolder, '_Stim-FiberSignals_ResponseIndividual_CI_SST.fig'));
end