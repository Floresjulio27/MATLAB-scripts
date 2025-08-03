function [] = Fig_Total_Sleep_Bouts_SST_project(rootFolder,saveState,delim)
%________________________________________________________________________________________________________________________
% Written by Julio Flores-Cuadra
% The Pennsylvania State University, Dept. of Biology
%
% Acknowledgments. This code was written following scripts from Dr. Kevin
% Turner and MD Shakhawat Hossain
%________________________________________________________________________________________________________________________

%% set-up and process data
path = [rootFolder delim 'Results_SST_project']; % Path to the folder containing Results_SleepData.mat
cd(path)
resultsStruct = 'Results_SleepData';
load(resultsStruct)

groups = {'Water', 'Alcohol'};
days = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};


% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{ee};  
    
    % Loop through days
    for xx = 1:length(days)
        currentDay = days{xx};   
        
        % Check if this group actually has data for this day
        if ~isfield(Results_SleepData.(currentGroup), currentDay)
            continue
        end
        
        FPanimalIDs = fieldnames(Results_SleepData.(currentGroup).(currentDay));
        if isempty(FPanimalIDs)
            continue
        end

        %Pre-allocate
        allNREM_bouts = [];
        allREM_bouts  = [];
        
        % Loop through animals for this (group,day)
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa};

            if ~isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'NREM_Bouts')
                continue
            end

            NREMpath = Results_SleepData.(currentGroup).(currentDay).(animalID).NREM_Bouts(:,3);
            allNREM_bouts = [allNREM_bouts; NREMpath];   %#ok<AGROW> %concatonate vertically the data

            if ~isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'NREM_Bouts')
                continue
            end
            
            REMpath = Results_SleepData.(currentGroup).(currentDay).(animalID).REM_Bouts(:,3);
            allREM_bouts = [allREM_bouts;  REMpath];    %#ok<AGROW>
        end
        
        % compute and store stats
        N = size(allNREM_bouts,1);
        data.(currentGroup).(currentDay).NREM_bouts_N_Exp  = N;
        data.(currentGroup).(currentDay).NREM_bouts_Mean   = mean(allNREM_bouts,   1);
        data.(currentGroup).(currentDay).NREM_bouts_SEM    = std(allNREM_bouts,0,1) / sqrt(N);

        N = size(allREM_bouts,1);
        data.(currentGroup).(currentDay).REM_bouts_N_Exp  = N;
        data.(currentGroup).(currentDay).REM_bouts_Mean   = mean(allREM_bouts,   1);
        data.(currentGroup).(currentDay).REM_bouts_SEM    = std(allREM_bouts,0,1) / sqrt(N);
    end
end

%% Arrange to plot Water vs Alcohol

nGroups = numel(groups);   % should be 2
nDays   = numel(days);     % e.g. 5

% Each row of these matrices will be one group (row 1 = Water, row 2 = Alcohol).
NREM_mean = nan(nGroups, nDays);
NREM_sem  = nan(nGroups, nDays);
REM_mean = nan(nGroups, nDays);
REM_sem  = nan(nGroups, nDays);

% (c) Loop over groups and days, pulling the scalar Mean and SEM out of 'data'
for gi = 1:nGroups
    currentGroup = groups{gi};       
    for di = 1:nDays
        currentDay = days{di};       
        
        % Check that the field actually exists:
        if isfield(data.(currentGroup), currentDay)
            % Assuming those fields are scalars (1×1):
            NREM_mean(gi,di) = data.(currentGroup).(currentDay).NREM_bouts_Mean;
            NREM_sem( gi,di) = data.(currentGroup).(currentDay).NREM_bouts_SEM;
            REM_mean(gi,di) = data.(currentGroup).(currentDay).REM_bouts_Mean;
            REM_sem( gi,di) = data.(currentGroup).(currentDay).REM_bouts_SEM;
        else
            % If you didn’t compute DayX for this group, leave as NaN or zero.
            NREM_mean(gi,di) = NaN;
            NREM_sem( gi,di) = NaN;
            REM_mean( gi,di) = NaN;
            REM_sem( gi,di) = NaN;
        end
    end
end


%% Plot NREM time bouts water vs. alcohol 
%Colors
Water_color = [0 0 0];
Alcohol_color = [0.8667, 0.1098, 0.4667]; %Magenta

figure;
hold on;
x = 1:nDays;

% Plot Water group 
Water_plot = errorbar(x, NREM_mean(1,:), NREM_sem( 1,:), 'o-', 'Color', Water_color, 'MarkerFaceColor', Water_color, 'MarkerEdgeColor', Water_color, 'LineWidth', 1.5,'CapSize', 8 );

% Plot Alcohol group
Alcohol_plot = errorbar(x, NREM_mean(2,:), NREM_sem( 2,:), 'o-', 'Color', Alcohol_color, 'MarkerFaceColor', Alcohol_color, 'MarkerEdgeColor', Alcohol_color, 'LineWidth', 1.5,'CapSize', 8 );

xticks(x);
xticklabels({'BL (Week 1)','DID1 (Week 2)','DID2 (Week 3)', 'DID3 (Week 4)', 'DID4 (Week 5)'});
xtickangle(45);   % rotate them 45° so they do not overlap (tip from google)
xlim([0.5, nDays + 0.5]);

% % Set reasonable y‐limits and y‐ticks
% %    (Adjust these numbers if your NREM bouts units are very different.)
% ylim([0, max(NREM_mean(:)+NREM_sem(:)) + 5]);  
% % or, if you know the max should be e.g. [0 – 100], you can hard‐code that.
% yticks(0:10:100);  

ylabel('Mean NREM Duration (min)', 'FontSize', 12);
sgtitle('Total NREM Duration Water vs. Alcohol');
legend([Water_plot, Alcohol_plot], {'Water','Alcohol'}, 'Location','northwest', 'FontSize', 11);

%% Histfit water
% Define bin edges/centers EXACTLY as you originally binned
%    (Change these edges to match your real bin definitions!)

waterIdx = strcmp(groups, 'Water');  % logical index
waterIdx = find(waterIdx); % 1 should be for water

% Extract the matrix (nBins × nDays) for Water → NREM
allNREM_Water = NREM_summary(:, :, waterIdx);  % size: [nBins × nDays]
totalNREM_allBins = sum(allNREM_Water, 2);            % [nBins×1]
binEdges   = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 300];  % example
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;  

% 3) Expand counts into a (pseudo) vector of durations
pseudoRaw = repelem(binCenters, totalNREM_allBins.');  

% 4) Plot with histfit
figure;
histfit(pseudoRaw, numel(binCenters));
xlim([0, max(pseudoRaw)*1.1])
%histfit(pseudoRaw, numel(binCenters), 'loglogistic')
xlabel('NREM Event Duration (s)');
ylabel('Count');
title('Water Group → NREM Durations (Day1–Day6) with Normal Fit');

%% Histfit alcohol
% Define bin edges/centers EXACTLY as you originally binned
%    (Change these edges to match your real bin definitions!)
binEdges   = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 300];  % example
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;  

allNREM_Alcohol = NREM_summary(:, :, alcoholIdx);  % size: [nBins × nDays]
totalNREM_allBins_Alcohol = sum(allNREM_Alcohol, 2);            % [nBins×1]

% 3) Expand counts into a (pseudo) vector of durations
pseudoRaw = repelem(binCenters,totalNREM_allBins_Alcohol.');  

% 4) Plot with histfit
figure;
histfit(pseudoRaw, numel(binCenters)); 
xlim([0, max(pseudoRaw)*1.1])
xlabel('NREM Event Duration (s)');
ylabel('Count');
title('Water Group → NREM Durations (Day1–Day6) with Normal Fit');
