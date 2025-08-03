function [] = Fig_Total_Sleep_Duration_SST_project(rootFolder,saveState,delim)
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
        allNREM_duration = [];
        allREM_duration  = [];
        
        % Loop through animals for this (group,day)
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa};

            if ~isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'TotalNREM_Min')
                continue
            end

            NREMpath = Results_SleepData.(currentGroup).(currentDay).(animalID).TotalNREM_Min;
            allNREM_duration = [allNREM_duration; NREMpath];   %#ok<AGROW> %concatonate vertically the data

            if ~isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'TotalREM_Min')
                continue
            end
            
            REMpath = Results_SleepData.(currentGroup).(currentDay).(animalID).TotalREM_Min;
            allREM_duration = [allREM_duration;  REMpath];    %#ok<AGROW>
        end
        
        % compute and store stats
        N = size(allNREM_duration,1);
        data.(currentGroup).(currentDay).NREM_duration_N_Exp  = N;
        data.(currentGroup).(currentDay).NREM_duration_Mean   = mean(allNREM_duration,   1);
        data.(currentGroup).(currentDay).NREM_duration_SEM    = std(allNREM_duration,0,1) / sqrt(N);

        N = size(allREM_duration,1);
        data.(currentGroup).(currentDay).REM_duration_N_Exp  = N;
        data.(currentGroup).(currentDay).REM_duration_Mean   = mean(allREM_duration,   1);
        data.(currentGroup).(currentDay).REM_duration_SEM    = std(allREM_duration,0,1) / sqrt(N);
    end
end

%% Arrange to plot Water vs Alcohol

nGroups = numel(groups);   % should be 2
nDays   = numel(days);     % e.g. 5

% Each row of these matrices will be one group (row 1 = Water, row 2 = Alcohol).
NREM_mean = nan(nGroups, nDays);
NREM_sem  = nan(nGroups, nDays);

% (c) Loop over groups and days, pulling the scalar Mean and SEM out of 'data'
for gi = 1:nGroups
    currentGroup = groups{gi};        % 'Water' or 'Alcohol'
    for di = 1:nDays
        currentDay = days{di};       % 'Day1', 'Day2', etc.
        
        % Check that the field actually exists:
        if isfield(data.(currentGroup), currentDay)
            % Assuming those fields are scalars (1×1):
            NREM_mean(gi,di) = data.(currentGroup).(currentDay).NREM_duration_Mean;
            NREM_sem( gi,di) = data.(currentGroup).(currentDay).NREM_duration_SEM;
            REM_mean(gi,di) = data.(currentGroup).(currentDay).REM_duration_Mean;
            REM_sem( gi,di) = data.(currentGroup).(currentDay).REM_duration_SEM;
        else
            % If you didn’t compute DayX for this group, leave as NaN or zero.
            NREM_mean(gi,di) = NaN;
            NREM_sem( gi,di) = NaN;
            REM_mean( gi,di) = NaN;
            REM_sem( gi,di) = NaN;
        end
    end
end


%% Plot NREM time duration water vs. alcohol 
%Colors
Water_color = [0 0 0];
Alcohol_color = [0.8667, 0.1098, 0.4667]; %Magenta

figure;
hold on;
x = 1:nDays;

% Plot Water group 
Water_plot = errorbar(x,NREM_mean(1,:), NREM_sem( 1,:), 'o-', 'Color', Water_color, 'MarkerFaceColor', Water_color, 'MarkerEdgeColor', Water_color, 'LineWidth', 1.5,'CapSize', 8 );

% Plot Alcohol group
Alcohol_plot = errorbar(x, NREM_mean(2,:), NREM_sem( 2,:), 'o-', 'Color', Alcohol_color, 'MarkerFaceColor', Alcohol_color, 'MarkerEdgeColor', Alcohol_color, 'LineWidth', 1.5,'CapSize', 8 );

xticks(x);
xticklabels({'BL (Week 1)','DID1 (Week 2)','DID2 (Week 3)', 'DID3 (Week 4)', 'DID4 (Week 5)'});
xtickangle(45);   % rotate them 45° so they do not overlap (tip from google)

% % Set reasonable y‐limits and y‐ticks
% %    (Adjust these numbers if your NREM duration units are very different.)
% ylim([0, max(NREM_mean(:)+NREM_sem(:)) + 5]);  
% % or, if you know the max should be e.g. [0 – 100], you can hard‐code that.
% yticks(0:10:100);  

ylabel('Mean NREM Duration (min)', 'FontSize', 12);
sgtitle('Total NREM Duration Water vs. Alcohol');
legend([Water_plot, Alcohol_plot], {'Water','Alcohol'}, 'Location','northwest', 'FontSize', 11);






