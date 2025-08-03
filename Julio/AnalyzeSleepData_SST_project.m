function [Results_SleepData] = AnalyzeSleepData_SST_project(animalID,group,setName,rootFolder,delim,Results_SleepData,days)
%________________________________________________________________________________________________________________________
% Written by Julio Flores-Cuadra
% The Pennsylvania State University, Dept. of Biology
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze sleep data
%________________________________________________________________________________________________________________________
%% Set folder
dataLocation = [rootFolder delim 'Data' delim group delim setName delim days delim animalID delim 'CombinedImaging'];
cd(dataLocation)

% Load scoring results
modelScoringResults = [animalID '_Manual_ScoringResults.mat'];
load(modelScoringResults)

%% Save and Get the percentage of time in each state
Results_SleepData.(group).(days).(animalID).numberOfScores = length(ScoringResults.alllabels);
Results_SleepData.(group).(days).(animalID).awakePercent   = round((sum(strcmp(ScoringResults.alllabels,'Not Sleep'))/length(ScoringResults.alllabels))*100,1);
Results_SleepData.(group).(days).(animalID).nremPercent    = round((sum(strcmp(ScoringResults.alllabels,'NREM Sleep'))/length(ScoringResults.alllabels))*100,1);
Results_SleepData.(group).(days).(animalID).remPercent     = round((sum(strcmp(ScoringResults.alllabels,'REM Sleep'))/length(ScoringResults.alllabels))*100,1);

%% Calculate the duration of NREM and REM events

bin_size     = 5;               % Each bin represents 5 seconds
state_values = [1, 2, 3];       % NREM = 1, REM = 2, Not Sleep = 3

% Define classification bins
%bin_edges  = [5, 10, 15, 20, 25, 30, Inf]; 5 seconds
bin_edges  = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, Inf]; 
%bin_edges  = [30, 45, 60, 75, 90, 105, 120, 135, 150, Inf]; 
%bin_labels = {'5-10s', '10-15s', '15-20s', '20-25s', '25-30s', '>30s'};
bin_labels = {'0-15s', '15-30s', '30-45s', '45-60s', '60-75s','75-90s', '90-105s','105-120s','', '135-150s', '>150s'};
%bin_labels = {'30-45s', '45-60s', '60-75s','75-90s', '90-105s','105-120s','120-135s', '135-150s', '>150s'};

% Initialize counts for NREM, REM, and NotSleep
nrem_counts     = zeros(1, length(bin_labels));
rem_counts      = zeros(1, length(bin_labels));
NotSleep_counts = zeros(1, length(bin_labels));

% Extract the labels
all_labels = ScoringResults.alllabels;

%% Calculate total duration time in sec and min 

% Convert labels to numerical values for processing
state_numeric = zeros(size(all_labels));
state_numeric(strcmp(all_labels, 'NREM Sleep')) = 1;
state_numeric(strcmp(all_labels, 'REM Sleep'))  = 2;
state_numeric(strcmp(all_labels, 'Not Sleep'))  = 3;

nBins_NREM     = sum(state_numeric == 1);
nBins_REM      = sum(state_numeric == 2);
nBins_NotSleep = sum(state_numeric == 3);

total_NREM_sec     = nBins_NREM     * bin_size;
total_REM_sec      = nBins_REM      * bin_size;
total_NotSleep_sec = nBins_NotSleep * bin_size;

% Duration in minutes
total_NREM_min     = total_NREM_sec     / 60;
total_REM_min      = total_REM_sec      / 60;
total_NotSleep_min = total_NotSleep_sec / 60;

% Save totals in seconds
Results_SleepData.(group).(days).(animalID).TotalNREM_Sec     = total_NREM_sec;
Results_SleepData.(group).(days).(animalID).TotalREM_Sec      = total_REM_sec;
Results_SleepData.(group).(days).(animalID).TotalNotSleep_Sec = total_NotSleep_sec;

% Save totals in minutes
Results_SleepData.(group).(days).(animalID).TotalNREM_Min     = total_NREM_min;
Results_SleepData.(group).(days).(animalID).TotalREM_Min      = total_REM_min;
Results_SleepData.(group).(days).(animalID).TotalNotSleep_Min = total_NotSleep_min;

%% Identify event durations for NREM and REM according to bin_labels

%Initialize empty arrays to collect bout information
NREM_bouts = [];  % Will be N×3: [start_time_sec, end_time_sec, duration_sec] for each NREM bout
REM_bouts  = [];  % Will be M×3: [start_time_sec, end_time_sec, duration_sec] for each REM bout

for i = 1:length(state_values)
    state = state_values(i);

    % Find transitions (start_indices/end_indices define consecutive runs)
    state_diff    = diff([0; state_numeric == state; 0]); 
    start_indices = find(state_diff == 1);
    end_indices   = find(state_diff == -1) - 1;

    % Compute durations (seconds per bout)
    durations = (end_indices - start_indices + 1) * bin_size;

    %Build a small "bout_matrix" = [start_sec, end_sec, duration_sec]
    %
    %  - Bin index j means:   bin j covers [(j−1)*bin_size … j*bin_size) seconds.
    %  - So if a bout starts at bin = 10, its start_time = (10−1)*5 = 45 sec.
    %  - If it ends at bin = 15, its end_time = 15*5 = 75 sec.
    %  - The bout’s duration is already (end_indices(j) − start_indices(j) + 1)*5.
    % 
    nBout   = length(start_indices);
    bout_matrix = zeros(nBout, 3);
    for j = 1:nBout
        bout_start_sec = (start_indices(j) - 1) * bin_size;
        bout_end_sec   = end_indices(j) * bin_size;
        bout_dur_sec   = durations(j);
        bout_matrix(j, :) = [bout_start_sec, bout_end_sec, bout_dur_sec];
    end

    % Classify durations into histogram bins 
    bin_indices = discretize(durations, bin_edges);

    % Accumulate counts AND append bout_matrix to the proper array
    if state == 1
        NREM_bouts = [NREM_bouts; bout_matrix];  %#ok<AGROW>
        for j = 1:length(bin_indices)
            if ~isnan(bin_indices(j))
                nrem_counts(bin_indices(j)) = nrem_counts(bin_indices(j)) + 1;
            end
        end

    elseif state == 2
        REM_bouts = [REM_bouts; bout_matrix];  %#ok<AGROW>
        for j = 1:length(bin_indices)
            if ~isnan(bin_indices(j))
                rem_counts(bin_indices(j)) = rem_counts(bin_indices(j)) + 1;
            end
        end

    elseif state == 3
        for j = 1:length(bin_indices)
            if ~isnan(bin_indices(j))
                NotSleep_counts(bin_indices(j)) = NotSleep_counts(bin_indices(j)) + 1;
            end
        end

    end
end

% Save NREM counts, REM counts, Awake counts, and Bin labels 
Results_SleepData.(group).(days).(animalID).NREM_counts = nrem_counts;
Results_SleepData.(group).(days).(animalID).REM_counts  = rem_counts;
Results_SleepData.(group).(days).(animalID).Awake_counts = NotSleep_counts;
Results_SleepData.(group).Bin_Labels = bin_labels;
% Each row of NREM_bouts is [start_sec, end_sec, duration_sec] for one NREM bout.
Results_SleepData.(group).(days).(animalID).NREM_Bouts = NREM_bouts;
% Each row of REM_bouts is [start_sec, end_sec, duration_sec] for one REM bout.
Results_SleepData.(group).(days).(animalID).REM_Bouts  = REM_bouts;

% Save data
cd([rootFolder delim 'Results_SST_project'])
save('Results_SleepData.mat','Results_SleepData')
cd([rootFolder delim 'Data'])
end
