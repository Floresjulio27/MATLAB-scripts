clear;clc;
load('AnalysisResults.mat');

% Define mice, stimulation, and variable
mice = {'JF037', 'JF038', 'JF039', 'JF040'};
weeks = {'Baseline1', 'DID1', 'DID2', 'DID3', 'DID4'};
mice = {'JF037', 'JF038','JF039', 'JF040'};
regions = {'PFC'};
stimulations = {'LPadSol', 'RPadSol', 'AudSol'};
%stimulations = {'RPadSol'};

% Initialize the structure where you want the variables of interest
WhiskerPuff_PFC = {};
counter = 1; %This is for creating the table. It is basically to give a number to the rows so it can move in every loop

% Loop through each week
for iWeek = 1:length(weeks)
    week = weeks{iWeek};
    folderPath = fullfile(week); % Path to the week's folder
    dataFile = fullfile(folderPath, 'AnalysisResults.mat'); % Path to the analysis file

    % Load the AnalysisResults for the current week
    if exist(dataFile, 'file')
        load(dataFile); % Loads AnalysisResults
    else
        fprintf('File not found for week: %s\n', week);
        continue;
    end

    for iMouse = 1:length(mice)
        mouse = mice{iMouse};
        for iRegion = 1:length(regions)
            region = regions{iRegion};
            % Determine the data fields based on the region
            if strcmp(region, 'PFC')
                dataFields = {'P_NE'};
            end
            for iDataField = 1:length(dataFields)
                dataField = dataFields{iDataField};
                for iStimulation = 1:length(stimulations)
                    stimulation = stimulations{iStimulation};
                    % Try to extract the data
                    try
                        %data = AnalysisResults.(mouse).Stim.(dataField).(stimulation).CBV.CBV; %Change here for CBV or GFP
                        data = AnalysisResults.(mouse).Stim.(dataField).(stimulation).GFP.GFP;
                    catch
                        % If data is missing, skip this combination
                        fprintf('Data missing for %s, %s, %s, %s\n', mouse, region, dataField, stimulation);
                        continue;
                    end

                    % Calculate means before and after the transition
                    mean_before = mean(data(391:451));
                    %mean_after = mean(data(452:601)); %5sec
                    mean_after = mean(data(452:511)); %2sec
                    diff_mean = mean_after - mean_before;

                    % Store the results
                    WhiskerPuff_PFC{counter, 1} = week;
                    WhiskerPuff_PFC{counter, 2} = mouse;
                    WhiskerPuff_PFC{counter, 3} = region;
                    WhiskerPuff_PFC{counter, 4} = dataField;
                    WhiskerPuff_PFC{counter, 5} = stimulation;
                    WhiskerPuff_PFC{counter, 6} = mean_before;
                    WhiskerPuff_PFC{counter, 7} = mean_after;
                    WhiskerPuff_PFC{counter, 8} = diff_mean;
                    counter = counter + 1;
                end
            end
        end
    end
end

regions = {'S1BF'};
WhiskerPuff_S1BF = {};
counter = 1;

% Loop through each week
for iWeek = 1:length(weeks)
    week = weeks{iWeek};
    folderPath = fullfile(week); % Path to the week's folder
    dataFile = fullfile(folderPath, 'AnalysisResults.mat'); % Path to the analysis file

    % Load the AnalysisResults for the current week
    if exist(dataFile, 'file')
        load(dataFile); % Loads AnalysisResults
    else
        fprintf('File not found for week: %s\n', week);
        continue;
    end

    for iMouse = 1:length(mice)
        mouse = mice{iMouse};
        for iRegion = 1:length(regions)
            region = regions{iRegion};
            % Determine the data fields based on the region
            if strcmp(region, 'S1BF')
                dataFields = {'P_ACh'};
            end
            for iDataField = 1:length(dataFields)
                dataField = dataFields{iDataField};
                for iStimulation = 1:length(stimulations)
                    stimulation = stimulations{iStimulation};
                    % Try to extract the data
                    try
                        %data = AnalysisResults.(mouse).Stim.(dataField).(stimulation).CBV.CBV;
                        data = AnalysisResults.(mouse).Stim.(dataField).(stimulation).GFP.GFP;
                    catch
                        % If data is missing, skip this combination
                        fprintf('Data missing for %s, %s, %s, %s\n', mouse, region, dataField, stimulation);
                        continue;
                    end

                    % Calculate means before and after the transition
                    mean_before = mean(data(391:451));
                    mean_after = mean(data(452:601));
                    diff_mean = mean_after - mean_before;

                    % Store the results
                    WhiskerPuff_S1BF{counter, 1} = week;
                    WhiskerPuff_S1BF{counter, 2} = mouse;
                    WhiskerPuff_S1BF{counter, 3} = region;
                    WhiskerPuff_S1BF{counter, 4} = dataField;
                    WhiskerPuff_S1BF{counter, 5} = stimulation;
                    WhiskerPuff_S1BF{counter, 6} = mean_before;
                    WhiskerPuff_S1BF{counter, 7} = mean_after;
                    WhiskerPuff_S1BF{counter, 8} = diff_mean;
                    counter = counter + 1;
                end
            end
        end
    end
end



% Convert the results to a table for better readability
WhiskerPuff_PFC_Table = cell2table(WhiskerPuff_PFC, 'VariableNames', {'Week', 'Mouse', 'Region', 'DataField', 'Stimulation', 'MeanBefore', 'MeanAfter', 'DiffMean'});

WhiskerPuff_S1BF_Table = cell2table(WhiskerPuff_S1BF, 'VariableNames', {'Week', 'Mouse', 'Region', 'DataField', 'Stimulation', 'MeanBefore', 'MeanAfter', 'DiffMean'});









