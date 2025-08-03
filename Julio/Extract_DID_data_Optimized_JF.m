clear; clc;

% Define the weeks (folder names)
weeks = {'Baseline1', 'DID1', 'DID2', 'DID3', 'DID4'};
mice = {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF050'};
regions = {'PFC', 'S1BF'};
transitions = {'AWAKEtoNREM', 'NREMtoAWAKE', 'NREMtoREM', 'REMtoAWAKE'};

% Initialize results table
results = {};
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
    
    % Process data for each mouse, region, and transition
    for iMouse = 1:length(mice)
        mouse = mice{iMouse};
        for iRegion = 1:length(regions)
            region = regions{iRegion};
            % Determine the data fields based on the region
            if strcmp(region, 'PFC')
                dataFields = {'GRAB_NE', 'NE_CBV'};
            elseif strcmp(region, 'S1BF')
                dataFields = {'GRAB_ACh', 'ACh_CBV'};
            end
            for iDataField = 1:length(dataFields)
                dataField = dataFields{iDataField};
                for iTransition = 1:length(transitions)
                    transition = transitions{iTransition};
                    % Try to extract the data
                    try
                        data = AnalysisResults.(mouse).Transitions.(transition).(dataField);
                    catch
                        % If data is missing, skip this combination
                        fprintf('Data missing for Week: %s, %s, %s, %s, %s\n', week, mouse, region, dataField, transition);
                        continue;
                    end
                    % Calculate means before and after the transition
                    mean_before = mean(data(1:900));
                    mean_after = mean(data(901:1800));
                    diff_mean = mean_after - mean_before;
                    % Store the results
                    results{counter, 1} = week;
                    results{counter, 2} = mouse;
                    results{counter, 3} = region;
                    results{counter, 4} = dataField;
                    results{counter, 5} = transition;
                    results{counter, 6} = mean_before;
                    results{counter, 7} = mean_after;
                    results{counter, 8} = diff_mean;
                    counter = counter + 1;
                end
            end
        end
    end
end

% Convert the results to a table for better readability
resultsTable = cell2table(results, 'VariableNames', {'Week', 'Mouse', 'Region', 'DataField', 'Transition', 'MeanBefore', 'MeanAfter', 'DiffMean'});

% Display the results table
disp(resultsTable);

save('AnalysisResults_DID.mat', 'resultsTable')
writetable(resultsTable, 'ResultsTable_AllWeeks.csv');


