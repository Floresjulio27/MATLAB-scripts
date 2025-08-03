clear; clc;

load('AnalysisResults.mat');

% Define mice, regions, and transitions
%mice = {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF050'};
week = {'Baseline', 'DID1', 'DID2'}; %Remember to change 
mice = {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF050'};
regions = {'PFC', 'S1BF'};
transitions = {'AWAKEtoNREM', 'NREMtoAWAKE', 'NREMtoREM', 'REMtoAWAKE'};

% Initialize results table
results = {};
counter = 1;

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
                    % If data is missing, skip this combination. Because
                    % there are some mice that do not have REM
                    fprintf('Data missing for %s, %s, %s, %s\n', mouse, region, dataField, transition);
                    continue;
                end
                % Calculate means before and after the transition
                mean_before = mean(data(1:900));
                mean_after = mean(data(901:1800));
                diff_mean = mean_after - mean_before;
                % Store the results
                results{counter, 1} = mouse;
                results{counter, 2} = region;
                results{counter, 3} = dataField;
                results{counter, 4} = transition;
                results{counter, 5} = mean_before;
                results{counter, 6} = mean_after;
                results{counter, 7} = diff_mean;
                counter = counter + 1;
            end
        end
    end
end

% Convert the results to a table for better readability
resultsTable = cell2table(results, 'VariableNames', {'Mouse', 'Region', 'DataField', 'Transition', 'MeanBefore', 'MeanAfter', 'DiffMean'});

% Display the results table
disp(resultsTable);


% % Initialize or load the existing structure if it exists
% if isfile('AnalysisResults_DID.mat')
%     load('AnalysisResults_DID.mat', 'AnalysisResults_DID'); 
% else
%     AnalysisResults_DID = struct(); 
% end
% 
% % Save the DiffMean results into AnalysisResults_DID, organized by mouse, week, region, and transition
% for iMouse = 1:length(mice)
%     mouse = mice{iMouse};
% 
%     % Check if the mouse already exists in the structure
%     if ~isfield(AnalysisResults_DID, mouse)
%         AnalysisResults_DID.(mouse) = struct(); % Initialize if the mouse field does not exist
%     end
% 
%     % Check if the current week exists in the structure
%     if ~isfield(AnalysisResults_DID.(mouse), week)
%         AnalysisResults_DID.(mouse).(week) = struct(); % Initialize if the week field does not exist
%     end
% 
%     for iRegion = 1:length(regions)
%         region = regions{iRegion};
% 
%         % Check if the region exists within the current week
%         if ~isfield(AnalysisResults_DID.(mouse).(week), region)
%             AnalysisResults_DID.(mouse).(week).(region) = struct(); % Initialize if the region field does not exist
%         end
% 
%         % Filter results for the current region
%         mouseResults = results(strcmp(results(:, 1), mouse), :);
%         regionResults = mouseResults(strcmp(mouseResults(:, 2), region), :);
% 
%         for iTransition = 1:length(transitions)
%             transition = transitions{iTransition};
% 
%             % Filter results for the current transition
%             transitionResults = regionResults(strcmp(regionResults(:, 4), transition), :);
% 
%             % Save the DiffMean value directly to the structure
%             if ~isempty(transitionResults)
%                 % Assuming only one DiffMean per combination of mouse, region, and transition
%                 diffMeanValue = transitionResults{1, 7}; % Extract DiffMean
%                 AnalysisResults_DID.(mouse).(week).(region).(transition) = diffMeanValue;
%             end
%         end
%     end
% end
% 
% % Save the updated structure 
% save('AnalysisResults_DID.mat', 'AnalysisResults_DID');



%% Plot some results 
