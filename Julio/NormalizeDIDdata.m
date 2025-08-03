clear; clc;
load('AnalysisResults_DID.mat');

% Define parameters for filtering
weekNames = {'Baseline1', 'DID1'};
transitions = {'AWAKEtoNREM', 'NREMtoAWAKE', 'NREMtoREM', 'REMtoAWAKE'};
region = 'PFC';
dataField = 'GRAB_NE';

% Initialize storage for differences
diffData = struct();

% Loop through weeks and transitions to extract data
for w = 1:length(weekNames)
    for t = 1:length(transitions)
        % Filter the table for current week and transition
        filteredData = table2struct(resultsTable(strcmp(resultsTable.Week, weekNames{w}) & ...
                                                 strcmp(resultsTable.Region, region) & ...
                                                 strcmp(resultsTable.DataField, dataField) & ...
                                                 strcmp(resultsTable.Transition, transitions{t}), :));
        
        % Extract DiffMean values or set to NaN if empty
        diffKey = sprintf('%s_%s', weekNames{w}, transitions{t});
        if isempty(filteredData)
            diffData.(diffKey) = NaN;
        else
            diffData.(diffKey) = arrayfun(@(x) x.DiffMean, filteredData)';
        end
    end
end

% Normalize data
normData = struct();
for t = 1:length(transitions)
    transitionKey = transitions{t};
    
    % Baseline normalization
    baselineKey = sprintf('Baseline1_%s', transitionKey);
    if all(isnan(diffData.(baselineKey))) % Check if baseline is NaN
        normData.(sprintf('NormBaseline_%s', transitionKey)) = NaN;
    else
        normData.(sprintf('NormBaseline_%s', transitionKey)) = ...
            diffData.(baselineKey) ./ diffData.(baselineKey);
    end
    
    % DID1 normalization
    didKey = sprintf('DID1_%s', transitionKey);
    if all(isnan(diffData.(didKey))) % Check if DID1 is NaN
        normData.(sprintf('NormDID1_%s', transitionKey)) = NaN;
    else
        normData.(sprintf('NormDID1_%s', transitionKey)) = ...
            diffData.(didKey) ./ diffData.(baselineKey);
    end
end

% Example: Access normalized AWAKEtoNREM data
disp('Normalized Baseline AWAKEtoNREM:');
disp(normData.NormBaseline_AWAKEtoNREM);

disp('Normalized DID1 AWAKEtoNREM:');
disp(normData.NormDID1_AWAKEtoNREM);

