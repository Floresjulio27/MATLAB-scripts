% Initialize or load the existing structure if it exists
if isfile('AnalysisResults_DID.mat')
    load('AnalysisResults_DID.mat', 'AnalysisResults_DID'); % Load existing data
else
    AnalysisResults_DID = struct(); % Create a new structure if it doesn't exist
end

% Save the DiffMean results into AnalysisResults_DID, organized by mouse, week, region, and transition
for iMouse = 1:length(mice)
    mouse = mice{iMouse};
    
    % Check if the mouse already exists in the structure
    if ~isfield(AnalysisResults_DID, mouse)
        AnalysisResults_DID.(mouse) = struct(); % Initialize if the mouse field does not exist
    end
    
    % Check if the current week exists in the structure
    if ~isfield(AnalysisResults_DID.(mouse), week)
        AnalysisResults_DID.(mouse).(week) = struct(); % Initialize if the week field does not exist
    end
    
    for iRegion = 1:length(regions)
        region = regions{iRegion};
        
        % Check if the region exists within the current week
        if ~isfield(AnalysisResults_DID.(mouse).(week), region)
            AnalysisResults_DID.(mouse).(week).(region) = struct(); % Initialize if the region field does not exist
        end
        
        % Filter results for the current region
        mouseResults = results(strcmp(results(:, 1), mouse), :);
        regionResults = mouseResults(strcmp(mouseResults(:, 2), region), :);
        
        for iTransition = 1:length(transitions)
            transition = transitions{iTransition};
            
            % Filter results for the current transition
            transitionResults = regionResults(strcmp(regionResults(:, 4), transition), :);
            
            % Save the DiffMean value directly to the structure
            if ~isempty(transitionResults)
                % Assuming only one DiffMean per combination of mouse, region, and transition
                diffMeanValue = transitionResults{1, 7}; % Extract DiffMean
                AnalysisResults_DID.(mouse).(week).(region).(transition) = diffMeanValue;
            end
        end
    end
end

% Save the updated structure 
save('AnalysisResults_DID.mat', 'AnalysisResults_DID');