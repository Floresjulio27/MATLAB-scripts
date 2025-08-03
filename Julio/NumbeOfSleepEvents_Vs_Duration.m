clear; clc;

%% Define the parent folder where each mouse folder is located
parent_folder = cd; % Change to the actual path

% List of mouse IDs (folder names)
%mice_ids = {'JF048', 'JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'};
%mice_ids = {'JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'};% Add all your mice IDs here
mice_ids = {'JF037', 'JF038', 'JF039', 'JF040', 'JF044', 'JF045', 'JF046'};

bin_size = 5; % Each bin represents 5 seconds
state_values = [1, 2]; % NREM = 1, REM = 2

% Define classification bins
bin_edges = [5, 10, 15, 20, 25, 30, Inf]; 
bin_labels = {'5-10s', '10-15s', '15-20s', '20-25s', '25-30s', '>30s'};

%% Initialize counts for NREM and REM
nrem_counts = zeros(1, length(bin_labels));
rem_counts = zeros(1, length(bin_labels));

for m = 1:length(mice_ids)
    mouse_id = mice_ids{m};
    mat_file_path = fullfile(parent_folder, mouse_id, 'CombinedImaging', strcat(mouse_id, '_Manual_ScoringResults.mat'));
    
    % Check if the file exists
    if exist(mat_file_path, 'file')
        % Load the .mat file
        load(mat_file_path, 'ScoringResults');
        
        % Extract the labels
        all_labels = ScoringResults.alllabels;

        % Convert labels to numerical values for processing
        state_numeric = zeros(size(all_labels)); % Initialize
        state_numeric(strcmp(all_labels, 'NREM Sleep')) = 1;
        state_numeric(strcmp(all_labels, 'REM Sleep')) = 2;

        % Identify event durations for NREM and REM
        for i = 1:length(state_values)
            state = state_values(i);
            
            % Find transitions
            state_diff = diff([0; state_numeric == state; 0]); % Add 0s at edges % Shak functions
            start_indices = find(state_diff == 1);
            end_indices = find(state_diff == -1) - 1;
            
            % Compute durations
            durations = (end_indices - start_indices + 1) * bin_size; % Convert bins to seconds
            
            % Classify durations into bins
            bin_indices = discretize(durations, bin_edges);
            
            % Accumulate counts
            if state == 1
                for j = 1:length(bin_indices)
                    if ~isnan(bin_indices(j))
                        nrem_counts(bin_indices(j)) = nrem_counts(bin_indices(j)) + 1;
                    end
                end
            elseif state == 2
                for j = 1:length(bin_indices)
                    if ~isnan(bin_indices(j))
                        rem_counts(bin_indices(j)) = rem_counts(bin_indices(j)) + 1;
                    end
                end
            end
        end
    else
        fprintf('File not found: %s\n', mat_file_path);
    end
end

%% Plot

% Grouped bar chart (NREM vs REM together)
figure;
bar_data = [nrem_counts; rem_counts]';
bar(bar_data, 'grouped');
xticklabels(bin_labels);
xlabel('Sleep Event Duration Categories');
ylabel('Number of Events');
title('Distribution of Sleep Event Durations');
legend({'NREM', 'REM'}, 'Location', 'northwest');
grid on;

% Separate bar chart for NREM events
figure;
bar(nrem_counts, 'b');
xticklabels(bin_labels);
xlabel('Sleep Event Duration Categories');
ylabel('Number of NREM Events');
title('NREM Sleep Event Duration Distribution');
grid on;

% Separate bar chart for REM events
figure;
bar(rem_counts, 'r');
xticklabels(bin_labels);
xlabel('Sleep Event Duration Categories');
ylabel('Number of REM Events');
title('REM Sleep Event Duration Distribution');
grid on;



