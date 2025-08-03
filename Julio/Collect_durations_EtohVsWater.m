%% Function to Collect Data from a Given Folder
function [nrem_durations, rem_durations] = Collect_durations_EtohVsWater(folder, mice_list, bin_size, state_values)
    nrem_durations = [];
    rem_durations = [];
    for m = 1:length(mice_list)
        mouse_id = mice_list{m};
        mat_file_path = fullfile(folder, mouse_id, 'CombinedImaging', strcat(mouse_id, '_Manual_ScoringResults.mat'));

        if exist(mat_file_path, 'file')
            load(mat_file_path, 'ScoringResults');
            all_labels = ScoringResults.alllabels;

            % Convert labels to numerical values
            state_numeric = zeros(size(all_labels));
            state_numeric(strcmp(all_labels, 'NREM Sleep')) = 1;
            state_numeric(strcmp(all_labels, 'REM Sleep')) = 2;

            % Identify event durations for NREM and REM
            for i = 1:length(state_values)
                state = state_values(i);

                % Find transitions
                state_diff = diff([0; state_numeric == state; 0]);
                start_indices = find(state_diff == 1);
                end_indices = find(state_diff == -1) - 1;

                % Compute durations
                durations = (end_indices - start_indices + 1) * bin_size;

                % Store durations based on state
                if state == 1
                    nrem_durations = [nrem_durations; durations];
                elseif state == 2
                    rem_durations = [rem_durations; durations];
                end
            end
        else
            fprintf('File not found: %s\n', mat_file_path);
        end
    end
end
