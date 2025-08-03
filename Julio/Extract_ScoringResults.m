% Load the original ScoringResults data
load('JF043_Manual_ScoringResults.mat');

% Define the dates of interest
%dates_of_interest = {'240608', '240615', '240622', '240701'};
dates_of_interest = {'240512'};

% Convert fileIDs to a string format to allow for easy matching
fileIDs_str = cellfun(@(x) string(x), ScoringResults.fileIDs, 'UniformOutput', false);
allfileIDs_str = cellfun(@(x) string(x), ScoringResults.allfileIDs, 'UniformOutput', false);

% Create a mask to filter fileIDs and allfileIDs containing the desired dates
mask = false(size(fileIDs_str));
mask_all = false(size(allfileIDs_str));

for i = 1:length(fileIDs_str)
    for j = 1:length(dates_of_interest)
        if contains(fileIDs_str{i}, dates_of_interest{j})
            mask(i) = true;
            break;
        end
    end
end

for i = 1:length(allfileIDs_str)
    for j = 1:length(dates_of_interest)
        if contains(allfileIDs_str{i}, dates_of_interest{j})
            mask_all(i) = true;
            break;
        end
    end
end

% Filter fileIDs, labels, allfileIDs, and alllabels based on the masks
filtered_fileIDs = ScoringResults.fileIDs(mask);
filtered_labels = ScoringResults.labels(mask);
filtered_allfileIDs = ScoringResults.allfileIDs(mask_all);
filtered_alllabels = ScoringResults.alllabels(mask_all);

% Create a new structure to store the extracted data
ScoringResults.fileIDs = filtered_fileIDs;
ScoringResults.labels = filtered_labels;
ScoringResults.allfileIDs = filtered_allfileIDs;
ScoringResults.alllabels = filtered_alllabels;

% Save the new structure to a new .mat file
save('JF043_ScoringResults_Extracted_day3', 'ScoringResults');

