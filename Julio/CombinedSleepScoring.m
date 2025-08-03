% Load the files
combinedData = load('JF046_Manual_ScoringResults_Combined.mat');
did1Data = load('JF046_Manual_ScoringResults_DID1.mat');

% Assuming the structure is named 'ScoringResults'
% Concatenate the fields of the structures

combinedData.ScoringResults.fileIDs = [combinedData.ScoringResults.fileIDs; did1Data.ScoringResults.fileIDs];
combinedData.ScoringResults.labels = [combinedData.ScoringResults.labels; did1Data.ScoringResults.labels];
combinedData.ScoringResults.allfileIDs = [combinedData.ScoringResults.allfileIDs; did1Data.ScoringResults.allfileIDs];
combinedData.ScoringResults.alllabels = [combinedData.ScoringResults.alllabels; did1Data.ScoringResults.alllabels];

% Save the updated data back into the combined file
save('JF046_Manual_ScoringResults2.0_Combined.mat', '-struct', 'combinedData');
