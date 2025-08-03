% Clear workspace, close figures, and clear command window
clearvars; close all; clc;

% Define the list of animal IDs
animalIDs = {'JF048', 'JF050', 'JF051', 'JF054'};
%animalIDs = {'JF041', 'JF043', 'JF048', 'JF050', 'JF051', 'JF054'};
% Preallocate arrays to store percentages
numAnimals = length(animalIDs);
awakePercents = NaN(1, numAnimals);
nremPercents = NaN(1, numAnimals);
remPercents = NaN(1, numAnimals);

% Loop over each animal to calculate sleep state percentages
for idx = 1:numAnimals
    animalID = animalIDs{idx};
    loadFile = strcat(animalID, '_Manual_ScoringResults.mat');
    
    if isfile(loadFile)
        load(loadFile);
        % Calculate percentages for sleep states
        [awakePercents(idx), nremPercents(idx), remPercents(idx)] = ...
            calculateSleepPercentages(ScoringResults.alllabels);
    else
        warning('File not found: %s', loadFile);
    end
end

% Remove NaN entries (files not found)
validIndices = ~isnan(awakePercents);
awakePercents = awakePercents(validIndices);
nremPercents = nremPercents(validIndices);
remPercents = remPercents(validIndices);

% Check if data was loaded for at least one animal
if isempty(awakePercents)
    error('No data loaded. Please check the animal IDs and corresponding files.');
end

% Compute average percentages across all animals
avgData = [mean(awakePercents), mean(nremPercents), mean(remPercents)];

% Verify total percentage sums to approximately 100%
totalPerc = sum(avgData);
if abs(totalPerc - 100) > 0.1
    warning('Total percentage does not sum to 100%% (current total: %.1f%%).', totalPerc);
end

% Create figure with subplots for the pie chart and ternary plot
figure;

% Pie chart of average sleep state distribution
subplot(2, 4, 1);  % Positioning for 2 rows, 4 columns, at 1st position
p = pie(avgData);
pText = findobj(p, 'Type', 'text');
percentValues = get(pText, 'String');
labels = {'Awake: ', 'NREM: ', 'REM: '};
newLabels = strcat(labels, percentValues);
[pText.String] = deal(newLabels{:});
title({'Sleep scoring label probability', 'Mean animal sleep scoring labels'});

% Ternary plot of individual animal data
ax2 = subplot(2, 4, 2);  % Positioning for ternary plot in 2nd slot
terplot();
% Plot individual animal sleep state percentages in ternary plot
[hd] = ternaryc(awakePercents/100, nremPercents/100, remPercents/100);
hlabels = terlabel('rfc-Awake', 'rfc-NREM', 'rfc-REM');
title({'[2c] Ternary plot of individual animals', ' ', ' '});

% Save the figure as a .png and .fig file
saveas(gcf, 'AverageSleepStateDistribution.png');
savefig(gcf, 'AverageSleepStateDistribution.fig');

% Function to calculate sleep percentages
function [awakePerc, nremPerc, remPerc] = calculateSleepPercentages(labels)
    numberOfScores = length(labels);
    awakePerc = round((sum(strcmp(labels, 'Not Sleep')) / numberOfScores) * 100, 1);
    nremPerc = round((sum(strcmp(labels, 'NREM Sleep')) / numberOfScores) * 100, 1);
    remPerc = round((sum(strcmp(labels, 'REM Sleep')) / numberOfScores) * 100, 1);
end
