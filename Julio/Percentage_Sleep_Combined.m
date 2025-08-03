% Clear variables, close figures, and clear command window
clearvars; close all; clc;

% Define the list of animal IDs
%animalIDs = {'JF037', 'JF041', 'JF043', 'JF048', 'JF050'};
%animalIDs = {'JF037', 'JF039', 'JF040', 'JF044', 'JF046'};
animalIDs = {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF044', 'JF045', 'JF046', 'JF081', 'JF082', 'JF083', 'JF084'};
%animalIDs = {'JF048', 'JF050'}; 
%animalIDs = {'JF041', 'JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'};
%animalIDs = {'JF041','JF043', 'JF048', 'JF050','JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'};
%animalIDs = {'JF048', 'JF051' 'JF054'};
%animalIDs = {'JF048', 'JF050', 'JF051', 'JF054'};
%animalIDs = {'JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'};
%animalIDs = {'JF041'};

% Preallocate arrays to store percentages
numAnimals = length(animalIDs);
awakePercents = zeros(1, numAnimals);
nremPercents = zeros(1, numAnimals);
remPercents = zeros(1, numAnimals);

% Loop over each animal
for idx = 1:numAnimals
    animalID = animalIDs{idx};
    
    % Load the corresponding manual scoring results
    %loadFile = strcat(animalID, '_ScoringResults_Extracted_Nofirst.mat');
    loadFile = strcat(animalID, '_Manual_ScoringResults.mat');
    if isfile(loadFile)
        load(loadFile);
        
        % Calculate percentages for different sleep states
        %[awakePercents(idx), nremPercents(idx), remPercents(idx)] = calculateSleepPercentages(ScoringResults_Extracted.alllabels);
        [awakePercents(idx), nremPercents(idx), remPercents(idx)] = calculateSleepPercentages(ScoringResults.alllabels);
    else
        warning('File not found: %s', loadFile);
        awakePercents(idx) = NaN;
        nremPercents(idx) = NaN;
        remPercents(idx) = NaN;
    end
end

% Remove entries with NaN values
validIndices = ~isnan(awakePercents);
awakePercents = awakePercents(validIndices);
nremPercents = nremPercents(validIndices);
remPercents = remPercents(validIndices);

% Check if data was loaded for at least one animal
if isempty(awakePercents)
    error('No data loaded. Please check the animal IDs and corresponding files.');
end

% Compute the average percentages across all animals
avgAwakePerc = mean(awakePercents);
avgNremPerc = mean(nremPercents);
avgRemPerc = mean(remPercents);

% Check if total percentage sums to approximately 100%
totalPerc = avgAwakePerc + avgNremPerc + avgRemPerc;
if abs(totalPerc - 100) > 0.1
    warning('Total percentage does not sum to 100%% (current total: %.1f%%).', totalPerc);
end

% Combine the average percentages into a matrix
avgData = [avgAwakePerc, avgNremPerc, avgRemPerc];

% Create the horizontal stacked bar chart
% figure; 
% b = barh(1, avgData, 'stacked');  % Plot a single horizontal stacked bar
% title('Average Sleep State Distribution Across All Animals', 'FontSize', 16, 'FontWeight', 'bold');
% xlabel('Percentage', 'FontSize', 14, 'FontWeight', 'bold');
% set(gca, 'ytick', [], 'FontSize', 12, 'FontWeight', 'bold'); % Remove y-axis ticks since there's only one bar
% xlim([0 100]); % Ensure the x-axis goes from 0 to 100%
% legend({'Awake', 'NREM Sleep', 'REM Sleep'}, 'Location', 'best');

% Define specific colors for each section of the pie chart
awakeColor = [0.4, 0.4, 0.4]; % Example: light blue for 'Awake'
nremColor = [0.6392, 0.8431, 0.9137]; % Example: light green for 'NREM'
remColor = [0.9686, 0.7137, 0.3333]; % Example: light pink for 'REM'

% Plot the pie chart
figure;
p1 = pie(avgData);

% Customize colors for each slice
set(p1(1), 'FaceColor', awakeColor); % 'Awake' color
set(p1(3), 'FaceColor', nremColor);  % 'NREM' color
set(p1(5), 'FaceColor', remColor);   % 'REM' color

% Adjust text labels
pText = findobj(p1, 'Type', 'text');
percentValues = get(pText, 'String');
txt = {'Awake: '; 'NREM: '; 'REM: '};
combinedtxt = strcat(txt, percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);

% Title for the pie chart
title({'Sleep scoring label probability', 'Mean animal sleep scoring labels'});


figure (2);
terplot();
% Plot individual animal sleep state percentages in ternary plot
[hd] = ternaryc(awakePercents/100, nremPercents/100, remPercents/100);
hlabels = terlabel('Awake', 'NREM', 'REM');
title({'[2c] Ternary plot of individual animals', ' ', ' '});

% Change colors
% colors = [0.2 0.6 1; 0.1 0.9 0.1; 1 0.5 0.1];  % Define custom colors: blue, green, orange
% for i = 1:length(b)
%     b(i).FaceColor = colors(i, :);
% end

% % Add text labels inside each bar section
% cumulativeData = [0, cumsum(avgData)];
% labels = {sprintf('%.1f%%', avgAwakePerc), sprintf('%.1f%%', avgNremPerc), sprintf('%.1f%%', avgRemPerc)};
% for i = 1:length(avgData)
%     textPos = (cumulativeData(i) + cumulativeData(i+1)) / 2;
%     text(textPos, 1, labels{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
% end

% Save the figure
saveas(gcf, 'AverageSleepStateDistribution.png');
savefig(gcf, 'AverageSleepStateDistribution.fig');

% Function to calculate sleep percentages
function [awakePerc, nremPerc, remPerc] = calculateSleepPercentages(labels)
    numberOfScores = length(labels);
    awakePerc = round((sum(strcmp(labels, 'Not Sleep')) / numberOfScores) * 100, 1);
    nremPerc = round((sum(strcmp(labels, 'NREM Sleep')) / numberOfScores) * 100, 1);
    remPerc = round((sum(strcmp(labels, 'REM Sleep')) / numberOfScores) * 100, 1);
end

