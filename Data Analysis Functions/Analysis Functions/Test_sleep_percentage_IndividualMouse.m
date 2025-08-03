zap;
% Define the animal ID
animalID = 'JF046';

% Load the corresponding manual scoring results
loadFile = strcat(animalID, '_Manual_ScoringResults.mat');
if isfile(loadFile)
    load(loadFile);
    
    % Calculate the total number of scores
    numberOfScores = length(ScoringResults.alllabels);
    
    % Calculate percentages for different sleep states
    iAwakePerc = round((sum(strcmp(ScoringResults.alllabels, 'Not Sleep')) / numberOfScores) * 100, 1);
    iNremPerc = round((sum(strcmp(ScoringResults.alllabels, 'NREM Sleep')) / numberOfScores) * 100, 1);
    iRemPerc = round((sum(strcmp(ScoringResults.alllabels, 'REM Sleep')) / numberOfScores) * 100, 1);
    
else
    warning('File not found: %s', loadFile);
end

% Combine the percentage data into a matrix where each percentage is a separate column
data = [iAwakePerc, iNremPerc, iRemPerc];

% Create the horizontal stacked bar chart
figure; 
b = barh(1, data, 'stacked');  % Plot a single horizontal stacked bar
title('Sleep State Distribution');
xlabel('Percentage');
set(gca, 'ytick', []); % Remove y-axis ticks since there's only one bar
xlim([0 100]); % Ensure the x-axis goes from 0 to 100%
legend({'Awake', 'NREM Sleep', 'REM Sleep'}, 'Location', 'best');

% Change colors
colors = [0.2 0.6 1; 0.1 0.9 0.1; 1 0.5 0.1];  % Define custom colors: blue, green, orange
for i = 1:length(b)
    b(i).FaceColor = colors(i, :);
end

% Add text labels inside each bar section
barWidth = cumsum(data);  % Cumulative sum to find end positions of each bar segment
textPositions = [0, barWidth(1:end-1)];  % Starting positions for the text
textMiddle = (barWidth - [0, barWidth(1:end-1)]) / 2 + textPositions;  % Middle of each segment
labels = {num2str(iAwakePerc), num2str(iNremPerc), num2str(iRemPerc)};  % Text labels as strings
for i = 1:length(data)
    text(textMiddle(i), 1, labels{i}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
end



