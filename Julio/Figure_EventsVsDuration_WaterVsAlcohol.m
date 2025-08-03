zap;
%% Load Data
data = load("NREM_counts_Day5_Alcohol.mat");
nrem_counts_alcohol = data.nrem_counts(:); % Ensure column vector

data = load("REM_counts_Day5_Alcohol.mat");
rem_counts_alcohol = data.rem_counts(:); % Ensure column vector

data = load("NREM_counts_Day5_Water.mat");
nrem_counts_water = data.nrem_counts(:); % Ensure column vector

data = load("REM_counts_Day5_Water.mat");
rem_counts_water = data.rem_counts(:); % Ensure column vector

bin_labels = {'5-10s', '10-15s', '15-20s', '20-25s', '25-30s', '>30s'};

%% Separate Plots for NREM and REM Comparisons

% NREM Comparison (Water vs. Alcohol)
figure;
bar([nrem_counts_water nrem_counts_alcohol], 'grouped'); % Use horizontal concatenation
xticklabels(bin_labels);
xlabel('Sleep Event Duration Categories');
ylabel('Number of NREM Events');
title('NREM Sleep Event Duration: Water vs. Alcohol');
legend({'Water', 'Alcohol'}, 'Location', 'northwest');
grid on;

% REM Comparison (Water vs. Alcohol)
figure;
bar([rem_counts_water rem_counts_alcohol], 'grouped'); % Use horizontal concatenation
xticklabels(bin_labels);
xlabel('Sleep Event Duration Categories');
ylabel('Number of REM Events');
title('REM Sleep Event Duration: Water vs. Alcohol');
legend({'Water', 'Alcohol'}, 'Location', 'northwest');
grid on;

