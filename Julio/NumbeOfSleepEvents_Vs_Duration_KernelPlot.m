clear; clc;

%% Define Parent Folders for Control and Experimental Groups
control_folder = 'I:\Combined_Cohort\Day 5';  % Change this to the actual path
experimental_folder = 'K:\GFP+Alcohol_all\DID4';  % Change this to the actual path

% List of mouse IDs
control_mice = {'JF048', 'JF051', 'JF054', 'JF067', 'JF068', 'JF072', 'JF073'}; % Control group mice
experimental_mice = {'JF037', 'JF038', 'JF039', 'JF040', 'JF044', 'JF045', 'JF046', 'JF050'}; % Experimental group mice

bin_size = 5; % Each bin represents 5 seconds
state_values = [1, 2]; % NREM = 1, REM = 2

% Initialize duration storage for both groups
nrem_control = [];
rem_control = [];
nrem_experimental = [];
rem_experimental = [];

%% Collect Data for Control and Experimental Groups
[nrem_control, rem_control] = Collect_durations_EtohVsWater(control_folder, control_mice, bin_size, state_values);
[nrem_experimental, rem_experimental] = Collect_durations_EtohVsWater(experimental_folder, experimental_mice, bin_size, state_values);

%% Kevin plot
%% Define histogram edges
edges = 0:10:max([nrem_control; nrem_experimental]); % Adjusted for sleep durations
%edges = 0:10:100 ; % Adjusted for sleep durations

%% Create Figure for Smoothed Kernel Density Plot
figure;
hold on;

% Compute smoothed curves
[curve_control] = SmoothHistogramBins(nrem_control, edges);
[curve_experimental] = SmoothHistogramBins(nrem_experimental, edges);

% Plot control distribution
before = findall(gca);
fnplt(curve_control);
added = setdiff(findall(gca), before);
set(added, 'Color', 'b', 'LineWidth', 2); % Control in blue

% Plot experimental distribution
before = findall(gca);
fnplt(curve_experimental);
added = setdiff(findall(gca), before);
set(added, 'Color', 'r', 'LineWidth', 2); % Experimental in red

% Labels and Formatting
title('Smoothed Distribution of NREM Event Durations');
xlabel('NREM Event Duration (seconds)');
ylabel('Probability');
xlim([0, max([nrem_control; nrem_experimental])]);
ylim([0, 0.5]); % Adjusted based on your original ylim
axis square;
set(gca, 'box', 'off');
grid on;
hold off;

%% REM kevins
figure;
hold on;

% Compute smoothed curves
[curve_control] = SmoothHistogramBins(rem_control, edges);
[curve_experimental] = SmoothHistogramBins(rem_experimental, edges);

% Plot control distribution
before = findall(gca);
fnplt(curve_control);
added = setdiff(findall(gca), before);
set(added, 'Color', 'b', 'LineWidth', 2); % Control in blue

% Plot experimental distribution
before = findall(gca);
fnplt(curve_experimental);
added = setdiff(findall(gca), before);
set(added, 'Color', 'r', 'LineWidth', 2); % Experimental in red

% Labels and Formatting
title('Smoothed Distribution of REM Event Durations');
xlabel('NREM Event Duration (seconds)');
ylabel('Probability');
xlim([0, max([rem_control; rem_experimental])]);
ylim([0, 0.5]); % Adjusted based on your original ylim
axis square;
set(gca, 'box', 'off');
grid on;
hold off;


%% Plot NREM Histogram (Control vs Experimental)
figure;
hold on;
bins = 0:10:max([nrem_control; nrem_experimental]); % Define common bins
histogram(nrem_control, bins, 'FaceColor', 'b', 'FaceAlpha', 1, 'EdgeColor', 'none'); % Control in blue
histogram(nrem_experimental, bins, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Experimental in green
xlabel('NREM Event Duration (seconds)');
ylabel('Frequency');
title('Histogram of NREM Durations: Control vs Experimental');
legend({'Control', 'Experimental'});
grid on;
hold off;

%% Plot NREM KDE (Control vs Experimental)
figure;
hold on;
if ~isempty(nrem_control)
    [f_nrem_ctrl, xi_nrem_ctrl] = ksdensity(nrem_control);
    plot(xi_nrem_ctrl, f_nrem_ctrl, 'b', 'LineWidth', 2); % Control in blue
end
if ~isempty(nrem_experimental)
    [f_nrem_exp, xi_nrem_exp] = ksdensity(nrem_experimental);
    plot(xi_nrem_exp, f_nrem_exp, 'g', 'LineWidth', 2); % Experimental in green
end
xlabel('NREM Event Duration (seconds)');
ylabel('Density');
title('Kernel Density Estimation of NREM Durations: Control vs Experimental');
legend({'Control', 'Experimental'});
grid on;
hold off;

%% Plot REM Histogram (Control vs Experimental)
figure;
hold on;
bins = 0:10:max([rem_control; rem_experimental]); % Define common bins
histogram(rem_control, bins, 'FaceColor', 'r', 'FaceAlpha', 1, 'EdgeColor', 'none'); % Control in red
histogram(rem_experimental, bins, 'FaceColor', 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Experimental in magenta
xlabel('REM Event Duration (seconds)');
ylabel('Frequency');
title('Histogram of REM Durations: Control vs Experimental');
legend({'Control', 'Experimental'});
grid on;
hold off;

%% Plot REM KDE (Control vs Experimental)
figure;
hold on;
if ~isempty(rem_control)
    [f_rem_ctrl, xi_rem_ctrl] = ksdensity(rem_control);
    plot(xi_rem_ctrl, f_rem_ctrl, 'r', 'LineWidth', 2); % Control in red
end
if ~isempty(rem_experimental)
    [f_rem_exp, xi_rem_exp] = ksdensity(rem_experimental);
    plot(xi_rem_exp, f_rem_exp, 'm', 'LineWidth', 2); % Experimental in magenta
end
xlabel('REM Event Duration (seconds)');
ylabel('Density');
title('Kernel Density Estimation of REM Durations: Control vs Experimental');
legend({'Control', 'Experimental'});
grid on;
hold off;



