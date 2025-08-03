zap;
load ("AnalysisResults.mat");

% List of mice identifiers
mice = {'JF037', 'JF038', 'JF039', 'JF048', 'JF050'};

% Initialize structures to hold NREM and REM data
NREMforCohen = struct();
REMforCohen = struct();

% Loop through each mouse and extract the data
for i = 1:length(mice)
    mouse_id = mice{i};
    
    % Access the GRAB_NE data for each mouse
    SST_NREMtoREM = AnalysisResults.(mouse_id).Transitions.NREMtoREM.GRAB_NE;
    
    % Extract the first 900 cells for NREM
    SST_NREM.(mouse_id) = SST_NREMtoREM(1:900);
    
    % Extract cells 901-1800 for REM
    SST_REM.(mouse_id) = SST_NREMtoREM(901:1800);
end

% Initialize a structure to hold Cohen's d values
cohen_d = struct();

% Loop through each mouse and calculate Cohen's d
for i = 1:length(mice)
    mouse_id = mice{i};
    
    % Extract NREM and REM data for the current mouse
    M1 = SST_NREM.(mouse_id);
    M2 = SST_REM.(mouse_id);
    
    % Calculate the means of NREM and REM
    mean_SST_NREM = mean(M1);
    mean_SST_REM = mean(M2);
    
    % Calculate the variances
    var_SST_NREM = var(M1);
    var_SST_REM = var(M2);
    
    % Number of observations
    n1_SST_NREM = length(M1);
    n2_SST_REM = length(M2);
    
    % Calculate pooled standard deviation
    pooled_std = sqrt(((n1_SST_NREM-1)*var_SST_NREM + (n2_SST_REM-1)*var_SST_REM) / (n1_SST_NREM + n2_SST_REM - 2));
    
    % Calculate Cohen's d
    cohen_d.(mouse_id) = (mean_SST_NREM - mean_SST_REM) / pooled_std;
end





