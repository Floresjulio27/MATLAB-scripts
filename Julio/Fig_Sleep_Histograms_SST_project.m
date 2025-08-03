function [] = Fig_Sleep_Histograms2.0_SST_project(rootFolder,saveState,delim)
%________________________________________________________________________________________________________________________
% Written by Julio Flores-Cuadra
% The Pennsylvania State University, Dept. of Biology
%
% Acknowledgments. This code was written following scripts from Dr. Kevin
% Turner and MD Shakhawat Hossain
%________________________________________________________________________________________________________________________

%% set-up and process data
path = [rootFolder delim 'Results_SST_project']; % Path to the folder containing Results_SleepData.mat
cd(path)
resultsStruct = 'Results_SleepData';
load(resultsStruct)

groups = {'Water', 'Alcohol'};
days = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};
bin_labels = Results_SleepData.Water.Bin_Labels; % These are the same for every group/day
nBins = numel(bin_labels);
nDays = numel(days);
nGroups = numel(groups);


% Preallocate summary arrays:
%   NREM_summary(:,:,g)  is an (nBins × nDays) matrix for group g’s NREM_counts
%   REM_summary(:,:,g)   is an (nBins × nDays) matrix for group g’s REM_counts
NREM_summary = zeros(nBins, nDays, nGroups);
REM_summary  = zeros(nBins, nDays, nGroups);

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{ee};  
    
    % Loop through days
    for xx = 1:length(days)
        currentDay = days{xx};   
        
        % Check if this group actually has data for this day
        if ~isfield(Results_SleepData.(currentGroup), currentDay)
            continue
        end
        
        FPanimalIDs = fieldnames(Results_SleepData.(currentGroup).(currentDay));
        if isempty(FPanimalIDs)
            continue
        end
        
        % Initialize per-(group,day) accumulators
        total_nrem = zeros(1, nBins);
        total_rem  = zeros(1, nBins);
        
        % Loop through animals for this (group,day)
        for aa = 1:numel(FPanimalIDs)
            animalID = FPanimalIDs{aa};
            
            if isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'NREM_counts')
                total_nrem = total_nrem + Results_SleepData.(currentGroup).(currentDay).(animalID).NREM_counts;
            end
            if isfield(Results_SleepData.(currentGroup).(currentDay).(animalID), 'REM_counts')
                total_rem  = total_rem  + Results_SleepData.(currentGroup).(currentDay).(animalID).REM_counts;
            end
        end
        
        % Store the sums into array 3D array
        NREM_summary(:, xx, ee) = total_nrem(:); %nbins x days x groups = 6x6x2
        REM_summary(:,  xx, ee) = total_rem(:);
    end
end

%% Plot NREM water group across days 
waterIdx = strcmp(groups, 'Water');  % logical index
waterIdx = find(waterIdx); % 1 should be for water

% Extract the matrix (nBins × nDays) for Water → NREM
allNREM_Water = NREM_summary(:, :, waterIdx);  % size: [nBins × nDays]

figure;
bar(allNREM_Water, 'grouped');
xticklabels(bin_labels);
xtickangle(45);
xlabel('Sleep‐Event Duration Bins');
ylabel('Number of NREM Events');
title('Water Group ▶︎ NREM Events (summed over all animals)');

% Legend: one entry per day
legend(days, 'Location', 'northeastoutside');

%% Plot REM water group across days
allREM_Water = REM_summary(:, :, waterIdx);  % size: [nBins × nDays]

figure;
bar(allREM_Water, 'grouped');
xticklabels(bin_labels);
xtickangle(45);
xlabel('Sleep‐Event Duration Bins');
ylabel('Number of REM Events');
title('Water Group ▶︎ REM Events (summed over all animals)');

legend(days, 'Location', 'northeastoutside');
grid on;

%% Water vs Alcohol across days (NREM
waterIdx   = find(strcmp(groups, 'Water'));
alcoholIdx = find(strcmp(groups, 'Alcohol'));
if isempty(alcoholIdx)
    error('No "Alcohol" group found in Results_SleepData.');
end

commonDays = intersect(fieldnames(Results_SleepData.Water), fieldnames(Results_SleepData.Alcohol));
nCommon    = numel(commonDays);
if nCommon == 0
    error('No overlapping day‐names between Water and Alcohol.');
end

figure('Name','NREM: Water vs Alcohol','NumberTitle','off');

for d = 1:nCommon
    dayName = commonDays{d};
    
    % Extract Water's NREM vector for this day (nBins×1)
    waterNREMsum = NREM_summary(:, strcmp(days, dayName), waterIdx);
    % Extract Alcohol's NREM vector for this day
    alcoholNREMsum = NREM_summary(:, strcmp(days, dayName), alcoholIdx);
    
    % In case one group has the field but no animals (sum would be all zeros)
    if isempty(waterNREMsum)
        waterNREMsum = zeros(nBins, 1);
    end
    if isempty(alcoholNREMsum)
        alcoholNREMsum = zeros(nBins, 1);
    end
    
    % Create grouped bar: columns = [Water, Alcohol], rows = bins
    subplot(nCommon, 1, d);
    barData = [waterNREMsum, alcoholNREMsum];  % size [nBins × 2]
    b = bar(barData, 'grouped');
    b(1).FaceColor = [0 0.5 1];   % Water = bluish
    b(2).FaceColor = [1 0 0];     % Alcohol = reddish
    
    xticklabels(bin_labels);
    xtickangle(45);
    ylabel('Count');
    title(sprintf('%s ▶︎ NREM: Water vs Alcohol', dayName));
    
    if d == 1
        legend({'Water','Alcohol'}, 'Location','northwest');
    end
    grid on;
    
    if d == nCommon
        xlabel('Sleep‐Event Duration Bins');
    end
end

sgtitle('NREM Event Distribution: Water vs Alcohol (per day)');

%% Water vs alcohol across days (REM)
figure('Name','REM: Water vs Alcohol','NumberTitle','off');

for d = 1:nCommon
    dayName = commonDays{d};
    
    % Extract Water's REM vector for this day
    waterREMsum = REM_summary(:, strcmp(days, dayName), waterIdx);
    % Extract Alcohol's REM vector for this day
    alcoholREMsum = REM_summary(:, strcmp(days, dayName), alcoholIdx);
    
    if isempty(waterREMsum)
        waterREMsum = zeros(nBins, 1);
    end
    if isempty(alcoholREMsum)
        alcoholREMsum = zeros(nBins, 1);
    end
    
    subplot(nCommon, 1, d);
    barData = [waterREMsum, alcoholREMsum];
    b = bar(barData, 'grouped');
    b(1).FaceColor = [0 0.5 1];   % Water = blue
    b(2).FaceColor = [1 0 0];     % Alcohol = red
    
    xticklabels(bin_labels);
    xtickangle(45);
    ylabel('Count');
    title(sprintf('%s ▶︎ REM: Water vs Alcohol', dayName));
    
    if d == 1
        legend({'Water','Alcohol'}, 'Location','northwest');
    end
    grid on;
    
    if d == nCommon
        xlabel('Sleep‐Event Duration Bins');
    end
end

sgtitle('REM Event Distribution: Water vs Alcohol (per day)');
