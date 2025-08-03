function [] = FigS01_WhiskerProbabilityBefandAftPuff(rootFolder,saveState,delim)

path = [rootFolder delim 'Results_SST_project'];
cd(path)
resultsStruct = 'Results_WhiskerProbability';
load(resultsStruct);

% Set up group, day
groups = {'Water', 'Alcohol'};
days = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};
TimeVector = Results_WhiskerProbability.Alcohol.Day1.JF037.TimeVector;

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through animals
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_WhiskerProbability.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_WhiskerProbability.(currentGroup).(currentDay));

        if isempty(FPanimalIDs)
            continue
        end

        % gather all animals' raw traces for this day
        allWhiskerProbability = [];

        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa};
            if ~isfield(Results_WhiskerProbability.(currentGroup).(currentDay), animalID)
                continue
            end
            S = Results_WhiskerProbability.(currentGroup).(currentDay).(animalID).RawWhiskerProbability;
            allWhiskerProbability = [allWhiskerProbability; S];   %#ok<AGROW>
        end

        % compute and store stats
        N = size(allWhiskerProbability,1);
        data.(currentGroup).(currentDay).WhiskerN  = N;
        data.(currentGroup).(currentDay).WhiskerProbabilityMean   = mean(allWhiskerProbability,   1);
        data.(currentGroup).(currentDay).WhiskerProbabilitySEM   = std(allWhiskerProbability,0,1) / sqrt(N);
    end
end

%Plot
Day1_color = [0 0 0];
Day5_color = [0.38, 0.56, 0.74];

% GFP day 1
x = TimeVector (:)';
y = data.Alcohol.Day1.WhiskerProbabilityMean(:)';
e = data.Alcohol.Day1.WhiskerProbabilitySEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', Day1_color, 'LineWidth', 2);

% GFP day 5 
x2 = TimeVector (:)';
y2 = data.Alcohol.Day5.WhiskerProbabilityMean(:)';
e2 = data.Alcohol.Day5.WhiskerProbabilitySEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', Day5_color, 'LineWidth', 2);

title('Probability of whisking before and after sensory stimulation.')
ylabel('Probability')
xlabel('Time (s)')
legend([p1 p2],'BL','DID4')
axis square; 
axis tight;
xlim([-15 15]);


