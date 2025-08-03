function [] = Fig_Total_Sleep_Bouts_Histfit_SST_project(rootFolder,saveState,delim)
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
days   = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};
% Loop through groups first, then days & animals,
% so that for each (group, day) you collect ALL bouts

for ee = 1:length(groups)
    currentGroup = groups{ee};
    
    for xx = 1:length(days)
        currentDay = days{xx};
        
        % Check if this group has data for currentDay
        if ~isfield(Results_SleepData.(currentGroup), currentDay)
            continue
        end
        
        FPanimalIDs = fieldnames(Results_SleepData.(currentGroup).(currentDay));
        if isempty(FPanimalIDs)
            continue
        end

        % Pre‐allocate per (group, day):
        allNREM_bouts = [];
        allREM_bouts  = [];

        % Loop through every animal in this (group, day)
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa};
            animalData = Results_SleepData.(currentGroup).(currentDay).(animalID);

            if isfield(animalData, 'NREM_Bouts')
                NREMpath = animalData.NREM_Bouts(:,3);
                allNREM_bouts = [allNREM_bouts; NREMpath];   %#ok<AGROW>
            end

            if isfield(animalData, 'REM_Bouts')
                REMpath = animalData.REM_Bouts(:,3);
                allREM_bouts = [allREM_bouts; REMpath];      %#ok<AGROW>
            end
        end
        
        % At this point, allNREM_bouts and allREM_bouts contain
        % every bout (all animals) for (currentGroup, currentDay).
        % Store them in your data struct, keeping groups separate:
        data.(currentGroup).(currentDay).NREM_bouts_All = allNREM_bouts;
        data.(currentGroup).(currentDay).REM_bouts_All  = allREM_bouts;
        
        % compute & store stats to calculate bouts lenght in the future:
        N = size(allNREM_bouts, 1);
        data.(currentGroup).(currentDay).NREM_bouts_N_Exp = N;
        data.(currentGroup).(currentDay).NREM_bouts_Mean  = mean(allNREM_bouts, 1);
        data.(currentGroup).(currentDay).NREM_bouts_SEM   = std(allNREM_bouts, 0, 1) / sqrt(N);

        N = size(allREM_bouts, 1);
        data.(currentGroup).(currentDay).REM_bouts_N_Exp  = N;
        data.(currentGroup).(currentDay).REM_bouts_Mean   = mean(allREM_bouts, 1);
        data.(currentGroup).(currentDay).REM_bouts_SEM    = std(allREM_bouts, 0, 1) / sqrt(N);
    end
end


%% Prep data for plotting
% Combine all days NREM bouts into one vector per group so we can do all
% water vs all alcohol
for ee = 1:length(groups)
    currentGroup = groups{ee};
    allNREM = [];
    allREM = [];
    for xx = 1:length(days)
        currentDay = days{xx};
        if isfield(data.(currentGroup), currentDay)
            allNREM = [allNREM; data.(currentGroup).(currentDay).NREM_bouts_All];
            allREM = [allREM; data.(currentGroup).(currentDay).REM_bouts_All];
        end
    end
    
    % store the combined vector
    data.(currentGroup).AllDays_NREM = allNREM;
    data.(currentGroup).AllDays_REM = allREM;
end

%% Plot Histfit water vs alcohol NREM
Water_color   = [0 0 0];
Alcohol_color = [0.8667, 0.1098, 0.4667];

ax1 = subplot(1,2,1)
hold on;

% Water
hW = histfit( data.Water.AllDays_NREM, 20, 'kernel' );
set( hW(1), 'FaceColor', Water_color,   'FaceAlpha', 0.3, 'EdgeColor', 'none' );
set( hW(2), 'Color',     Water_color,   'LineWidth', 2 );

% Alcohol
hA = histfit( data.Alcohol.AllDays_NREM, 20, 'kernel' );
set( hA(1), 'FaceColor', Alcohol_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none' );
set( hA(2), 'Color',     Alcohol_color, 'LineWidth', 2 );

xlabel('NREM bout duration (s)');
xlim([0 250])
ylabel('NREM Count');
legend('Water histogram','Water fit','Alcohol histogram','Alcohol fit');
title('NREM bout distributions (all days combined)');
%% Plot Histfit water vs alcohol REM
Water_color   = [0 0 0];
Alcohol_color = [0.8667, 0.1098, 0.4667];

ax2 = subplot(1,2,2)
hold on;

% Water
hW = histfit( data.Water.AllDays_REM, 20, 'kernel' );
set( hW(1), 'FaceColor', Water_color,   'FaceAlpha', 0.3, 'EdgeColor', 'none' );
set( hW(2), 'Color',     Water_color,   'LineWidth', 2 );

% Alcohol
hA = histfit( data.Alcohol.AllDays_REM, 20, 'kernel' );
set( hA(1), 'FaceColor', Alcohol_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none' );
set( hA(2), 'Color',     Alcohol_color, 'LineWidth', 2 );

xlabel('REM bout duration (s)');
xlim([0 300])
ylabel('Count');
legend('Water histogram','Water fit','Alcohol histogram','Alcohol fit');
title('REM bout distributions (all days combined)');
hold off;

%% plot with bins NREM 

% figure;
% binWidth = 20;
% maxEdge  = 200;
% edges    = [0:binWidth:maxEdge, Inf];
% 
% % Plot Water
% hW = histogram(data.Water.AllDays_NREM, edges, ...
%                'Normalization','count', ...
%                'FaceColor',Water_color,'FaceAlpha',0.3,'EdgeColor','none');
% hold on
% % Plot Alcohol
% hA = histogram(data.Alcohol.AllDays_NREM, edges, ...
%                'Normalization','count', ...
%                'FaceColor',Alcohol_color,'FaceAlpha',0.3,'EdgeColor','none');
% 
% % Build labels for each bin:
% nBins = numel(edges)-1;
% labels = cell(1,nBins);
% for i = 1:nBins
%     if isinf(edges(i+1))
%         labels{i} = '>200';
%     else
%         labels{i} = sprintf('%d–%d', edges(i), edges(i+1));
%     end
% end
% 
% % Place ticks at the center of each bin:
% centers = (edges(1:end-1) + edges(2:end))/2;
% xticks(centers);
% xticklabels(labels);
% 
% xlabel('NREM bout duration (s)');
% ylabel('Count');
% legend('Water','Alcohol');
% title('Bouts grouped with final “>200 s” bin');

end 
