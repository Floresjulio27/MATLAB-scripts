function [] = Fig_Alcohol_MeanGFP_relationship(rootFolder,saveState,delim)
%________________________________________________________________________________________________________________________
% Written by Julio Flores-Cuadra
% The Pennsylvania State University, Dept. of Biology
%
%________________________________________________________________________________________________________________________

%% set-up and process data
path = [rootFolder delim 'Results_SST_project']; % Path to the folder containing Results
cd(path)
resultsStruct_GFP = 'Results_MeanGFP';
load(resultsStruct_GFP)
resultsStruct_alcohol = "Results_AlcoholDrinking";
load(resultsStruct_alcohol)


groups = {'Alcohol'};
days = {'Day5'};
behavior = {'Rest', 'Whisk', 'Stim', 'NREM', 'REM'};
hemisphere = {'MeanNE', 'MeanACh'}; % NE is PFC and ACh is S1BF

% %% Create summary table for total etoh consumed
% 
% % Get list of AnimalIDs from TotalEtoh sub‑struct
% animalList = fieldnames(AlcoholDrinking.TotalEtoh);
% nAnimals   = numel(animalList);
% 
% Total_Ethanol = nan(nAnimals,1);
% 
% for i = 1:nAnimals
%     AnimalID = animalList{i};
%     % total ethanol consumed across all DID cycles
%     Total_Ethanol(i) = AlcoholDrinking.TotalEtoh.(AnimalID).TotalEtohConsumed;
% end
% 
% Alcohol_Table = table(animalList, Total_Ethanol,'VariableNames', {'AnimalID','Total_Ethanol'} );

%% Create summary table for etoh consumen on binge day 4
Animals = {'JF037','JF038','JF039','JF040','JF082','JF083','JF084'};
nAnimals = numel(Animals);
% Get list of AnimalIDs from TotalEtoh sub‑struct
Total_Ethanol_DID4 = nan(nAnimals,1);

for i = 1:nAnimals
    AnimalID = Animals{i};
    % total ethanol consumed across all DID cycles
    Total_Ethanol_DID4(i) = AlcoholDrinking.DID4.(AnimalID).BingeDay;
end

Alcohol_Table = table(Animals(:), Total_Ethanol_DID4,'VariableNames', {'AnimalID','Total_Ethanol'} );

%% Extract mean PFC and S1BF data from DID4 and create table 
behaviors = {'Rest','NREM','REM','Stim','Whisk'};
groups    = {'Alcohol'};

PFC_data   = nan(nAnimals, numel(behaviors));
S1BF_data  = nan(nAnimals, numel(behaviors));
%% you need to calculate the mean of all the events first to have one measurement per mouse. 
for g = 1:length(groups)
    grp = groups{g};

    for b = 1:length(behaviors)
        currentBehavior = behaviors{b};

        for i = 1:nAnimals
            AnimalID = Alcohol_Table.AnimalID{i};

            if ~isfield(Results_MeanGFP.(grp).Day5.(AnimalID).MeanGFP, currentBehavior)
                continue
            end

            allAChGFP = [];
            allNEGFP  = [];

            S = Results_MeanGFP.(grp).Day5.(AnimalID).MeanGFP.(currentBehavior).GFP;

            if isfield(S, 'MeanNE'), allAChGFP = [allAChGFP; S.MeanACh]; end
            if isfield(S, 'MeanACh'),  allNEGFP  = [allNEGFP;  S.MeanNE];  end

            % store group stats (mean ± SEM across animals)
            N = size(allAChGFP,1);
            data.(AnimalID).(currentBehavior).N = N;
            data.(AnimalID).(currentBehavior).TotalMeanS1BF = mean(allAChGFP, 1);
            data.(AnimalID).(currentBehavior).TotalMeanPFC = mean(allNEGFP, 1);

            %Store the GFP data acording to hemisphere/ animalID and
            %behavior
            PFC_data(i,b)  =  data.(AnimalID).(currentBehavior).TotalMeanPFC;
            S1BF_data(i,b) = data.(AnimalID).(currentBehavior).TotalMeanS1BF;
        end
    end
end

% Build a table of GFP data
PFC_table   = array2table(PFC_data,  'VariableNames', strcat('PFC_',  behaviors));
S1BF_table  = array2table(S1BF_data, 'VariableNames', strcat('S1BF_', behaviors));
MeanGFP_table  = [table(Alcohol_Table.AnimalID,'VariableNames',{'AnimalID'}) PFC_table S1BF_table];

%% Merge alcohol talble and GFP table 
merged = innerjoin(Alcohol_Table, MeanGFP_table, 'Keys','AnimalID');

%% Calculate correlation. Do not take into account NaN values

disp('Correlations vs Total_Ethanol');
for b = 1:numel(behaviors)
    currentBehavior = behaviors{b};

    % Extract values from the table
    ethanolVals = merged.Total_Ethanol;
    pfcVals = merged{:,sprintf('PFC_%s',currentBehavior)};
    s1bfVals = merged{:,sprintf('S1BF_%s',currentBehavior)};

    % Remove NaNs for PFC
    validPFC = ~isnan(ethanolVals) & ~isnan(pfcVals);
    if any(validPFC)
        [rPFC, pPFC] = corr(ethanolVals(validPFC), pfcVals(validPFC), 'Type', 'Pearson');
    else
        rPFC = NaN; pPFC = NaN;
    end

    % Remove NaNs for S1BF
    validS1BF = ~isnan(ethanolVals) & ~isnan(s1bfVals);
    if any(validS1BF)
        [rS1BF, pS1BF] = corr(ethanolVals(validS1BF), s1bfVals(validS1BF), 'Type', 'Pearson');
    else
        rS1BF = NaN; pS1BF = NaN;
    end

    fprintf('%s:  PFC r=%.3f (p=%.3f),  S1BF r=%.3f (p=%.3f)\n', currentBehavior, rPFC,pPFC, rS1BF,pS1BF);
end


%% Scatter plot for visualizating relationship

%PFC data 
fig1 = figure;
subplot(2,3,1)
scatter(merged.Total_Ethanol, merged.PFC_Rest, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST PFC (DID4) during Rest');
title('Total Alcohol conumption vs Mean GFP SST PFC during Rest ');
lsline;

subplot(2,3,2)
scatter(merged.Total_Ethanol, merged.PFC_Whisk, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST PFC (DID4) during Whisk');
title('Total Alcohol conumption vs Mean GFP SST PFC during Whisk ');
lsline;

subplot(2,3,3)
scatter(merged.Total_Ethanol, merged.PFC_Stim, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST PFC (DID4) during Stim');
title('Total Alcohol conumption vs Mean GFP SST PFC during Stim ');
lsline;

subplot(2,3,4)
scatter(merged.Total_Ethanol, merged.PFC_NREM, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST PFC (DID4) during NREM');
title('Total Alcohol conumption vs Mean GFP SST PFC during NREM ');
lsline;

subplot(2,3,5)
scatter(merged.Total_Ethanol, merged.PFC_REM, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST PFC (DID4) during REM');
title('Total Alcohol conumption vs Mean GFP SST PFC during REM ');
xlim([0 12])
lsline;

%S1BF data 
fig2 = figure;
subplot(2,3,1)
scatter(merged.Total_Ethanol, merged.S1BF_Rest, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST S1BF (DID4) during Rest');
title('Total Alcohol conumption vs Mean GFP SST S1BF during Rest ');
lsline;

subplot(2,3,2)
scatter(merged.Total_Ethanol, merged.S1BF_Whisk, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST S1BF (DID4) during Whisk');
title('Total Alcohol conumption vs Mean GFP SST S1BF during Whisk ');
lsline;

subplot(2,3,3)
scatter(merged.Total_Ethanol, merged.S1BF_Stim, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST S1BF (DID4) during Stim');
title('Total Alcohol conumption vs Mean GFP SST S1BF during Stim ');
lsline;

subplot(2,3,4)
scatter(merged.Total_Ethanol, merged.S1BF_NREM, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST S1BF (DID4) during NREM');
title('Total Alcohol conumption vs Mean GFP SST S1BF during NREM ');
lsline;

subplot(2,3,5)
scatter(merged.Total_Ethanol, merged.S1BF_REM, 'filled');
xlabel('Total Ethanol conusmption (g/kg)');
ylabel('Mean GFP (%) SST S1BF (DID4) during REM');
title('Total Alcohol conumption vs Mean GFP SST S1BF during REM ');
lsline;

%% Save figure
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figures Water vs Alcohol' delim 'Correlations'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(fig1, fullfile(dirpath, 'Total Alcohol conumption vs Mean GFP SST PFC - All Behaviors'));
    savefig(fig2, fullfile(dirpath, 'Total Alcohol conumption vs Mean GFP SST S1BF - All Behaviors'));
end
