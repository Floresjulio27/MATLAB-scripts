function [] = Fig_FP_Transition_SST_PFC_GFP_Alcohol_ANIMALS(rootFolder,saveState,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%Modified by Julio Flores-Cuadra. Re adapted from Shakhawat Hossain code
%________________________________________________________________________________________________________________________

%% set-up and process data
path = [rootFolder delim 'Results_SST_project']; %Path works way better
cd(path)
resultsStruct = 'Results_Transitions';
load(resultsStruct)

groups = {'Water', 'Alcohol'};
days = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through animals
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_Transitions.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_Transitions.(currentGroup).(currentDay));

        if isempty(FPanimalIDs)
            continue
        end

        % Loop through days
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa,1}; %it is organizing it as a column vector


            % Loop through transitions
            for cc = 1:length(transitions)
                transition = transitions{1,cc};
                if ~isfield(Results_Transitions.(currentGroup).(currentDay).(animalID).Transitions, transition)
                    continue;
                end

                %Concatenate Raw Data
                %Barrel cortex
                data.(currentGroup).(currentDay).(animalID).(transition).AChCBVRaw = Results_Transitions.(currentGroup).(currentDay).(animalID).Transitions.(transition).ACh_CBVRaw;
                data.(currentGroup).(currentDay).(animalID).(transition).AChGFPRaw = Results_Transitions.(currentGroup).(currentDay).(animalID).Transitions.(transition).GRAB_AChRaw;
                %PFC
                data.(currentGroup).(currentDay).(animalID).(transition).NECBVRaw = Results_Transitions.(currentGroup).(currentDay).(animalID).Transitions.(transition).NE_CBVRaw;
                data.(currentGroup).(currentDay).(animalID).(transition).NEGFPRaw = Results_Transitions.(currentGroup).(currentDay).(animalID).Transitions.(transition).GRAB_NERaw;
            end
        end
    end
end

% Compute Day-by-Day averages (and SEM) across days/animals. This is other way to organize it. 

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through days
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_Transitions.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_Transitions.(currentGroup).(currentDay));

        if isempty(FPanimalIDs)
            continue
        end

        % Loop through transitions
        for cc = 1:length(transitions)
            transition = transitions{1,cc};

            %% Pooled data from days

            % % gather all animals' raw traces for this day
            % allAChCBV = [];
            % allNECBV  = [];
            % allAChGFP = [];
            % allNEGFP  = [];
            % for ai = 1:length(FPanimalIDs)
            %     animalID = FPanimalIDs{ai};
            %     if ~isfield(data.(currentGroup).(currentDay).(animalID), transition)
            %         continue
            %     end
            %     S = data.(currentGroup).(currentDay).(animalID).(transition);
            %     allAChCBV = [allAChCBV; S.AChCBVRaw];   %#ok<AGROW> %concatonate vertically the data
            %     allNECBV  = [allNECBV;  S.NECBVRaw];    %#ok<AGROW>
            %     allAChGFP = [allAChGFP; S.AChGFPRaw];   %#ok<AGROW>
            %     allNEGFP  = [allNEGFP;  S.NEGFPRaw];    %#ok<AGROW>
            % end
            % 
            % % compute and store stats
            % N = size(allAChCBV,1);
            % data.(currentGroup).(currentDay).(transition).P_AChCBV_N_Exp  = N;
            % data.(currentGroup).(currentDay).(transition).P_AChCBV_Mean   = mean(allAChCBV,   1);
            % data.(currentGroup).(currentDay).(transition).P_AChCBV_SEM    = std(allAChCBV,0,1) / sqrt(N);
            % 
            % N = size(allNECBV,1);
            % data.(currentGroup).(currentDay).(transition).P_NECBV_N_Exp   = N;
            % data.(currentGroup).(currentDay).(transition).P_NECBV_Mean    = mean(allNECBV,    1);
            % data.(currentGroup).(currentDay).(transition).P_NECBV_SEM     = std(allNECBV,0,1)  / sqrt(N);
            % 
            % N = size(allAChGFP,1);
            % data.(currentGroup).(currentDay).(transition).P_AChGFP_N_Exp  = N;
            % data.(currentGroup).(currentDay).(transition).P_AChGFP_Mean   = mean(allAChGFP,   1);
            % data.(currentGroup).(currentDay).(transition).P_AChGFP_SEM    = std(allAChGFP,0,1) / sqrt(N);
            % 
            % N = size(allNEGFP,1);
            % data.(currentGroup).(currentDay).(transition).P_NEGFP_N_Exp   = N;
            % data.(currentGroup).(currentDay).(transition).P_NEGFP_Mean    = mean(allNEGFP,    1);
            % data.(currentGroup).(currentDay).(transition).P_NEGFP_SEM     = std(allNEGFP,0,1)  / sqrt(N);

            %% Average within mice and then across days

            fields = {'AChCBVRaw','NECBVRaw','AChGFPRaw','NEGFPRaw'};
            finalData = {'P_AChCBV','P_NECBV','P_AChGFP','P_NEGFP'};

            for ff = 1:length(fields)
                fieldsi = fields{ff};
                finalDatai = finalData{ff};

                 mouseMeans = [];
                 for ai = 1:length(FPanimalIDs)
                     animalID = FPanimalIDs{ai};
                     if ~isfield(data.(currentGroup).(currentDay).(animalID), transition)
                         continue
                     end
                    rawEvents = data.(currentGroup).(currentDay).(animalID).(transition).(fieldsi);
                    %average across events → 1×T
                    mouseMeans(end+1,:) = mean(rawEvents,1);  %#ok<AGROW>
                 end
                 %2) group stats across mice
                 N    = size(mouseMeans,1);
                 Tmean   = mean(mouseMeans,1);               
                 gSEM = std(mouseMeans,0,1)/sqrt(N);     

                 %3) store back
                 data.(currentGroup).(currentDay).(transition).([finalDatai '_N_Exp']) = N;
                 data.(currentGroup).(currentDay).(transition).([finalDatai '_Mean'])  = Tmean;
                 data.(currentGroup).(currentDay).(transition).([finalDatai '_SEM'])   = gSEM;
            end 
        end
    end
end
T1 = -30 + (1/30):(1/30):30;
%% Awake to NREM

%Colors
gfp_Day1_color = [0 0 0];
gfp_Day2_color = [0.38, 0.56, 0.74];
gfp_Day3_color = [0.42, 0.56, 0.14];
gfp_Day4_color = [0.71, 0.52, 0.52];
gfp_Day5_color = [0.55, 0.45, 0.70];
gfp_Day6_color = [0.42, 0.56, 0.14];

%SST Ca2+ activity in PFC/Alcohol
sgtitle ('ArousalStateTransition SST PFC GFP Alcohol pooled')
%sgtitle ('ArousalStateTransition SST PFC GFP Alcohol')

ax1 = subplot(2,2,1);

% GFP day 1
x = T1;
y = data.Alcohol.Day1.AWAKEtoNREM.P_NEGFP_Mean;
e = data.Alcohol.Day1.AWAKEtoNREM.P_NEGFP_SEM;
hold on;
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = T1;
y2 = data.Alcohol.Day2.AWAKEtoNREM.P_NEGFP_Mean;
e2 = data.Alcohol.Day2.AWAKEtoNREM.P_NEGFP_SEM;
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3
x3 = T1;
y3 = data.Alcohol.Day3.AWAKEtoNREM.P_NEGFP_Mean;
e3 = data.Alcohol.Day3.AWAKEtoNREM.P_NEGFP_SEM;
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = T1;
y4 = data.Alcohol.Day4.AWAKEtoNREM.P_NEGFP_Mean;
e4 = data.Alcohol.Day4.AWAKEtoNREM.P_NEGFP_SEM;
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);

% GFP day 5
x5 = T1;
y5 = data.Alcohol.Day5.AWAKEtoNREM.P_NEGFP_Mean;
e5 = data.Alcohol.Day5.AWAKEtoNREM.P_NEGFP_SEM;
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6 
% x6 = T1;
% y6 = data.Alcohol.Day6.AWAKEtoNREM.P_NEGFP_Mean;
% e6 = data.Alcohol.Day6.AWAKEtoNREM.P_NEGFP_SEM;
% fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
% p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);

%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5','Location','northeast','FontSize',5)
legend([p1 p2 p3 p4 p5],'Week1','Week2', 'Week3', 'Week4', 'Week5','Location','northeast','FontSize',5)
xlim([-30,30])
ylim([-0.2 0.6])
title('Awake to NREM')
ylabel('\DeltaF/F SST Ca2+ activity PFC (%)')
xlabel('Time (s)')
set(gca,'box','off')
axis square
%% NREM to Awake

ax2 = subplot(2,2,2);

% GFP day 1
x = T1;
y = data.Alcohol.Day1.NREMtoAWAKE.P_NEGFP_Mean;
e = data.Alcohol.Day1.NREMtoAWAKE.P_NEGFP_SEM;
hold on;
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = T1;
y2 = data.Alcohol.Day2.NREMtoAWAKE.P_NEGFP_Mean;
e2 = data.Alcohol.Day2.NREMtoAWAKE.P_NEGFP_SEM;
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3
x3 = T1;
y3 = data.Alcohol.Day3.NREMtoAWAKE.P_NEGFP_Mean;
e3 = data.Alcohol.Day3.NREMtoAWAKE.P_NEGFP_SEM;
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = T1;
y4 = data.Alcohol.Day4.NREMtoAWAKE.P_NEGFP_Mean;
e4 = data.Alcohol.Day4.NREMtoAWAKE.P_NEGFP_SEM;
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);

% GFP day 5
x5 = T1;
y5 = data.Alcohol.Day5.NREMtoAWAKE.P_NEGFP_Mean;
e5 = data.Alcohol.Day5.NREMtoAWAKE.P_NEGFP_SEM;
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6 
% x6 = T1;
% y6 = data.Alcohol.Day6.NREMtoAWAKE.P_NEGFP_Mean;
% e6 = data.Alcohol.Day6.NREMtoAWAKE.P_NEGFP_SEM;
% fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
% p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);

%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5','Location','northeast','FontSize',5)
legend([p1 p2 p3 p4 p5],'Week1','Week2', 'Week3', 'Week4', 'Week5','Location','northeast','FontSize',5)
xlim([-30,30])
ylim([-0.2 0.6])
title('NREM to Awake')
ylabel('\DeltaF/F SST Ca2+ activity PFC (%)')
xlabel('Time (s)')
set(gca,'box','off')
axis square
%% NREM to REM

ax3 = subplot(2,2,3);

% GFP day 1
x = T1;
y = data.Alcohol.Day1.NREMtoREM.P_NEGFP_Mean;
e = data.Alcohol.Day1.NREMtoREM.P_NEGFP_SEM;
hold on;
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = T1;
y2 = data.Alcohol.Day2.NREMtoREM.P_NEGFP_Mean;
e2 = data.Alcohol.Day2.NREMtoREM.P_NEGFP_SEM;
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3
x3 = T1;
y3 = data.Alcohol.Day3.NREMtoREM.P_NEGFP_Mean;
e3 = data.Alcohol.Day3.NREMtoREM.P_NEGFP_SEM;
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = T1;
y4 = data.Alcohol.Day4.NREMtoREM.P_NEGFP_Mean;
e4 = data.Alcohol.Day4.NREMtoREM.P_NEGFP_SEM;
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);
 
% GFP day 5
x5 = T1;
y5 = data.Alcohol.Day5.NREMtoREM.P_NEGFP_Mean;
e5 = data.Alcohol.Day5.NREMtoREM.P_NEGFP_SEM;
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6 
% x6 = T1;
% y6 = data.Alcohol.Day6.NREMtoREM.P_NEGFP_Mean;
% e6 = data.Alcohol.Day6.NREMtoREM.P_NEGFP_SEM;
% fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
% p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);

% legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5','Location','northeast','FontSize',5)
legend([p1 p2 p3 p4 p5],'Week1','Week2', 'Week3', 'Week4', 'Week5','Location','northeast','FontSize',5)
xlim([-30,30])
ylim([-1.0 2.2])
title('NREM to REM')
ylabel('\DeltaF/F SST Ca2+ activity PFC (%)')
xlabel('Time (s)')
set(gca,'box','off')
axis square
%% REM to AWAKE

ax4 = subplot(2,2,4);

% GFP day 1
x = T1;
y = data.Alcohol.Day1.REMtoAWAKE.P_NEGFP_Mean;
e = data.Alcohol.Day1.REMtoAWAKE.P_NEGFP_SEM;
hold on;
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = T1;
y2 = data.Alcohol.Day2.REMtoAWAKE.P_NEGFP_Mean;
e2 = data.Alcohol.Day2.REMtoAWAKE.P_NEGFP_SEM;
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3
x3 = T1;
y3 = data.Alcohol.Day3.REMtoAWAKE.P_NEGFP_Mean;
e3 = data.Alcohol.Day3.REMtoAWAKE.P_NEGFP_SEM;
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = T1;
y4 = data.Alcohol.Day4.REMtoAWAKE.P_NEGFP_Mean;
e4 = data.Alcohol.Day4.REMtoAWAKE.P_NEGFP_SEM;
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);

% GFP day 5
x5 = T1;
y5 = data.Alcohol.Day5.REMtoAWAKE.P_NEGFP_Mean;
e5 = data.Alcohol.Day5.REMtoAWAKE.P_NEGFP_SEM;
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6 
% x6 = T1;
% y6 = data.Alcohol.Day6.REMtoAWAKE.P_NEGFP_Mean;
% e6 = data.Alcohol.Day6.REMtoAWAKE.P_NEGFP_SEM;
% fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
% p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);

% legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5','Location','northeast','FontSize',5)
legend([p1 p2 p3 p4 p5],'Week1','Week2', 'Week3', 'Week4', 'Week5','Location','northeast','FontSize',5)
xlim([-30,30])
ylim([-1.0 2.2])
title('REM to Awake')
ylabel('\DeltaF/F SST Ca2+ activity PFC (%)')
xlabel('Time (s)')
set(gca,'box','off')
axis square
%% Save figure
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figures Water vs Alcohol' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(fullfile(dirpath, 'ArousalStateTransition_SST_PFC_GFP_Alcohol_ANIMALS'));
end
