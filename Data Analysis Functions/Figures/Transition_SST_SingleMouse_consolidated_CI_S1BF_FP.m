function [AnalysisResults] = Transition_SST_SingleMouse_consolidated_CI_S1BF_FP(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data

    if firstHrs == "false"
         transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
         %transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
         %transitions = {'NREMtoREM','REMtoAWAKE'};
    elseif firstHrs == "true"
        transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    end
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};

        % the size of the matrix
        MatLength = size(AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw,1);
        if aa == 1
            DataLength = 0;
        else
            DataLength = size(data.(transition).AChCBVRaw,1);
        end           

        % data.(transition).EMG(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;

        % data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        % data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        % data.(transition).LH_Cort(:,:,DataLength+1:DataLength+MatLength) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;

        data.(transition).AChCBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw;
        data.(transition).AChGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_AChRaw;
        
        data.(transition).NECBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).NE_CBVRaw;
        data.(transition).NEGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NERaw;
    end
end
% take average for each behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    % data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).AChCBV_N_Exp = size(data.(transition).AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChCBV_Mean = mean(data.(transition).AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AChCBV_SEM = std(data.(transition).AChCBVRaw,1)/sqrt(data.(transition).AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AChCBV_CI95 = tinv([0.025 0.975], data.(transition).AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChCBV_yCI95 = bsxfun(@times, data.(transition).AChCBV_SEM, data.(transition).AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
    
    data.(transition).NECBV_N_Exp = size(data.(transition).NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NECBV_Mean = mean(data.(transition).NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NECBV_SEM = std(data.(transition).NECBVRaw,1)/sqrt(data.(transition).NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NECBV_CI95 = tinv([0.025 0.975], data.(transition).NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NECBV_yCI95 = bsxfun(@times, data.(transition).NECBV_SEM, data.(transition).NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).AChGFP_N_Exp = size(data.(transition).AChGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChGFP_Mean = mean(data.(transition).AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AChGFP_SEM = std(data.(transition).AChGFPRaw,1)/sqrt(data.(transition).AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AChGFP_CI95 = tinv([0.025 0.975], data.(transition).AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChGFP_yCI95 = bsxfun(@times, data.(transition).AChGFP_SEM, data.(transition).AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).NEGFP_N_Exp = size(data.(transition).NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NEGFP_Mean = mean(data.(transition).NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_SEM = std(data.(transition).NEGFPRaw,1)/sqrt(data.(transition).NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_CI95 = tinv([0.025 0.975], data.(transition).NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NEGFP_yCI95 = bsxfun(@times, data.(transition).NEGFP_SEM, data.(transition).NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
%% Fig. 4 Plot CBV and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs CBV');
sgtitle('Blood Volume and Neuromodulators Transition')
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);

%% Awake to NREM
figure (2)
ax1 = subplot(2,2,1);
AwaketoNREM_CBVupperBound = data.AWAKEtoNREM.AChCBV_Mean + data.AWAKEtoNREM.AChCBV_yCI95;
AwaketoNREM_CBVlowerBound = data.AWAKEtoNREM.AChCBV_Mean - data.AWAKEtoNREM.AChCBV_yCI95;
AwaketoNREM_GCaMPupperBound = data.AWAKEtoNREM.AChGFP_Mean + data.AWAKEtoNREM.AChGFP_yCI95;
AwaketoNREM_GCaMPlowerBound = data.AWAKEtoNREM.AChGFP_Mean - data.AWAKEtoNREM.AChGFP_yCI95;

% CBV and SST Ca2+ activity
p1 = plot(T1,data.AWAKEtoNREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

x1 = [T1, fliplr(T1)];
inBetween1 = [AwaketoNREM_CBVlowerBound, fliplr(AwaketoNREM_CBVupperBound)];
fill(x1, inBetween1, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylabel('\DeltaF/F(%)')
xlim([-30,30])

yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on

x2 = [T1, fliplr(T1)];
inBetween2 = [AwaketoNREM_GCaMPlowerBound, fliplr(AwaketoNREM_GCaMPupperBound)];
fill(x2, inBetween2, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'CBV S1BF','SST S1BF')
ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = [0 0.4470 0.7410];
ax1.TickLength = [0.03,0.03];
%ax1.YAxis(1).Limits = [-1.0 1.0];
%ax1.YAxis(2).Limits = [-1.0 1.0];
%Axes control GFP
ax1.YAxis(1).Limits = [-0.8 0.8];
ax1.YAxis(2).Limits = [-0.8 0.8];
axis square

%% NREM to Awake
ax2 = subplot(2,2,2);
NREMtoAWAKE_CBVupperBound = data.NREMtoAWAKE.AChCBV_Mean + data.NREMtoAWAKE.AChCBV_yCI95;
NREMtoAWAKE_CBVlowerBound = data.NREMtoAWAKE.AChCBV_Mean - data.NREMtoAWAKE.AChCBV_yCI95;
NREMtoAWAKE_GCaMPupperBound = data.NREMtoAWAKE.AChGFP_Mean + data.NREMtoAWAKE.AChGFP_yCI95;
NREMtoAWAKE_GCaMPlowerBound = data.NREMtoAWAKE.AChGFP_Mean - data.NREMtoAWAKE.AChGFP_yCI95;

% CBV and SST Ca2+ activity
p1 = plot(T1,data.NREMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

x3 = [T1, fliplr(T1)];
inBetween3 = [NREMtoAWAKE_CBVlowerBound, fliplr(NREMtoAWAKE_CBVupperBound)];
fill(x3, inBetween3, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylabel('\DeltaF/F(%)')
xlim([-30,30])

yyaxis right
p2 = plot(T1,data.NREMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on

x4 = [T1, fliplr(T1)];
inBetween4 = [NREMtoAWAKE_GCaMPlowerBound, fliplr(NREMtoAWAKE_GCaMPupperBound)];
fill(x4, inBetween4, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'CBV S1BF','SST S1BF')
ax2.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax2.YAxis(2).Color = [0 0.4470 0.7410];
ax2.TickLength = [0.03,0.03];
% ax2.YAxis(1).Limits = [-1.0 1.0];
% ax2.YAxis(2).Limits = [-1.0 1.0];
%Axes control GFP
ax2.YAxis(1).Limits = [-0.8 0.8];
ax2.YAxis(2).Limits = [-0.8 0.8];
axis square

%% NREM to REM
ax3 = subplot(2,2,3);
NREMtoREM_CBVupperBound = data.NREMtoREM.AChCBV_Mean + data.NREMtoREM.AChCBV_yCI95;
NREMtoREM_CBVlowerBound = data.NREMtoREM.AChCBV_Mean - data.NREMtoREM.AChCBV_yCI95;
NREMtoREM_GCaMPupperBound = data.NREMtoREM.AChGFP_Mean + data.NREMtoREM.AChGFP_yCI95;
NREMtoREM_GCaMPlowerBound = data.NREMtoREM.AChGFP_Mean - data.NREMtoREM.AChGFP_yCI95;

% CBV and SST Ca2+ activity
p1 = plot(T1,data.NREMtoREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

x5 = [T1, fliplr(T1)];
inBetween5 = [NREMtoREM_CBVlowerBound, fliplr(NREMtoREM_CBVupperBound)];
fill(x5, inBetween5, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylabel('\DeltaF/F(%)')
xlim([-30,30])

yyaxis right
p2 = plot(T1,data.NREMtoREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on

x6 = [T1, fliplr(T1)];
inBetween6 = [NREMtoREM_GCaMPlowerBound, fliplr(NREMtoREM_GCaMPupperBound)];
fill(x6, inBetween6, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'CBV S1BF','SST S1BF')
ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax3.YAxis(2).Color = [0 0.4470 0.7410];
ax3.TickLength = [0.03,0.03];
% ax3.YAxis(1).Limits = [-1.0,2.5];
% ax3.YAxis(2).Limits = [-1.0,2.5];
%Axes control GFP
ax3.YAxis(1).Limits = [-0.5 1.5];
ax3.YAxis(2).Limits = [-0.5 1.0];
axis square

%% REM to Awake
ax4 = subplot(2,2,4);
REMtoAWAKE_CBVupperBound = data.REMtoAWAKE.AChCBV_Mean + data.REMtoAWAKE.AChCBV_yCI95;
REMtoAWAKE_CBVlowerBound = data.REMtoAWAKE.AChCBV_Mean - data.REMtoAWAKE.AChCBV_yCI95;
REMtoAWAKE_GCaMPupperBound = data.REMtoAWAKE.AChGFP_Mean + data.REMtoAWAKE.AChGFP_yCI95;
REMtoAWAKE_GCaMPlowerBound = data.REMtoAWAKE.AChGFP_Mean - data.REMtoAWAKE.AChGFP_yCI95;

% CBV and SST Ca2+ activity
p1 = plot(T1,data.REMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

x7 = [T1, fliplr(T1)];
inBetween7 = [REMtoAWAKE_CBVlowerBound, fliplr(REMtoAWAKE_CBVupperBound)];
fill(x7, inBetween7, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylabel('\DeltaF/F(%)')
xlim([-30,30])

yyaxis right
p2 = plot(T1,data.REMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on

x8 = [T1, fliplr(T1)];
inBetween8 = [REMtoAWAKE_GCaMPlowerBound, fliplr(REMtoAWAKE_GCaMPupperBound)];
fill(x8, inBetween8, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'CBV S1BF','SST S1BF')
ax4.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax4.YAxis(2).Color = [0 0.4470 0.7410];
ax4.TickLength = [0.03,0.03];
% ax4.YAxis(1).Limits = [-1.0,2.5];
% ax4.YAxis(2).Limits = [-1.0,2.5];
%Axes control GFP
ax4.YAxis(1).Limits = [-1.0 2.0];
ax4.YAxis(2).Limits = [-1.0 1.0];
axis square


%% Save figures
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_D,[dirpath ManipulationType '_' 'Transition_CBV_vs_GFP_Sleep_CI_SST_S1BF_noUP']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_' 'Transition_CBV_vs_GFP_Sleep_CI_SST_S1BF_noUP'])
end
%close


end