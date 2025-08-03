function [AnalysisResults] = Fig1_S2_Stim_SST_ResponseIndividual_CI_ToCompare(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
% load('S:\NEACh\AnalysisResults_firstHrs.mat');
% AnalysisResults = AnalysisResults_firstHrs;

FigInital = strfind(rootFolder,'\');
if isempty(FigInital)
    ManipulationType = rootFolder;
else
    ManipulationType = rootFolder(FigInital(end)+1:end);
end

%%%original
%FigInital = strfind(rootFolder,'\');
%ManipulationType = rootFolder(FigInital(end)+1:end);
%% set-up and process data
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.P_NE.(solenoidName).count(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).count;
            % mean
            data.P_NE.(solenoidName).CBV(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).CBV.CBV;
            data.P_ACh.(solenoidName).CBV(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).CBV.CBV;
            data.P_NE.(solenoidName).GFP(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).GFP.GFP;
            data.P_ACh.(solenoidName).GFP(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).GFP.GFP;
            %std
            data.P_NE.(solenoidName).CBV_std(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).CBV.CBVStD;
            data.P_ACh.(solenoidName).CBV_std(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).CBV.CBVStD;
            data.P_NE.(solenoidName).GFP_std(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).GFP.GFPStD;
            data.P_ACh.(solenoidName).GFP_std(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).GFP.GFPStD;

            % neural activity
            data.cortical.(solenoidName).cortMUA(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).MUA.corticalData;
            data.cortical.(solenoidName).cortGam(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).Gam.corticalData;
            data.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
            data.cortical.(solenoidName).cortS(:,:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS;
            data.cortical.(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS(49:end,20:23);
            data.cortical.(solenoidName).T(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.T;
            data.cortical.(solenoidName).F(:,aa) = AnalysisResults.AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.F;
        end
end

for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            % the size of the matrix
            MatLength = size(AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).CBV.CBVRaw,1);
            if aa == 1
                DataLength = 0;
            else
                DataLength = size(data.P_NE.(solenoidName).CBVRaw,1);
            end
            % mean
            data.P_NE.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).CBV.CBVRaw;
            data.P_ACh.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).CBV.CBVRaw;
            data.P_NE.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.AnalysisResults.(animalID).Stim.P_NE.(solenoidName).GFP.GFPRaw;
            data.P_ACh.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.AnalysisResults.(animalID).Stim.P_ACh.(solenoidName).GFP.GFPRaw;
        end
end

% concatenate the data from the contra and ipsi data
% contra
data.Contra.count = data.P_NE.LPadSol.count;
data.Contra.P_AChCBV = data.P_ACh.RPadSol.CBV;
data.Contra.P_AChGFP = data.P_ACh.RPadSol.GFP;
data.Contra.P_NECBV = data.P_NE.LPadSol.CBV;
data.Contra.P_NEGFP = data.P_NE.LPadSol.GFP;

data.Contra.P_AChCBV_std = data.P_ACh.RPadSol.CBV_std;
data.Contra.P_AChGFP_std = data.P_ACh.RPadSol.GFP_std;
data.Contra.P_NECBV_std = data.P_NE.LPadSol.CBV_std;
data.Contra.P_NEGFP_std = data.P_NE.LPadSol.GFP_std;

data.Contra.cortMUA = data.cortical.LPadSol.cortMUA;
data.Contra.cortGam = data.cortical.LPadSol.cortGam;
data.Contra.timeVector = data.cortical.LPadSol.timeVector;
data.Contra.cortS = data.cortical.LPadSol.cortS;
data.Contra.cortS_Gam = data.cortical.LPadSol.cortS_Gam;
data.Contra.T = data.cortical.LPadSol.T;
data.Contra.F = data.cortical.LPadSol.F;

%ipsi
data.Ipsi.count = data.P_NE.RPadSol.count;
data.Ipsi.P_AChCBV_std = data.P_ACh.LPadSol.CBV_std; 
data.Ipsi.P_AChGFP_std = data.P_ACh.LPadSol.GFP_std;
data.Ipsi.P_NECBV_std = data.P_NE.RPadSol.CBV_std;
data.Ipsi.P_NEGFP_std = data.P_NE.RPadSol.GFP_std; 

data.Ipsi.P_AChCBV = data.P_ACh.LPadSol.CBV; 
data.Ipsi.P_AChGFP = data.P_ACh.LPadSol.GFP;
data.Ipsi.P_NECBV = data.P_NE.RPadSol.CBV;
data.Ipsi.P_NEGFP = data.P_NE.RPadSol.GFP; 

data.Ipsi.cortMUA = data.cortical.RPadSol.cortMUA; 
data.Ipsi.cortGam = data.cortical.RPadSol.cortGam; 
data.Ipsi.timeVector = data.cortical.RPadSol.timeVector;
data.Ipsi.cortS = data.cortical.RPadSol.cortS; 
data.Ipsi.cortS_Gam = data.cortical.RPadSol.cortS_Gam;
data.Ipsi.T = data.cortical.RPadSol.T;
data.Ipsi.F = data.cortical.RPadSol.F;

%auditory
data.Auditory.count = data.P_NE.AudSol.count;

data.Auditory.P_NECBV = data.P_NE.AudSol.CBV;
data.Auditory.P_AChCBV = data.P_ACh.AudSol.CBV;
data.Auditory.P_NEGFP = data.P_NE.AudSol.GFP;
data.Auditory.P_AChGFP = data.P_ACh.AudSol.GFP;

data.Auditory.P_NECBV_std = data.P_NE.AudSol.CBV_std;
data.Auditory.P_AChCBV_std = data.P_ACh.AudSol.CBV_std;
data.Auditory.P_NEGFP_std = data.P_NE.AudSol.GFP_std;
data.Auditory.P_AChGFP_std = data.P_ACh.AudSol.GFP_std;

data.Auditory.cortMUA = data.cortical.AudSol.cortMUA;
data.Auditory.cortGam = data.cortical.AudSol.cortGam;
data.Auditory.timeVector = data.cortical.AudSol.timeVector;
data.Auditory.cortS = data.cortical.AudSol.cortS;
data.Auditory.cortS_Gam = data.cortical.AudSol.cortS_Gam;
data.Auditory.T = data.cortical.AudSol.T;
data.Auditory.F = data.cortical.AudSol.F;
% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    
    data.(compDataType).P_AChmean_CBV = mean(data.(compDataType).P_AChCBV,2);
    data.(compDataType).P_AChstd_CBV = std(data.(compDataType).P_AChCBV,0,2);

    data.(compDataType).P_NEmean_CBV = mean(data.(compDataType).P_NECBV,2);
    data.(compDataType).P_NEstd_CBV = std(data.(compDataType).P_NECBV,0,2);

    data.(compDataType).P_AChmean_GFP = mean(data.(compDataType).P_AChGFP,2);
    data.(compDataType).P_AChstd_GFP = std(data.(compDataType).P_AChGFP,0,2);

    data.(compDataType).P_NEmean_GFP = mean(data.(compDataType).P_NEGFP,2);
    data.(compDataType).P_NEstd_GFP = std(data.(compDataType).P_NEGFP,0,2);

    data.(compDataType).P_AChmean_CBV_std = mean(data.(compDataType).P_AChCBV_std,2);
    data.(compDataType).P_NEmean_CBV_std = mean(data.(compDataType).P_NECBV_std,2);

    data.(compDataType).P_AChmean_GFP_std = mean(data.(compDataType).P_AChGFP_std,2);
    data.(compDataType).P_NEmean_GFP_std = mean(data.(compDataType).P_NEGFP_std,2);


    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),3);
    data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end
%% raw data
% concatenate the data from the contra and ipsi data
% contra
data.Contra.P_AChCBVRaw = data.P_ACh.RPadSol.CBVRaw;
data.Contra.P_AChGFPRaw = data.P_ACh.RPadSol.GFPRaw;
data.Contra.P_NECBVRaw = data.P_NE.LPadSol.CBVRaw;
data.Contra.P_NEGFPRaw = data.P_NE.LPadSol.GFPRaw;

%ipsi
data.Ipsi.P_AChCBVRaw = data.P_ACh.LPadSol.CBVRaw; 
data.Ipsi.P_AChGFPRaw = data.P_ACh.LPadSol.GFPRaw;
data.Ipsi.P_NECBVRaw = data.P_NE.RPadSol.CBVRaw;
data.Ipsi.P_NEGFPRaw = data.P_NE.RPadSol.GFPRaw; 

%auditory
data.Auditory.P_NECBVRaw = data.P_NE.AudSol.CBVRaw;
data.Auditory.P_AChCBVRaw = data.P_ACh.AudSol.CBVRaw;
data.Auditory.P_NEGFPRaw = data.P_NE.AudSol.GFPRaw;
data.Auditory.P_AChGFPRaw = data.P_ACh.AudSol.GFPRaw;

% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).P_AChCBV_N_Exp = size(data.(compDataType).P_AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_AChCBV_Mean = mean(data.(compDataType).P_AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChCBV_SEM = std(data.(compDataType).P_AChCBVRaw,1)/sqrt(data.(compDataType).P_AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChCBV_CI95 = tinv([0.025 0.975], data.(compDataType).P_AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_AChCBV_yCI95 = bsxfun(@times, data.(compDataType).P_AChCBV_SEM, data.(compDataType).P_AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_NECBV_N_Exp = size(data.(compDataType).P_NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_NECBV_Mean = mean(data.(compDataType).P_NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NECBV_SEM = std(data.(compDataType).P_NECBVRaw,1)/sqrt(data.(compDataType).P_NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NECBV_CI95 = tinv([0.025 0.975], data.(compDataType).P_NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_NECBV_yCI95 = bsxfun(@times, data.(compDataType).P_NECBV_SEM, data.(compDataType).P_NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_AChGFP_N_Exp = size(data.(compDataType).P_AChGFPRaw,1);
    data.(compDataType).P_AChGFP_Mean = mean(data.(compDataType).P_AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChGFP_SEM = std(data.(compDataType).P_AChGFPRaw,1)/sqrt(data.(compDataType).P_AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChGFP_CI95 = tinv([0.025 0.975], data.(compDataType).P_AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_AChGFP_yCI95 = bsxfun(@times, data.(compDataType).P_AChGFP_SEM, data.(compDataType).P_AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_NEGFP_N_Exp = size(data.(compDataType).P_NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_NEGFP_Mean = mean(data.(compDataType).P_NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NEGFP_SEM = std(data.(compDataType).P_NEGFPRaw,1)/sqrt(data.(compDataType).P_NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NEGFP_CI95 = tinv([0.025 0.975], data.(compDataType).P_NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_NEGFP_yCI95 = bsxfun(@times, data.(compDataType).P_NEGFP_SEM, data.(compDataType).P_NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
end
%% plot the comparison of GCaMP and TRITC with SST Ca2+ activity and TRITC
summaryFigureN = figure ;
sgtitle('Stimulus evoked responses in fiber photometry signals')

%This what happens to the LH when you puff the right whiskers
%LH (S1BF) contra stim 

%CBV
ax1 = subplot(2,3,1);
plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean + data.Contra.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean - data.Contra.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Contra stim S1BF fiber')
ylabel('\DeltaF/F S1BF CBV (%)')
ax1.YLim = [-1 2];

%GFP
yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean + data.Contra.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean - data.Contra.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F S1BF SST Ca2+ acitvity (%)')

ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YLim = [-1 3];
xlim([-5 15])
axis square

%This is what happens to the LH when you puff the left whiskers 
%LH (SIBF) ispi stim

%CBV
ax2 = subplot(2,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChCBV_Mean + data.Ipsi.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChCBV_Mean - data.Ipsi.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Ipsi stim S1BF fiber')
ylabel('\DeltaF/F S1BF CBV (%)')
ax2.YLim = [-1.5 1];

%GFP
yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChGFP_Mean + data.Ipsi.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_AChGFP_Mean - data.Ipsi.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F S1BF SST Ca+ activity (%)')
ax2.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax2.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-1.5 1];
xlim([-5 15])
axis square

%LH (SIBF) auditory stim

%CBV
ax3 = subplot(2,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChCBV_Mean + data.Auditory.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChCBV_Mean - data.Auditory.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Auditory stim S1BF fiber ')
ylabel('\DeltaF/F S1BF CBV (%)')
ax3.YLim = [-3 3];

%GFP
yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChGFP_Mean + data.Auditory.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.P_AChGFP_Mean - data.Auditory.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F S1BF SSt Ca2+ activity (%)')
ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax3.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YLim = [-3 3];
xlim([-5 15])
axis square

%This is what happens to the RH when you stimulate the left whiskers
% RH(PFC) contra stim

%CBV
ax4 = subplot(2,3,4);
plot(data.Contra.mean_timeVector,data.Contra.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.P_NECBV_Mean + data.Contra.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.P_NECBV_Mean - data.Contra.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Contra stim PFC fiber')
ylabel('\DeltaF/F PFC CBV (%)')
ax4.YLim = [-0.5 5.0];

%GFP
yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.P_NEGFP_Mean + data.Contra.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.P_NEGFP_Mean - data.Contra.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F PFC SST Ca2+ activity (%)')
ax4.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax4.YAxis(2).Color = [0.4660 0.6740 0.1880];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
ax4.YLim = [-0.5 3.5];
xlim([-5 15])
axis square

% This is what happens to the RH when you puff the right whiskers
% RH(PFC) ispi stim

%CBV
ax5 = subplot(2,3,5);
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NECBV_Mean + data.Ipsi.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NECBV_Mean - data.Ipsi.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Ipsi stim PFC fiber')
ylabel('\DeltaF/F PFC CBV (%)')
ax5.YLim = [-0.5 5.0];

%GFP
yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NEGFP_Mean + data.Ipsi.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.P_NEGFP_Mean - data.Ipsi.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F PFC SST Ca2+ activity (%)')
ax5.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax5.YAxis(2).Color = [0.4660 0.6740 0.1880];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
ax5.YLim = [-0.5 3.5];
xlim([-5 15])
axis square

% RH (PFC) auditory stim

%CBV
ax6 = subplot(2,3,6);
plot(data.Auditory.mean_timeVector,data.Auditory.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.P_NECBV_Mean + data.Auditory.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.P_NECBV_Mean - data.Auditory.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Auditory stim PFC fiber')
ylabel('\DeltaF/F PFC CBV (Z)')
ax6.YLim = [-3 3];

%GFP
yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.P_NEGFP_Mean + data.Auditory.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.P_NEGFP_Mean - data.Auditory.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F PFC SST Ca2+ activity (%)')
ax6.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax6.YAxis(2).Color = [0.4660 0.6740 0.1880];
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
ax6.YLim = [-3 3];
xlim([-5 15])
axis square
%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end    
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath ManipulationType '_' 'Stim-FiberSignals_ResponseIndividual_CI_SST']); %animalID   
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_'  'Stim-FiberSignals_ResponseIndividual_CI_SST_UP']) % animalID _NoSleep NoSleep_
    close 
end
end