
%% set-up and process data
IOSanimalIDs = {'FS54','FS60'};

% load('All_Transitions.mat');
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%%
transitionsFileStruct = dir('*_Results_Transition.mat');
transitionsFiles = {transitionsFileStruct.name}';

All_Transitions = struct();

for a = 1:length(transitionsFiles)
    transitionsFileID = transitionsFiles{a};
    
    % Load the .mat file
    transitionsloadedData = load(transitionsFileID);

    % Extract the AnimalID from the filename
    % Assuming the filename format is like 'fs54_RestingBaselines.mat'
    % Adjust the regexp according to your filename format
    animalID = regexp(transitionsFileID, 'FS\d+', 'match');
    if isempty(animalID)
        continue; % Skip if no AnimalID is found
    end
    animalID = animalID{1};
    
    % Assuming the structure inside the .mat file has a consistent name,
    % otherwise you might need to adjust this part.
    structName = fieldnames(transitionsloadedData);
    data = transitionsloadedData.(structName{1});
    
    % Add the loaded data to the All_RestingBaseline structure under a new field
    All_Transitions.(animalID) = data;
end

save('All_both_Transitions.mat','All_Transitions');

%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
% awake to nrem 
All_AWAKEtoNREM_GCaMP_animal1 = cat(1, All_Transitions.FS54.transition1.FS54.Transitions.AWAKEtoNREM.both_GCaMP,All_Transitions.FS54.transition2.FS54.Transitions.AWAKEtoNREM.both_GCaMP);
All_AWAKEtoNREM_GCaMP_animal2 = cat(1, All_Transitions.FS60.transition1.FS60.Transitions.AWAKEtoNREM.both_GCaMP,All_Transitions.FS60.transition2.FS60.Transitions.AWAKEtoNREM.both_GCaMP);

All_GCaMP_AWAKEtoNREM = cat(1,All_AWAKEtoNREM_GCaMP_animal1, All_AWAKEtoNREM_GCaMP_animal2);
All_mean_GCaMP_AWAKEtoNREM = (mean(All_GCaMP_AWAKEtoNREM)-1)*100;
All_std_GCaMP_AWAKEtoNREM = (std(All_GCaMP_AWAKEtoNREM))*100;

All_AWAKEtoNREM_HbT_animal1 = cat(1, All_Transitions.FS54.transition1.FS54.Transitions.AWAKEtoNREM.both_HbT,All_Transitions.FS54.transition2.FS54.Transitions.AWAKEtoNREM.both_HbT);
All_AWAKEtoNREM_HbT_animal2 = cat(1, All_Transitions.FS60.transition1.FS60.Transitions.AWAKEtoNREM.both_HbT,All_Transitions.FS60.transition2.FS60.Transitions.AWAKEtoNREM.both_HbT);

All_HbT_AWAKEtoNREM = cat(1,All_AWAKEtoNREM_HbT_animal1, All_AWAKEtoNREM_HbT_animal2);
All_mean_HbT_AWAKEtoNREM = mean(All_HbT_AWAKEtoNREM);
All_std_HbT_AWAKEtoNREM = std(All_HbT_AWAKEtoNREM);


%  nrem to awake NREMtoAWAKE
All_NREMtoAWAKE_GCaMP_animal1 = All_Transitions.FS54.transition1.FS54.Transitions.NREMtoAWAKE.both_GCaMP;
All_NREMtoAWAKE_GCaMP_animal2 = cat(1,All_Transitions.FS60.transition1.FS60.Transitions.NREMtoAWAKE.both_GCaMP,All_Transitions.FS60.transition2.FS60.Transitions.NREMtoAWAKE.both_GCaMP);
All_GCaMP_NREMtoAWAKE = cat(1,All_NREMtoAWAKE_GCaMP_animal1, All_NREMtoAWAKE_GCaMP_animal2);
All_mean_GCaMP_NREMtoAWAKE = (mean(All_GCaMP_NREMtoAWAKE)- 1)*100;
All_std_GCaMP_NREMtoAWAKE = (std(All_GCaMP_NREMtoAWAKE))*100;

All_NREMtoAWAKE_HbT_animal1 = All_Transitions.FS54.transition1.FS54.Transitions.NREMtoAWAKE.both_HbT;
All_NREMtoAWAKE_HbT_animal2 = cat(1, All_Transitions.FS60.transition1.FS60.Transitions.NREMtoAWAKE.both_HbT,All_Transitions.FS60.transition2.FS60.Transitions.NREMtoAWAKE.both_HbT);

All_HbT_NREMtoAWAKE = cat(1,All_NREMtoAWAKE_HbT_animal1, All_NREMtoAWAKE_HbT_animal2);
All_mean_HbT_NREMtoAWAKE = mean(All_HbT_NREMtoAWAKE);
All_std_HbT_NREMtoAWAKE = std(All_HbT_NREMtoAWAKE);

% nrem to rem
All_NREMtoREM_GCaMP_animal1 = All_Transitions.FS54.transition1.FS54.Transitions.NREMtoREM.both_GCaMP;
All_NREMtoREM_GCaMP_animal2 = cat(1,All_Transitions.FS60.transition1.FS60.Transitions.NREMtoREM.both_GCaMP,All_Transitions.FS60.transition2.FS60.Transitions.NREMtoREM.both_GCaMP);

All_GCaMP_NREMtoREM = cat(1,All_NREMtoREM_GCaMP_animal1, All_NREMtoREM_GCaMP_animal2);
All_mean_GCaMP_NREMtoREM = (mean(All_GCaMP_NREMtoREM)-1)*100;
All_std_GCaMP_NREMtoREM = (std(All_GCaMP_NREMtoREM))*100;

All_NREMtoREM_HbT_animal1 =  All_Transitions.FS54.transition1.FS54.Transitions.NREMtoREM.both_HbT;
All_NREMtoREM_HbT_animal2 = cat(1,All_Transitions.FS60.transition1.FS60.Transitions.NREMtoREM.both_HbT,All_Transitions.FS60.transition2.FS60.Transitions.NREMtoREM.both_HbT);

All_HbT_NREMtoREM = cat(1,All_NREMtoREM_HbT_animal1, All_NREMtoREM_HbT_animal2);
All_mean_HbT_NREMtoREM = mean(All_HbT_NREMtoREM);
All_std_HbT_NREMtoREM = std(All_HbT_NREMtoREM);

%  rem to awake REMtoAWAKE
All_REMtoAWAKE_GCaMP_animal1 = All_Transitions.FS54.transition1.FS54.Transitions.REMtoAWAKE.both_GCaMP;
All_REMtoAWAKE_GCaMP_animal2 = cat(1,All_Transitions.FS60.transition1.FS60.Transitions.REMtoAWAKE.both_GCaMP,All_Transitions.FS60.transition2.FS60.Transitions.REMtoAWAKE.both_GCaMP);

All_GCaMP_REMtoAWAKE = cat(1,All_REMtoAWAKE_GCaMP_animal1, All_REMtoAWAKE_GCaMP_animal2);
All_mean_GCaMP_REMtoAWAKE = (mean(All_GCaMP_REMtoAWAKE)-1)*100;
All_std_GCaMP_REMtoAWAKE = (std(All_GCaMP_REMtoAWAKE))*100;

All_REMtoAWAKE_HbT_animal1 =  All_Transitions.FS54.transition1.FS54.Transitions.REMtoAWAKE.both_HbT;
All_REMtoAWAKE_HbT_animal2 = cat(1,All_Transitions.FS60.transition1.FS60.Transitions.REMtoAWAKE.both_HbT,All_Transitions.FS60.transition2.FS60.Transitions.REMtoAWAKE.both_HbT);

All_HbT_REMtoAWAKE = cat(1,All_REMtoAWAKE_HbT_animal1, All_REMtoAWAKE_HbT_animal2);
All_mean_HbT_REMtoAWAKE = mean(All_HbT_REMtoAWAKE);
All_std_HbT_REMtoAWAKE = std(All_HbT_REMtoAWAKE);

T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;

% 
% awake2nrem_meanGCaMP = (data.AWAKEtoNREM.meanGCaMP - 1)*100;
% awake2nrem_stdError = data.AWAKEtoNREM.stdGCaMP / sqrt(2);
% awake2nrem_stdGCaMP = (awake2nrem_stdError)*100;
% 
% nrem2awake_meanGCaMP = (data.NREMtoAWAKE.meanGCaMP - 1)*100;
% nrem2awake_stdError = data.NREMtoAWAKE.stdGCaMP / sqrt(2);
% nrem2awake_stdGCaMP = (nrem2awake_stdError)*100;
% 
% nrem2rem_meanGCaMP = (data.NREMtoREM.meanGCaMP -1)*100;
% nrem2rem_stdError = data.NREMtoREM.stdGCaMP / sqrt(2);
% nrem2rem_stdGCaMP = (nrem2rem_stdError)*100;
% 
% rem2awake_meanGCaMP = (data.REMtoAWAKE.meanGCaMP - 1)*100;
% rem2awake__stdError = data.REMtoAWAKE.stdGCaMP / sqrt(2);
% rem2awake_stdGCaMP = (rem2awake__stdError)*100;




%% Awake to NREM
figure (2)
ax1 = subplot(1,2,1:2);
AwaketoNREM_HbTupperBound = All_mean_HbT_AWAKEtoNREM + All_std_HbT_AWAKEtoNREM;
AwaketoNREM_HbTlowerBound = All_mean_HbT_AWAKEtoNREM - All_std_HbT_AWAKEtoNREM;
AwaketoNREM_GCaMPupperBound = All_mean_GCaMP_AWAKEtoNREM + All_std_GCaMP_AWAKEtoNREM;
AwaketoNREM_GCaMPlowerBound = All_mean_GCaMP_AWAKEtoNREM - All_std_GCaMP_AWAKEtoNREM;

% HbT and EMG
p1 = plot(T2,All_mean_HbT_AWAKEtoNREM,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
% plot(T2,All_mean_HbT_AWAKEtoNREM + data.AWAKEtoNREM.HbT_yCI95,'-','color',colors('dark candy apple red'),'LineWidth',0.1);

x2 = [T2 , fliplr(T2)];
inBetween = [AwaketoNREM_HbTlowerBound, fliplr(AwaketoNREM_HbTupperBound)];
fill(x2, inBetween, colors('dark candy apple red'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_HbT_AWAKEtoNREM + All_std_HbT_AWAKEtoNREM,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
% plot(T2,All_mean_HbT_AWAKEtoNREM - All_std_HbT_AWAKEtoNREM,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
% ylim([-20,95])
yyaxis right
p2 = plot(T2,All_mean_GCaMP_AWAKEtoNREM,'-','color',colors('dark spring green'),'LineWidth',2);
hold on
% plot(T2,data.AWAKEtoNREM.meanGCaMP + data.AWAKEtoNREM.GCaMP_yCI95,'-','color',colors('dark spring green'),'LineWidth',0.1);

x3 = [T2 , fliplr(T2)];
inBetween = [AwaketoNREM_GCaMPlowerBound, fliplr(AwaketoNREM_GCaMPupperBound)];
fill(x3, inBetween, colors('dark spring green'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_GCaMP_AWAKEtoNREM + All_std_GCaMP_AWAKEtoNREM,'-','color',colors('dark spring green'),'LineWidth',0.1);
% plot(T2,All_mean_GCaMP_AWAKEtoNREM - All_std_GCaMP_AWAKEtoNREM,'-','color',colors('dark spring green'),'LineWidth',0.1);

title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'HbT','GCaMP')
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('dark spring green');
% ylim([-5,10])
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Limits = [-20,95];
ax1.YAxis(2).Limits =[-5,10];
axis square
% % cort neural
% ax2 = subplot(3,1,3);
% Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
% axis xy
% % c1 = colorbar;
% ylabel('\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-1,2])
% xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% % hippocampal neural
% ax3 = subplot(6,2,5);
% semilog_imagesc_eLife2020(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,200])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
%% [4c] NREM to Awake

NREMtoAwake_HbTupperBound = All_mean_HbT_NREMtoAWAKE + All_std_HbT_NREMtoAWAKE;
NREMtoAwake_HbTlowerBound = All_mean_HbT_NREMtoAWAKE - All_std_HbT_NREMtoAWAKE;
NREMtoAwake_GCaMPupperBound = All_mean_GCaMP_NREMtoAWAKE + All_std_GCaMP_NREMtoAWAKE;
NREMtoAwake_GCaMPlowerBound = All_mean_GCaMP_NREMtoAWAKE - All_std_GCaMP_NREMtoAWAKE;
figure (3)
ax4 = subplot(1,2,1:2);
% HbT and EMG
plot(T2,All_mean_HbT_NREMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
% plot(T2,data.NREMtoAWAKE.meanHbT + data.NREMtoAWAKE.HbT_yCI95,'-','color',colors('dark candy apple red'),'LineWidth',0.1);

x4 = [T2 , fliplr(T2)];
inBetween = [NREMtoAwake_HbTlowerBound, fliplr(NREMtoAwake_HbTupperBound)];
fill(x4, inBetween, colors('dark candy apple red'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_HbT_NREMtoAWAKE + All_std_HbT_NREMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
% plot(T2,All_mean_HbT_NREMtoAWAKE - All_std_HbT_NREMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T2,All_mean_GCaMP_NREMtoAWAKE ,'-','color',colors('dark spring green'),'LineWidth',2);
hold on
% plot(T2,data.NREMtoAWAKE.meanGCaMP + data.NREMtoAWAKE.GCaMP_yCI95,'-','color',colors('dark spring green'),'LineWidth',0.1);

x5 = [T2 , fliplr(T2)];
inBetween = [NREMtoAwake_GCaMPlowerBound, fliplr(NREMtoAwake_GCaMPupperBound)];
fill(x5, inBetween, colors('dark spring green'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_GCaMP_NREMtoAWAKE + All_std_GCaMP_NREMtoAWAKE,'-','color',colors('rich black'),'LineWidth',0.1);
% plot(T2,All_mean_GCaMP_NREMtoAWAKE - All_std_GCaMP_NREMtoAWAKE,'-','color',colors('rich black'),'LineWidth',0.1);
title(' NREM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = colors('dark spring green');
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
ax4.YAxis(1).Limits = [-20,95];
ax4.YAxis(2).Limits =[-5,10];
axis square

% % cort neural
% ax5 = subplot(3,1,3);
% Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
% axis xy
% % c3 = colorbar;
% ylabel('\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-1,2])
% xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% % % hippocampal neural
% ax6 = subplot(6,2,6);
% semilog_imagesc_eLife2020(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,200])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
%% [4d] NREM to REM

NREMtoREM_HbTupperBound = All_mean_HbT_NREMtoREM + All_std_HbT_NREMtoREM;
NREMtoREM_HbTlowerBound = All_mean_HbT_NREMtoREM - All_std_HbT_NREMtoREM;
NREMtoREM_GCaMPupperBound = All_mean_GCaMP_NREMtoREM + All_std_GCaMP_NREMtoREM;
NREMtoREM_GCaMPlowerBound = All_mean_GCaMP_NREMtoREM - All_std_GCaMP_NREMtoREM;

figure (4)
ax7 = subplot(1,2,1:2);
% HbT and EMG
plot(T2,All_mean_HbT_NREMtoREM,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
% plot(T2,data.NREMtoREM.meanHbT + data.NREMtoREM.HbT_yCI95,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
% 
x6 = [T2 , fliplr(T2)];
inBetween = [NREMtoREM_HbTlowerBound, fliplr(NREMtoREM_HbTupperBound)];
fill(x6, inBetween, colors('dark candy apple red'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % plot(T2,data.NREMtoREM.meanHbT + data.NREMtoREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
% % plot(T2,data.NREMtoREM.meanHbT - data.NREMtoREM.stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
% % ylim([35,90])
yyaxis right
plot(T2,All_mean_GCaMP_NREMtoREM ,'-','color',colors('dark spring green'),'LineWidth',2);
hold on
% plot(T2,data.NREMtoREM.meanGCaMP + data.NREMtoREM.GCaMP_yCI95,'-','color',colors('dark spring green'),'LineWidth',0.1);
% 
x7 = [T2 , fliplr(T2)];
inBetween = [NREMtoREM_GCaMPlowerBound, fliplr(NREMtoREM_GCaMPupperBound)];
fill(x7, inBetween, colors('dark spring green'), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% % plot(T2,nrem2rem_meanGCaMP + nrem2rem_stdGCaMP,'-','color',colors('dark spring green'),'LineWidth',0.1);
% % plot(T2,nrem2rem_meanGCaMP - nrem2rem_stdGCaMP,'-','color',colors('dark spring green'),'LineWidth',0.1);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('dark candy apple red');
ax7.YAxis(2).Color = colors('dark spring green');
% ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
ax7.YAxis(1).Limits = [-20,95];
ax7.YAxis(2).Limits =[-5,10];
axis square

% cort neural
% ax8 = subplot(3,1,3);
% Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
% axis xy
% % c5 = colorbar;
% ylabel('\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,300])
% xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];
% % hippocampal neural
% ax9 = subplot(6,2,11);
% semilog_imagesc_eLife2020(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,300])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax9.TickLength = [0.03,0.03];
%% REM to Awake
REMtoAWAKE_HbTupperBound = All_mean_HbT_REMtoAWAKE + All_std_HbT_REMtoAWAKE;
REMtoAWAKE_HbTlowerBound = All_mean_HbT_REMtoAWAKE - All_std_HbT_REMtoAWAKE;
REMtoAWAKE_GCaMPupperBound = All_mean_GCaMP_REMtoAWAKE + All_std_GCaMP_REMtoAWAKE;
REMtoAWAKE_GCaMPlowerBound = All_mean_GCaMP_REMtoAWAKE - All_std_GCaMP_REMtoAWAKE;

figure (5)
ax10 = subplot(1,2,1:2);
plot(T2,All_mean_HbT_REMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
% plot(T2,data.REMtoAWAKE.meanHbT + data.REMtoAWAKE.HbT_yCI95,'-','color',colors('dark candy apple red'),'LineWidth',0.1);

x8 = [T2 , fliplr(T2)];
inBetween = [REMtoAWAKE_HbTlowerBound, fliplr(REMtoAWAKE_HbTupperBound)];
fill(x8, inBetween, colors('dark candy apple red'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_HbT_REMtoAWAKE + All_std_HbT_REMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
% plot(T2,All_mean_HbT_REMtoAWAKE - All_std_HbT_REMtoAWAKE,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta[HbT] (\muM)')
xlim([-30,30])
% ylim([0,90])
yyaxis right
plot(T2,All_mean_GCaMP_REMtoAWAKE ,'-','color',colors('dark spring green'),'LineWidth',2);
% plot(T2,data.REMtoAWAKE.meanGCaMP + data.REMtoAWAKE.GCaMP_yCI95,'-','color',colors('dark spring green'),'LineWidth',0.1);

hold on
x9 = [T2 , fliplr(T2)];
inBetween = [REMtoAWAKE_GCaMPlowerBound, fliplr(REMtoAWAKE_GCaMPupperBound)];
fill(x9, inBetween, colors('dark spring green'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(T2,All_mean_GCaMP_REMtoAWAKE + All_std_GCaMP_REMtoAWAKE,'-','color',colors('dark spring green'),'LineWidth',0.1);
% plot(T2,All_mean_GCaMP_REMtoAWAKE - All_std_GCaMP_REMtoAWAKE,'-','color',colors('dark spring green'),'LineWidth',0.1);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('dark candy apple red');
ax10.YAxis(2).Color = colors('dark spring green');
ax10.TickLength = [0.03,0.03];
ax10.YAxis(1).Limits = [-20,95];
ax10.YAxis(2).Limits =[-5,10];
% % cort neural
% ax11 = subplot(3,1,3);
% semilog_imagesc_eLife2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
% axis xy
% % c7 = colorbar;
% ylabel('\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-1,2])
% xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax11.TickLength = [0.03,0.03];
% % hippocampal neural
% ax12 = subplot(6,2,12);
% semilog_imagesc_eLife2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
% c8 = colorbar;
% ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,300])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];
% %% axes positionns
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% % ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% % ax6Pos = get(ax6,'position');
% ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
% % ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% % ax12Pos = get(ax12,'position');
% ax2Pos(3:4) = ax1Pos(3:4);
% % ax3Pos(3:4) = ax1Pos(3:4);
% ax5Pos(3:4) = ax4Pos(3:4);
% % ax6Pos(3:4) = ax4Pos(3:4);
% ax8Pos(3:4) = ax7Pos(3:4);
% % ax9Pos(3:4) = ax7Pos(3:4);
% ax11Pos(3:4) = ax10Pos(3:4);
% % ax12Pos(3:4) = ax10Pos(3:4);
% set(ax2,'position',ax2Pos);
% % set(ax3,'position',ax3Pos);
% set(ax5,'position',ax5Pos);
% % set(ax6,'position',ax6Pos);
% set(ax8,'position',ax8Pos);
% % set(ax9,'position',ax9Pos);
% set(ax11,'position',ax11Pos);
% % set(ax12,'position',ax12Pos);

