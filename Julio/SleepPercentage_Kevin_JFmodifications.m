% rootfolder = cd;
%% extract data from each animal's sleep scoring results
    for aa = 1:length(animalIDs)
        animalID = animalIDs{1,aa};
        dataLoc = [rootFolder '/' animalID '/CombinedImaging/'];
        cd(dataLoc)
        scoringResults = [animalID '_Manual_ScoringResults.mat'];
        load(scoringResults,'-mat')
        numberOfScores(aa,1) = length(ScoringResults.alllabels); %#ok<*AGROW,*NASGU>
        indAwakePerc(aa,1) = round((sum(strcmp(ScoringResults.alllabels,'Not Sleep'))/length(ScoringResults.alllabels))*100,1);
        indNremPerc(aa,1) = round((sum(strcmp(ScoringResults.alllabels,'NREM Sleep'))/length(ScoringResults.alllabels))*100,1);
        indRemPerc(aa,1) = round((sum(strcmp(ScoringResults.alllabels,'REM Sleep'))/length(ScoringResults.alllabels))*100,1);
        allCatLabels = vertcat(allCatLabels,ScoringResults.alllabels);
    end
    labels = {'Awake','NREM','REM'};
    % mean percentage of each state between animals
    meanAwakePerc = mean(indAwakePerc,1);
    stdAwakePerc = std(indAwakePerc,0,1);
    meanNremPerc = mean(indNremPerc,1);
    stdNremPerc = std(indNremPerc,0,1);
    meanRemPerc = mean(indRemPerc,1);
    stdRemPerc = std(indRemPerc,0,1);
    meanPercs = horzcat(meanAwakePerc,meanNremPerc,meanRemPerc);
    % percentage of each state for all labels together
    allAwakePerc = round((sum(strcmp(allCatLabels,'Not Sleep'))/length(allCatLabels))*100,1);
    allNremPerc = round((sum(strcmp(allCatLabels,'NREM Sleep'))/length(allCatLabels))*100,1);
    allRemPerc = round((sum(strcmp(allCatLabels,'REM Sleep'))/length(allCatLabels))*100,1);
    meanAllPercs = horzcat(allAwakePerc,allNremPerc,allRemPerc);
    % total time per animal behavioral states
    labelTime = 5;   % seconds
    IOS_indTotalTimeHours = ((numberOfScores*labelTime)/60)/60;
    IOS_allTimeHours = sum(IOS_indTotalTimeHours);
    IOS_meanTimeHours = mean(IOS_indTotalTimeHours,1);
    IOS_stdTimeHours = std(IOS_indTotalTimeHours,0,1);
    allTimeDays = sum(IOS_indTotalTimeHours)/24;
    totalTimeAwake = IOS_indTotalTimeHours.*(indAwakePerc/100);
    meanAwakeHours = mean(totalTimeAwake,1);
    stdAwakeHours = std(totalTimeAwake,0,1);
    totalTimeNREM = IOS_indTotalTimeHours.*(indNremPerc/100);
    meanNREMHours = mean(totalTimeNREM);
    stdNREMHours = std(totalTimeNREM,0,1);
    totalTimeREM = IOS_indTotalTimeHours.*(indRemPerc/100);
    meanREMHours = mean(totalTimeREM);
    stdREMHours = std(totalTimeREM,0,1);
    % update analysis structure
    AnalysisResults.IOSarousalProbability.indAwakePerc = indAwakePerc;
    AnalysisResults.IOSarousalProbability.indNremPerc = indNremPerc;
    AnalysisResults.IOSarousalProbability.indRemPerc = indRemPerc;
    AnalysisResults.IOSarousalProbability.IOS_indTotalTimeHours = IOS_indTotalTimeHours;
    AnalysisResults.IOSarousalProbability.IOS_allTimeHours = IOS_allTimeHours;
    AnalysisResults.IOSarousalProbability.IOS_meanTimeHours = IOS_meanTimeHours;
    AnalysisResults.IOSarousalProbability.IOS_stdTimeHours = IOS_stdTimeHours;
    AnalysisResults.IOSarousalProbability.allTimeDays = allTimeDays;
    AnalysisResults.IOSarousalProbability.meanPercs = meanPercs;
    AnalysisResults.IOSarousalProbability.totalTimeAwake = totalTimeAwake;
    AnalysisResults.IOSarousalProbability.meanAwakeHours = meanAwakeHours;
    AnalysisResults.IOSarousalProbability.stdAwakeHours = stdAwakeHours;
    AnalysisResults.IOSarousalProbability.totalTimeNREM = totalTimeNREM;
    AnalysisResults.IOSarousalProbability.meanNREMHours = meanNREMHours;
    AnalysisResults.IOSarousalProbability.stdNREMHours = stdNREMHours;
    AnalysisResults.IOSarousalProbability.totalTimeREM = totalTimeREM;
    AnalysisResults.IOSarousalProbability.meanREMHours = meanREMHours;
    AnalysisResults.IOSarousalProbability.stdREMHours = stdREMHours;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. 2 (part 2)
summaryFigureB = figure('Name','Fig2 (b-i)');
sgtitle('Figure 2 - Turner et al. 2020')
%% [2b] percentage of arousal-state scores
ax1 = subplot(2,4,1);
p1 = pie(meanPercs);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'rfc-Awake: ';'rfc-NREM: ';'rfc-REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'[2b] Sleep scoring label probability','Mean animal sleep scoring labels'})
%% [2c] ternary plot
ax2 = subplot(2,4,2);
terplot();
[hd] = ternaryc(indAwakePerc/100,indNremPerc/100,indRemPerc/100);
hlabels = terlabel('rfc-Awake','rfc-NREM','rfc-REM');
title({'[2c] Ternary plot of ind animals',' ',' '})
% %% [2d] arousal-state probability over time
% ax3 = subplot(2,4,3);
% xinds1 = (1:numBins)/(numBins/3);
% % awake
% scatter(xinds1,binnedAwakeProbability,25,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcAwake);
% [awakeHypExpCurve,~] = fit(xinds1',binnedAwakeProbability','exp1','StartPoint',[0,0]);
% awakeHypExpFit = awakeHypExpCurve(xinds1);
% hold on
% plot(xinds1,awakeHypExpFit,'color',colorRfcAwake,'LineWidth',2);
% % nrem
% scatter(xinds1,binnedNREMProbability,25,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcNREM);
% [nremHypExpCurve,~] = fit(xinds1',binnedNREMProbability','exp1','StartPoint',[0,0]);
% nremHypExpFit = nremHypExpCurve(xinds1);
% plot(xinds1,nremHypExpFit,'color',colorRfcNREM,'LineWidth',2);
% % rem
% scatter(xinds1,binnedREMProbability,25,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcREM);
% [remHypExpCurve,~] = fit(xinds1',binnedREMProbability','exp1','StartPoint',[0,0]);
% remHypExpFit = remHypExpCurve(xinds1);
% plot(xinds1,remHypExpFit,'color',colorRfcREM,'LineWidth',2);
% title({'[2d] Imaging timeline','Awake probability'})
% xlabel('Time (Hr)')
% ylabel('Probability')
% ylim([0,1])
% set(gca,'box','off')
% axis square
% ax3.TickLength = [0.03,0.03];
% %% [2e] true rest event probability 
% ax4 = subplot(2,4,4);
% xinds2 = 0:length(bins) - 1;
% xinds3 = 0:0.01:length(bins) - 1;
% scatter(xinds2,restEventProbability,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRfcAwake);
% [restExpCurve,~] = fit(xinds2',restEventProbability,'exp2');
% restExpFit = restExpCurve(xinds3);
% hold on
% plot(xinds3,restExpFit,'k','LineWidth',2);
% xticks([1,3,5,7,9,11])
% xticklabels({'10','20','30','40','50','60'})
% title({'[2e] Awake probability','of ''Rest'' events'})
% xlabel('Duration (s)')
% ylabel('Probability')
% xlim([0,12])
% ylim([0,1])
% set(gca,'box','off')
% axis square
% ax4.TickLength = [0.03,0.03];
%% [2f] EMG during different arousal states
% ax5 = subplot(2,4,5);
% edges = -2.5:0.5:2.5;
% [curve1] = SmoothHistogramBins(data.BehavioralDistributions.Awake.catEMG,edges);
% [curve2] = SmoothHistogramBins(data.BehavioralDistributions.NREM.catEMG,edges);
% [curve3] = SmoothHistogramBins(data.BehavioralDistributions.REM.catEMG,edges);
% before = findall(gca);
% fnplt(curve1);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcAwake)
% hold on
% before = findall(gca);
% fnplt(curve2);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcNREM)
% before = findall(gca);
% fnplt(curve3);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcREM)
% title({'[2f] EMG power','arousal-state distribution'})
% xlabel('EMG log10(pwr)')
% ylabel('Probability')
% xlim([-2.5,2.5])
% ylim([0,0.5])
% axis square
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
%% [2g] whisker variance distribution during different arousal-states
% ax6 = subplot(2,4,6);
% edges = -3:0.75:3;
% [curve1] = SmoothLogHistogramBins(data.BehavioralDistributions.Awake.catWhisk,edges);
% [curve2] = SmoothLogHistogramBins(data.BehavioralDistributions.NREM.catWhisk,edges);
% [curve3] = SmoothLogHistogramBins(data.BehavioralDistributions.REM.catWhisk,edges);
% before = findall(gca);
% fnplt(curve1);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcAwake)
% hold on
% before = findall(gca);
% fnplt(curve2);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcNREM)
% before = findall(gca);
% fnplt(curve3);
% added = setdiff(findall(gca),before);
% set(added,'Color',colorRfcREM)
% title({'[2g] Variance of whisker angle','arousal-state distribution'})
% xlabel('Whisker angle (deg^2)')
% ylabel('Probability')
% xlim([-3,3])
% ylim([0,0.35])
% axis square
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureB,[dirpath 'Fig2_B']);
    set(summaryFigureB,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig2_B'])
end


