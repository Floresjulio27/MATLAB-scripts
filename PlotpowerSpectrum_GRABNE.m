function [AnalysisResults] = PlotpowerSpectrum_GRABNE(rootFolder,saveFigs,delim,AnalysisResults)

animalIDs = {'GRABNE001'} ;% fieldnames(AnalysisResults);

behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
dataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_gammaBandPower'};
% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.PowerSpec.(behavField).(dataType).S = [];
        data.PowerSpec.(behavField).(dataType).f = [];
    end
end
% concatenate power spectra during different arousal states for each animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).S) == false
                data.PowerSpec.(behavField).(dataType).S = cat(2,data.PowerSpec.(behavField).(dataType).S,AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).S);
                data.PowerSpec.(behavField).(dataType).f = cat(1,data.PowerSpec.(behavField).(dataType).f,AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).f);
            end
        end
    end
end
% mean and standard error for arousal state power
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.PowerSpec.(behavField).(dataType).meanS = mean(data.PowerSpec.(behavField).(dataType).S,2);
        data.PowerSpec.(behavField).(dataType).semS = std(data.PowerSpec.(behavField).(dataType).S,0,2)./sqrt(size(data.PowerSpec.(behavField).(dataType).S,2));
        data.PowerSpec.(behavField).(dataType).meanf = mean(data.PowerSpec.(behavField).(dataType).f,1);
    end
end
%% power spectrum during arousal states
    powerSpec_Fig = figure('Name','PowerSpectrum_Figure');
    sgtitle('Power Spectrum in different arousal stages')
       SubplotNeed = ceil(length(dataTypes)/3);
       splotNo = 1;
for PPS = 1:length(dataTypes)
    dataType = dataTypes{1,PPS};

    subplot(3,SubplotNeed,splotNo);
    splotNo = splotNo + 1;

    L1 = loglog(data.PowerSpec.Rest.(dataType).meanf,data.PowerSpec.Rest.(dataType).meanS,'color','green','LineWidth',2);
    hold on
    loglog(data.PowerSpec.Rest.(dataType).meanf,data.PowerSpec.Rest.(dataType).meanS + data.PowerSpec.Rest.(dataType).semS,'color','green','LineWidth',0.5);
    loglog(data.PowerSpec.Rest.(dataType).meanf,data.PowerSpec.Rest.(dataType).meanS - data.PowerSpec.Rest.(dataType).semS,'color','green','LineWidth',0.5);
    rectangle('Position',[0.004,0.00004,0.1 - 0.004,0.001],'FaceColor','w','EdgeColor','w')
    
    L2 = loglog(data.PowerSpec.NREM.(dataType).meanf,data.PowerSpec.NREM.(dataType).meanS,'color','black','LineWidth',2);
    loglog(data.PowerSpec.NREM.(dataType).meanf,data.PowerSpec.NREM.(dataType).meanS + data.PowerSpec.NREM.(dataType).semS,'color','black','LineWidth',0.5);
    loglog(data.PowerSpec.NREM.(dataType).meanf,data.PowerSpec.NREM.(dataType).meanS - data.PowerSpec.NREM.(dataType).semS,'color','black','LineWidth',0.5);
%     rectangle('Position',[0.004,0.00004,1/30 - 0.004,0.001],'FaceColor','w','EdgeColor','w')
    
    L3 = loglog(data.PowerSpec.REM.(dataType).meanf,data.PowerSpec.REM.(dataType).meanS,'color','red','LineWidth',2);
    loglog(data.PowerSpec.REM.(dataType).meanf,data.PowerSpec.REM.(dataType).meanS + data.PowerSpec.REM.(dataType).semS,'color','red','LineWidth',0.5);
    loglog(data.PowerSpec.REM.(dataType).meanf,data.PowerSpec.REM.(dataType).meanS - data.PowerSpec.REM.(dataType).semS,'color','red','LineWidth',0.5);
    rectangle('Position',[0.004,0.00004,1/60 - 0.004,0.001],'FaceColor','w','EdgeColor','w')
    
%     L4 = loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS,'color',[0.4940 0.1840 0.5560],'LineWidth',2);
%     loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS + data.PowerSpec.Awake.(dataType).semS,'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
%     loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS - data.PowerSpec.Awake.(dataType).semS,'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
    
%     L5 = loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS,'color','magenta','LineWidth',2);
%     loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS + data.PowerSpec.Asleep.(dataType).semS,'color','magenta','LineWidth',0.5);
%     loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS - data.PowerSpec.Asleep.(dataType).semS,'color','magenta','LineWidth',0.5);
    
    L6 = loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS,'color','blue','LineWidth',2);
    loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS + data.PowerSpec.All.(dataType).semS,'color','blue','LineWidth',0.5);
    loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS - data.PowerSpec.All.(dataType).semS,'color','blue','LineWidth',0.5);
    xline(1/10)
    xline(1/30)
    xline(1/60)
    title([dataType])
    ylabel('power (a.u.)')
    xlabel('Freq (Hz)')
    if splotNo == length(dataTypes) + 1
    legend([L1,L2,L3,L6],'Rest','NREM','REM','All','Location','best')
    end
    axis square
    xlim([0.003,1])
%     ylim([0.000035,0.0012])
    set(gca,'box','off')
%     ax4.TickLength = [0.03,0.03];
end
    %% save the figure
    if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(powerSpec_Fig,[dirpath animalID  'Fig_PowerSpectrum']);
        set(powerSpec_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath animalID 'Fig_PowerSpectrum'])
    end
    close
