function [AnalysisResults] = PlotCOherence_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs)
%% coherence

animalIDs = {'GRABNE001','GRABNE002'} ; %fieldnames(AnalysisResults);

    if firstHrs == "false"
         behavFields = {'Rest','NREM','REM','Awake'};
    elseif firstHrs == "true"
        behavFields = {'Rest','NREM','Awake'};
    end

dataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};
SdataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};

% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
                data.Coherr.(behavField).(dataType).(CName) = [];
                data.Coherr.(behavField).(dataType).(fName) = [];
%                 data.Coherr.(behavField).(dataType).animalID = {};
                data.Coherr.(behavField).(dataType).behavField = {};
        end
    end
end
% concatenate coherence during different arousal states for each animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
            % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).C) == false
                    data.Coherr.(behavField).(dataType).(CName) = cat(2,data.Coherr.(behavField).(dataType).(CName),AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).C);
                    data.Coherr.(behavField).(dataType).(fName) = cat(1,data.Coherr.(behavField).(dataType).(fName),AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).f);
%                     data.Coherr.(behavField).(dataType).animalID = cat(1,data.Coherr.(behavField).(dataType).animalID,animalID,animalID);
                    data.Coherr.(behavField).(dataType).behavField = cat(1,data.Coherr.(behavField).(dataType).behavField,behavField,behavField);
                end
            end
        end
    end
end
% mean and standard error for arousal state coherence
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
                meanCName = ['mean' SdataType 'C'];
                semCName = ['sem' SdataType 'C'];
                meanfName = ['mean' SdataType 'f'];
                data.Coherr.(behavField).(dataType).(meanCName) = mean(data.Coherr.(behavField).(dataType).(CName),2);
                data.Coherr.(behavField).(dataType).(semCName) = std(data.Coherr.(behavField).(dataType).(CName),0,2)./sqrt(size(data.Coherr.(behavField).(dataType).(CName),2));
                data.Coherr.(behavField).(dataType).(meanfName) = mean(data.Coherr.(behavField).(dataType).(fName),1);
        end
    end
end
%% plot coherence during arousal states

 for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        Coherence_Fig = figure('Name','Coherence_Figure');
        sgtitle([(dataType) 'Coherence in different arousal stages'])
        SubplotNeed = ceil(length(SdataTypes)/3);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
            meanCName = ['mean' SdataType 'C'];
            semCName = ['sem' SdataType 'C'];
            meanfName = ['mean' SdataType 'f'];
         subplot(3,SubplotNeed,splotNo);

            L1 = semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName),'color','green','LineWidth',2);
                    hold on
                    semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName) + data.Coherr.Rest.(dataType).(semCName),'color','green','LineWidth',0.5);
                    semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName) - data.Coherr.Rest.(dataType).(semCName),'color','green','LineWidth',0.5);
                    rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
            L2 = semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName),'color','black','LineWidth',2);
                    hold on
                    semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName) + data.Coherr.NREM.(dataType).(semCName),'color','black','LineWidth',0.5);
                    semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName) - data.Coherr.NREM.(dataType).(semCName),'color','black','LineWidth',0.5);
                    rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
           if firstHrs == "false"
            L3 = semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName),'color','red','LineWidth',2);
                hold on        
                semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName) + data.Coherr.REM.(dataType).(semCName),'color','red','LineWidth',0.5);
                semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName) - data.Coherr.REM.(dataType).(semCName),'color','red','LineWidth',0.5);
                rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
           elseif firstHrs == "true"
           end
%             semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName),'color',[0.4940 0.1840 0.5560],'LineWidth',2);
%             semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName) + data.Coherr.Awake.(dataType).(semCName),'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
%             semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName) - data.Coherr.Awake.(dataType).(semCName),'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
%     %         semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName),'color',colors('custom asleep'),'LineWidth',2);
    %         semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName) + data.Coherr.Asleep.(dataType).(semCName),'color',colors('custom asleep'),'LineWidth',0.5);
    %         semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName) - data.Coherr.Asleep.(dataType).(semCName),'color',colors('custom asleep'),'LineWidth',0.5);
    %         semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName),'color',colors('custom all'),'LineWidth',2);
    %         semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName) + data.Coherr.All.(dataType).(semCName),'color',colors('custom all'),'LineWidth',0.5);
    %         semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName) - data.Coherr.All.(dataType).(semCName),'color',colors('custom all'),'LineWidth',0.5);
            xline(0.02,'color','b');
            xline(0.35,'color','r');
            title({ (dataType) },...
                {[(SdataType) ' coherence']})
            ylabel('Coherence')
            xlabel('Freq (Hz)')
            if strcmp(dataType,SdataType) == 1
                if firstHrs == "false"
                    legend([L1,L2,L3],'Rest','NREM','REM','Location','best')
                elseif firstHrs == "true"
                    legend([L1,L2],'Rest','NREM','Location','best')
                end
            end
            axis square
            xlim([0.001,1])
            ylim([0,0.5])
            set(gca,'box','off')
            splotNo = splotNo +1;
        end
        %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Coherence_Fig,[dirpath 'Fig_Coherence_' (dataType)]);
        set(Coherence_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_Coherence_' (dataType)])
    end
    close
 end

