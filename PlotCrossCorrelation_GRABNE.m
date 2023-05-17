function [AnalysisResults] = PlotCrossCorrelation_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs)
%% cross correlation

animalIDs = {'GRABNE001','GRABNE002'} ; %fieldnames(AnalysisResults);

behavFields = {'Rest','NREM','REM','Asleep','All'};%, 'Alert'};

    if firstHrs == "false"
        behavFields = {'Rest','NREM','REM','Asleep','All'};
    elseif firstHrs == "true"
        behavFields = {'Rest','NREM','Asleep','All'};
    end

SdataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};
dataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};

% concatenate the cross-correlation during different arousal states for each animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                lagsName = [SdataType 'lags'];
                xcValsName = [SdataType 'xcVals'];

                        data.XCorr.(behavField).(dataType).dummCheck = 1;
                        if isfield(data.XCorr.(behavField).(dataType),(xcValsName)) == false
                            data.XCorr.(behavField).(dataType).(xcValsName) = [];
                            % lags and stats fields
                            data.XCorr.(behavField).(dataType).(lagsName) = [];
                            data.XCorr.(behavField).(dataType).animalID = {};
                            data.XCorr.(behavField).(dataType).behavField = {};
                        end
                        % concatenate cross correlation during each arousal state, find peak + lag time
                        if isfield(AnalysisResults.(animalID).CrossCorrelation,(behavField)) == true
                            if isempty(AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).xcVals) == false
                                % peak + lag time
                                data.XCorr.(behavField).(dataType).(xcValsName) = cat(1,data.XCorr.(behavField).(dataType).(xcValsName),AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).xcVals);
                                % lags and stats fields
                                data.XCorr.(behavField).(dataType).(lagsName) = cat(1,data.XCorr.(behavField).(dataType).(lagsName),AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).lags);
%                                 data.XCorr.(behavField).(dataType).animalID = cat(1,data.XCorr.(behavField).(dataType).animalID,animalID);
                                data.XCorr.(behavField).(dataType).behavField = cat(1,data.XCorr.(behavField).(dataType).behavField,behavField);
                            end
                        end
            end
        end
    end
end
samplingRate = 30;
% mean and standard error/standard deviation of cross correlation values
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
         for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                lagsName = [SdataType 'lags'];
                xcValsName = [SdataType 'xcVals'];
                meanlagsName = ['mean' SdataType 'lags'];
                meanxcValsName = ['mean' SdataType 'xcVals'];
                stdxcValsName = ['std' SdataType 'xcVals'];
                semxcValsName = ['sem' SdataType 'xcVals'];
                peakxcValsNameVal = ['peak' SdataType 'Val'];
                peakxcValsNameStd = ['peak' SdataType 'Std'];
                peakxcValsNameLag = ['peak' SdataType 'Lag'];

                % Lags time vector
                data.XCorr.(behavField).(dataType).(meanlagsName) = mean(data.XCorr.(behavField).(dataType).(lagsName),1);
                % XC mean/sem
                data.XCorr.(behavField).(dataType).(meanxcValsName) = mean(data.XCorr.(behavField).(dataType).(xcValsName),1);
                data.XCorr.(behavField).(dataType).(stdxcValsName) = std(data.XCorr.(behavField).(dataType).(xcValsName),0,1);
                data.XCorr.(behavField).(dataType).(semxcValsName) = std(data.XCorr.(behavField).(dataType).(xcValsName),0,1)./sqrt(size(data.XCorr.(behavField).(dataType).(xcValsName),1));
                % find peak lag time and value/std at that time
                [~,idx] = max(abs(data.XCorr.(behavField).(dataType).(meanxcValsName)));
                data.XCorr.(behavField).(dataType).(peakxcValsNameVal) = data.XCorr.(behavField).(dataType).((meanxcValsName))(1,idx);
                data.XCorr.(behavField).(dataType).(peakxcValsNameStd) = data.XCorr.(behavField).(dataType).((stdxcValsName))(1,idx);
                data.XCorr.(behavField).(dataType).(peakxcValsNameLag) = data.XCorr.(behavField).(dataType).((meanlagsName))(1,idx)/samplingRate;
         end
    end
end
%% cross correlation [rest, NREM, REM]

for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        Crosscorrelation_Fig = figure('Name','CrossCorrelation_Figure');
        sgtitle([(dataType) 'CrossCorrelation in different arousal stages'])
        SubplotNeed = ceil(length(SdataTypes)/3);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
                meanlagsName = ['mean' SdataType 'lags'];
                meanxcValsName = ['mean' SdataType 'xcVals'];
%                 stdxcValsName = ['std' SdataType 'xcVals'];
                semxcValsName = ['sem' SdataType 'xcVals'];
%                 peakxcValsNameVal = ['peak' SdataType 'Val'];
%                 peakxcValsNameStd = ['peak' SdataType 'Std'];
%                 peakxcValsNameLag = ['peak' SdataType 'Lag'];
            
            subplot(3,SubplotNeed,splotNo);
        
            freq = 30;
            lagSec = 5;
            L1 = plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName),'color','green','LineWidth',2);
            hold on
            plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName) + data.XCorr.Rest.(dataType).(semxcValsName),'color','green','LineWidth',0.5);
            plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName) - data.XCorr.Rest.(dataType).(semxcValsName),'color','green','LineWidth',0.5);
            
            L2 = plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName),'color','black','LineWidth',2);
            plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName) + data.XCorr.NREM.(dataType).(semxcValsName),'color','black','LineWidth',0.5);
            plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName) - data.XCorr.NREM.(dataType).(semxcValsName),'color','black','LineWidth',0.5);
            if firstHrs == "false"
                L3= plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName),'color','red','LineWidth',2);
                plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName) + data.XCorr.REM.(dataType).(semxcValsName),'color','red','LineWidth',0.5);
                plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName) - data.XCorr.REM.(dataType).(semxcValsName),'color','red','LineWidth',0.5);
            elseif firstHrs == "true"
            end
%             L4 = plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName),'color','c','LineWidth',2);
%             plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName) + data.XCorr.Asleep.(dataType).(semxcValsName),'color','c','LineWidth',0.5);
%             plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName) - data.XCorr.Asleep.(dataType).(semxcValsName),'color','c','LineWidth',0.5);
%             
            L5= plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName),'color','magenta','LineWidth',2);
            plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName) + data.XCorr.All.(dataType).(semxcValsName),'color','magenta','LineWidth',0.5);
            plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName) - data.XCorr.All.(dataType).(semxcValsName),'color','magenta','LineWidth',0.5);
            
            
            xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
            xticklabels({'-5','-2.5','0','2.5','5'})
            xlim([-lagSec*freq,lagSec*freq])
%             ylim([-0.35,0.15])
            xlabel('Lags (s)')
            ylabel('Correlation')
            title({[(dataType) ' '],...
                [(SdataType) ' XCorr']})
            if strcmp(SdataType,dataType) == 1
                 if firstHrs == "false"
                      legend([L1,L2,L3,L5],'Rest','NREM','REM','All','Location','best')


                elseif firstHrs == "true"
                      legend([L1,L2,L5],'Rest','NREM','All','Location','best')
                end   
            end
            axis square
            set(gca,'box','off')

            splotNo = splotNo +1;
        end
        %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Crosscorrelation_Fig,[dirpath 'Fig_Crosscorrelation_' (dataType)]);
        set(Crosscorrelation_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_Crosscorrelation_' (dataType)])
    end
    close
 end


