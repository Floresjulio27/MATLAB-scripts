% close all;
% clear all; 
% clc
function TDTFiberPhotometry_bilateralGCaMP7s_Rerun(filepath,rawDataFilespath)

    for trialNum = 1:1:length(rawDataFilespath)
            %% generate file save path for each trials
            load(rawDataFilespath{trialNum,1});
            Params = FiberData.params;
            Params.decay_Freq=0.005;% 0.001
            ForceSensor = FiberData.pressureSensor;
            rawDataLH = [FiberData.LH.rawData.F405,FiberData.LH.rawData.F465,FiberData.LH.rawData.F560];
            rawDataRH = [FiberData.RH.rawData.F405,FiberData.RH.rawData.F465,FiberData.RH.rawData.F560];
            %% Correct exponential decay
            [expCorrected_LH, predicted_LH] = Correct_ExponentialDecay_HighPass( rawDataLH,ForceSensor,Params,'LH');
            [expCorrected_RH, predicted_RH] = Correct_ExponentialDecay_HighPass( rawDataRH,ForceSensor,Params,'RH');
            %% Low pass the data 
            lowPassData_LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F405);
            lowPassData_LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F465);
            lowPassData_LH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F560);
            lowPassData_RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F405);
            lowPassData_RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F465);
            lowPassData_RH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F560);
            %% rescale the data from 0 to 1
            for qN=1:size(lowPassData_LH,2)
            rescaleData_LH.F405(:,qN)=rescale(lowPassData_LH.F405(:,qN),0,1); %rescale all data between 0 to 1
            rescaleData_LH.F465(:,qN)=rescale(lowPassData_LH.F465(:,qN),0,1); %rescale all data between 0 to 1
            rescaleData_LH.F560(:,qN)=rescale(lowPassData_LH.F560(:,qN),0,1); %rescale all data between 0 to 1
            end
            for pN=1:size(lowPassData_RH,2)
            rescaleData_RH.F405(:,pN)=rescale(lowPassData_RH.F405(:,pN),0,1); %rescale all data between 0 to 1
            rescaleData_RH.F465(:,pN)=rescale(lowPassData_RH.F465(:,pN),0,1); %rescale all data between 0 to 1
            rescaleData_RH.F560(:,pN)=rescale(lowPassData_RH.F560(:,pN),0,1); %rescale all data between 0 to 1
            end
            %% remove motion related signal
            [motionCorrect_LH] =  Correct_motionartifact(rescaleData_LH,Params);  
            [motionCorrect_RH] =  Correct_motionartifact(rescaleData_RH,Params); 
            %% hemodynamic Corrections
            Params.correctionConstant = -0.23; % used from Kyles Data -0.23
            CBVCorrection_LH = motionCorrect_LH;
            CBVCorrection_RH = motionCorrect_RH;
            CBVCorrection_LH.F465 = CBVCorrection_LH.F465-(Params.correctionConstant*CBVCorrection_LH.F560);
            CBVCorrection_RH.F465 = CBVCorrection_RH.F465-(Params.correctionConstant*CBVCorrection_RH.F560);
            %
            figTime=(1:length(CBVCorrection_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals LH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals RH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.tiff'],'tiff')
            close 
            %% Correct 465 and 560 using 405
            [dF_LH] =  control405Correction(CBVCorrection_LH,Params,'LH');
            [dF_RH] =  control405Correction(CBVCorrection_RH,Params,'RH');
            %% plot Data comparison
            figTime=(1:length(expCorrected_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F465),'constant')); 
            title('Exp Corrected GCaMP7s Signals'); xlabel('Time (min)'); legend({'LH GCaMP7s','RH GCaMP7s'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F560),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F560),'constant')); 
            title('Exp Corrected TRITC Signals'); xlabel('Time (min)'); legend({'LH TRITC','RH TRITC'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'ExpCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'ExpCorrectionCompare.tiff'],'tiff')
            close 
            
            %
            figTime=(1:length(dF_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_RH.F465),'constant')); 
            title('Isosbestic Correction GCaMP7s Signals'); xlabel('Time (min)'); legend({'LH GCaMP7s','RH GCaMP7s'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_LH.F560),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_RH.F560),'constant')); 
            title('Isosbestic Correction TRITC Signals'); xlabel('Time (min)'); legend({'LH TRITC','RH TRITC'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrectionCompare.tiff'],'tiff')
            close 
            %% house keeping
            timeN = 210; % remove first 15 seconds;
            sampleN = floor(timeN*Params.DataFs);
          
           % add some dummy values to match the index
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    predicted_LH.(fieldNames{corI}) = [predicted_LH.(fieldNames{corI})(1:sampleN); predicted_LH.(fieldNames{corI});predicted_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    predicted_RH.(fieldNames{corI}) = [predicted_RH.(fieldNames{corI})(1:sampleN); predicted_RH.(fieldNames{corI});predicted_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    expCorrected_LH.(fieldNames{corI}) = [expCorrected_LH.(fieldNames{corI})(1:sampleN); expCorrected_LH.(fieldNames{corI})(1:end); expCorrected_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    expCorrected_RH.(fieldNames{corI}) = [expCorrected_RH.(fieldNames{corI})(1:sampleN); expCorrected_RH.(fieldNames{corI})(1:end); expCorrected_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    lowPassData_LH.(fieldNames{corI}) = [lowPassData_LH.(fieldNames{corI})(1:sampleN); lowPassData_LH.(fieldNames{corI})(1:end); lowPassData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    lowPassData_RH.(fieldNames{corI}) = [lowPassData_RH.(fieldNames{corI})(1:sampleN); lowPassData_RH.(fieldNames{corI})(1:end); lowPassData_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    rescaleData_LH.(fieldNames{corI}) = [rescaleData_LH.(fieldNames{corI})(1:sampleN); rescaleData_LH.(fieldNames{corI})(1:end); rescaleData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    rescaleData_RH.(fieldNames{corI}) = [rescaleData_RH.(fieldNames{corI})(1:sampleN); rescaleData_RH.(fieldNames{corI})(1:end); rescaleData_RH.(fieldNames{corI})(end-sampleN-1:end)];
            end
            
            fieldNames = {'F465','F560'};
            for corI=1:1:length(fieldNames)
                    motionCorrect_LH.(fieldNames{corI}) = [motionCorrect_LH.(fieldNames{corI})(1:sampleN); motionCorrect_LH.(fieldNames{corI})(1:end); motionCorrect_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    motionCorrect_RH.(fieldNames{corI}) = [motionCorrect_RH.(fieldNames{corI})(1:sampleN); motionCorrect_RH.(fieldNames{corI})(1:end); motionCorrect_RH.(fieldNames{corI})(end-sampleN-1:end)];
                    
                    CBVCorrection_LH.(fieldNames{corI}) = [CBVCorrection_LH.(fieldNames{corI})(1:sampleN); CBVCorrection_LH.(fieldNames{corI})(1:end); CBVCorrection_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    CBVCorrection_RH.(fieldNames{corI}) = [CBVCorrection_RH.(fieldNames{corI})(1:sampleN); CBVCorrection_RH.(fieldNames{corI})(1:end); CBVCorrection_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    dF_LH.(fieldNames{corI}) = [dF_LH.(fieldNames{corI})(1:sampleN); dF_LH.(fieldNames{corI})(1:end); dF_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    dF_RH.(fieldNames{corI}) = [dF_RH.(fieldNames{corI})(1:sampleN); dF_RH.(fieldNames{corI})(1:end); dF_RH.(fieldNames{corI})(end-sampleN-1:end)];

            end
            %% save the data            

            FiberData.LH.expCorrected = expCorrected_LH;
            FiberData.RH.expCorrected = expCorrected_RH;
            FiberData.LH.predicted = predicted_LH;
            FiberData.RH.predicted = predicted_RH;
            FiberData.LH.lowPassData = lowPassData_LH;
            FiberData.RH.lowPassData = lowPassData_RH;
            FiberData.LH.rescaleData= rescaleData_LH;
            FiberData.RH.rescaleData = rescaleData_RH;
            FiberData.LH.motionCorrect= motionCorrect_LH;
            FiberData.RH.motionCorrect = motionCorrect_RH;
            FiberData.LH.TRITCCorrection = CBVCorrection_LH;
            FiberData.RH.TRITCCorrection = CBVCorrection_RH;

            FiberData.LH.dF = dF_LH;
            FiberData.RH.dF = dF_RH;
%             Params.outliersPos.LH = TF_Outliers_LH;
%             Params.outliersPos.RH = TF_Outliers_RH;
            FiberData.params = Params;
            clearvars -except -regexp Data$ path$ Params$ trialNum$
            savepathmat = [rawDataFilespath{trialNum}];
            save(savepathmat,"FiberData",'-v7.3')
    end