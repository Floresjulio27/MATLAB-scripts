% close all;
% clear all; 
% clc
function TDTFiberPhotometry_bilateralGCaMP7s_OnlyHemodynamicCorrect(filepath,rawDataFilespath)

    for trialNum = 1:1:length(rawDataFilespath)
            %% generate file save path for each trials
            load(rawDataFilespath{trialNum,1});
            Params = FiberData.params;

            %reload motion corrected signals
            % remove some data to perform a better fit
            timeN = 210; % remove first 210 seconds;
            sampleN = floor(timeN*Params.DataFs);


            motionCorrect_LH.F405 = FiberData.LH.motionCorrect.F405;
            motionCorrect_LH.F465 = FiberData.LH.motionCorrect.F465(sampleN+2:end-sampleN-1);
            motionCorrect_LH.F560 = FiberData.LH.motionCorrect.F560(sampleN+2:end-sampleN-1);

            motionCorrect_RH.F405 = FiberData.RH.motionCorrect.F405;
            motionCorrect_RH.F465 = FiberData.RH.motionCorrect.F465(sampleN+2:end-sampleN-1);
            motionCorrect_RH.F560 = FiberData.RH.motionCorrect.F560(sampleN+2:end-sampleN-1);
            
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
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,motionCorrect_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals LH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,motionCorrect_RH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals RH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');
            linkaxes(h,'x');
            if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
            end
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.tiff'],'tiff')
            close 
            %% Correct 465 and 560 using 405
            [dF_LH] =  control405Correction(CBVCorrection_LH,Params,'LH');
            [dF_RH] =  control405Correction(CBVCorrection_RH,Params,'RH');
            %% house keeping
            timeN = 210; % readd data that was removed before . seconds;
            sampleN = floor(timeN*Params.DataFs);
          
           % add some dummy values to match the index

            fieldNames = {'F465','F560'};
            for corI=1:1:length(fieldNames)
                    
                    CBVCorrection_LH.(fieldNames{corI}) = [CBVCorrection_LH.(fieldNames{corI})(1:sampleN); CBVCorrection_LH.(fieldNames{corI})(1:end); CBVCorrection_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    CBVCorrection_RH.(fieldNames{corI}) = [CBVCorrection_RH.(fieldNames{corI})(1:sampleN); CBVCorrection_RH.(fieldNames{corI})(1:end); CBVCorrection_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    dF_LH.(fieldNames{corI}) = [dF_LH.(fieldNames{corI})(1:sampleN); dF_LH.(fieldNames{corI})(1:end); dF_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    dF_RH.(fieldNames{corI}) = [dF_RH.(fieldNames{corI})(1:sampleN); dF_RH.(fieldNames{corI})(1:end); dF_RH.(fieldNames{corI})(end-sampleN-1:end)];

            end
            %% save the data            
            FiberData.LH.TRITCCorrection = CBVCorrection_LH;
            FiberData.RH.TRITCCorrection = CBVCorrection_RH;

            FiberData.LH.dF = dF_LH;
            FiberData.RH.dF = dF_RH;
      
            FiberData.params = Params;

            clearvars -except -regexp Data$ path$ Params$ trialNum$
            savepathmat = [rawDataFilespath{trialNum}];
            save(savepathmat,"FiberData",'-v7.3')
    end