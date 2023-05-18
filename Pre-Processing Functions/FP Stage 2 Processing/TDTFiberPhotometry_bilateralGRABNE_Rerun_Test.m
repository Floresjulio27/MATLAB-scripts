function TDTFiberPhotometry_bilateralGRABNE_Rerun_Test(fiberDataFilespath)

    
    for trialNum = 1:1:length(fiberDataFilespath)
            %% generate file save path for each trials
            Fiberpath = [fiberDataFilespath{trialNum}];
            load(Fiberpath);
            Params = FiberData.params;
            %% 
            if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
            end 
            %%
            OR_RawData_Ach(:,1) = FiberData.Ach.rawData.F405;
            OR_RawData_Ach(:,2) = FiberData.Ach.rawData.F465;
            OR_RawData_Ach(:,3) = FiberData.Ach.rawData.F560;
            OR_RawData_NE(:,1) = FiberData.NE.rawData.F405;
            OR_RawData_NE(:,2) = FiberData.NE.rawData.F465;
            OR_RawData_NE(:,3) = FiberData.NE.rawData.F560;
            
            ForceSensor = FiberData.pressureSensor;
            %% Correct exponential decay
            [expCorrected_RH, ~,Params] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_NE,ForceSensor,Params,'NE');
            [expCorrected_LH, ~,Params] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_Ach,ForceSensor,Params,'Ach');
            %% Low pass the data 
            lowPassData_LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F405);
            lowPassData_LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F465);
            lowPassData_LH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F560);
            lowPassData_RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F405);
            lowPassData_RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F465);
            lowPassData_RH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F560);
            %% zScore the signal
            % determine the resting state baseline to calculate ZScore
            if isfield(Params,'baselineStartTime') == false 
                %{
                FrcSensor = ForceSensor((210*Params.DataFs)+1:end-(210*Params.DataFs));
                ImageCheck = figure;
                ImageCheck.WindowState = 'maximized';
                drawnow;
                L(1) = subplot(3,1,1);
                plot((1:length(FrcSensor))/Params.DataFs,FrcSensor);
                hold on 
                yyaxis right
                plot((1:length(lowPassData_RH.F405))/Params.DataFs,lowPassData_RH.F405);
                legend('Force Sensor', 'RH 405')
    
                xlim([0 length(FrcSensor)/Params.DataFs]);
                ylabel('forceSensor');
                L(2) = subplot(3,1,2);
                plot((1:length(lowPassData_LH.F560))/Params.DataFs,lowPassData_LH.F560);
                ylabel('dF LH 560');            
                xlim([0 length(lowPassData_RH.F560)/Params.DataFs]);
    
                L(3) = subplot(3,1,3);
                plot((1:length(lowPassData_RH.F560))/Params.DataFs,lowPassData_RH.F560);
                ylabel('dF RH 560');
                xlim([0 length(lowPassData_RH.F560)/Params.DataFs]);
                linkaxes(L,'x');
                commandwindow;
                Params.baselineStartTime = input('Input the start time for resting data: '); disp(' ')
                Params.baselineEndTime = input('Input the end time for resting data: '); disp(' ')
                close(ImageCheck) 
                %}
                Params.baselineStartTime = 500;
                Params.baselineEndTime = 2500;
            end

            [zScored_LH] =  ZScoreFiberData(lowPassData_LH,Params,'LH');
            [zScored_RH] =  ZScoreFiberData(lowPassData_RH,Params,'RH');
            %% remove motion related signal
            [motionCorrect_LH] =  Correct_motionartifact(zScored_LH,Params);  
            [motionCorrect_RH] =  Correct_motionartifact(zScored_RH,Params); 
            %% remove isosbestic changes related signal
            [isosCorrect_LH] = control405Correction_Update(motionCorrect_LH,Params,'LH');
            [isosCorrect_RH] = control405Correction_Update(motionCorrect_RH,Params,'RH');
            %% hemodynamic Corrections
            [Params.CorrectSlope_LH , Corrected_F465_LH ] = HemodynamicsCalculation(isosCorrect_LH,Params,Fiberpath,'LH');
            [Params.CorrectSlope_RH , Corrected_F465_RH ] = HemodynamicsCalculation(isosCorrect_RH,Params,Fiberpath,'RH');
%{
%             CorrectSlope_RH = HemodynamicsConstantCalculation(lowPassData_RH,Fiberpath,'RH'); % calculate the slope to correct for CBV
%             CorrectSlope_RH = HemodynamicsConstantCalculation(isosCorrect_RH,Fiberpath,'RH'); % calculate the slope to correct for CBV
%             CorrectSlope_LH = HemodynamicsConstantCalculation(isosCorrect_LH,Fiberpath); % calculate the slope to correct for CBV

            disp(['The hemodynamic correction slope is = ' num2str(CorrectSlope_RH)]); disp(' ')
            Params.correctionConstant_calculated = CorrectSlope_RH;

            ChangeSlope = input('Do you want to change the slope (y/n)? ','s'); disp(' ')
            if ChangeSlope == ("y")
                CorrectSlope_RH = input('What is the new slope? '); %-0.23;
            end

            Params.correctionConstant_used = CorrectSlope_RH;
%}
            CBVCorrection_LH = isosCorrect_LH;
            CBVCorrection_RH = isosCorrect_RH;

            CBVCorrection_LH.F465 = Corrected_F465_LH;
            CBVCorrection_RH.F465 = Corrected_F465_RH;

   
            figTime=(1:length(CBVCorrection_LH.F465))/(Params.DataFs*60);

            ImageCheck = figure;
            ImageCheck.WindowState = 'minimized';

            k(1) = subplot(4,1,1);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_RH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_RH.F560))
            legend('RH 465','RH 560')
            title('Pre-Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\delta F rescaled');

            k(2) = subplot(4,1,3);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_LH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_LH.F560))
            legend('LH 465','LH 560')
            title('Before Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\delta F rescaled');

%             CBVCorrection_LH.F465 = CBVCorrection_LH.F465-(Params.correctionConstant*CBVCorrection_LH.F560);
%             CBVCorrection_RH.F465 = CBVCorrection_RH.F465-(Params.correctionConstant*CBVCorrection_RH.F560);

            k(3)=subplot(4,1,2); 
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F560))
            legend('RH 465','RH 560')
            title('After Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F rescaled');

            k(4)=subplot(4,1,4); 
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F560))
            legend('LH 465','LH 560')
            title('After Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F rescaled');

            linkaxes(k,'x');

            saveas(gcf,['../Figures/Corrections/' Params.savepath 'HemoCorrection.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'HemoCorrection.tiff'],'tiff')
            close(ImageCheck)


            %% plot Data comparison
         
            %{
            figTime=(1:length(dFF0_z_LH.F465))/(Params.DataFs*60);
            ImageCheck = figure;
            ImageCheck.WindowState = 'minimized';
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dFF0_z_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dFF0_z_RH.F465),'constant')); 
            title('Final GFP Signals'); xlabel('Time (min)'); legend({'LH Ach','RH NE'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dFF0_z_LH.F560),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dFF0_z_RH.F560),'constant')); 
            title('Final Rhodamine Signals'); xlabel('Time (min)'); legend({'LH Rhodamine','RH Rhodamine'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'FinalSignalCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'FinalSignalCompare.tiff'],'tiff')
            close(ImageCheck) 
            %}
            %% house keeping
            timeN = 210; % remove first 15 seconds;
            sampleN = floor(timeN*Params.DataFs);
          
           % add some dummy values to match the index
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)

                    expCorrected_LH.(fieldNames{corI}) = [expCorrected_LH.(fieldNames{corI})(1:sampleN); expCorrected_LH.(fieldNames{corI})(1:end); expCorrected_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    expCorrected_RH.(fieldNames{corI}) = [expCorrected_RH.(fieldNames{corI})(1:sampleN); expCorrected_RH.(fieldNames{corI})(1:end); expCorrected_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    lowPassData_LH.(fieldNames{corI}) = [lowPassData_LH.(fieldNames{corI})(1:sampleN); lowPassData_LH.(fieldNames{corI})(1:end); lowPassData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    lowPassData_RH.(fieldNames{corI}) = [lowPassData_RH.(fieldNames{corI})(1:sampleN); lowPassData_RH.(fieldNames{corI})(1:end); lowPassData_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    zScored_LH.(fieldNames{corI}) = [zScored_LH.(fieldNames{corI})(1:sampleN); zScored_LH.(fieldNames{corI})(1:end); zScored_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    zScored_RH.(fieldNames{corI}) = [zScored_RH.(fieldNames{corI})(1:sampleN); zScored_RH.(fieldNames{corI})(1:end); zScored_RH.(fieldNames{corI})(end-sampleN-1:end)];
            end
            
            fieldNames = {'F465','F560'};
            for corI=1:1:length(fieldNames)
                    motionCorrect_LH.(fieldNames{corI}) = [motionCorrect_LH.(fieldNames{corI})(1:sampleN); motionCorrect_LH.(fieldNames{corI})(1:end); motionCorrect_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    motionCorrect_RH.(fieldNames{corI}) = [motionCorrect_RH.(fieldNames{corI})(1:sampleN); motionCorrect_RH.(fieldNames{corI})(1:end); motionCorrect_RH.(fieldNames{corI})(end-sampleN-1:end)];
                    
                    CBVCorrection_LH.(fieldNames{corI}) = [CBVCorrection_LH.(fieldNames{corI})(1:sampleN); CBVCorrection_LH.(fieldNames{corI})(1:end); CBVCorrection_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    CBVCorrection_RH.(fieldNames{corI}) = [CBVCorrection_RH.(fieldNames{corI})(1:sampleN); CBVCorrection_RH.(fieldNames{corI})(1:end); CBVCorrection_RH.(fieldNames{corI})(end-sampleN-1:end)];

%                     dF_LH.(fieldNames{corI}) = [dF_LH.(fieldNames{corI})(1:sampleN); dF_LH.(fieldNames{corI})(1:end); dF_LH.(fieldNames{corI})(end-sampleN-1:end)];
%                     dF_RH.(fieldNames{corI}) = [dF_RH.(fieldNames{corI})(1:sampleN); dF_RH.(fieldNames{corI})(1:end); dF_RH.(fieldNames{corI})(end-sampleN-1:end)];

%                     dFF0_z_LH.(fieldNames{corI}) = [dFF0_z_LH.(fieldNames{corI})(1:sampleN); dFF0_z_LH.(fieldNames{corI})(1:end); dFF0_z_LH.(fieldNames{corI})(end-sampleN-1:end)];
%                     dFF0_z_RH.(fieldNames{corI}) = [dFF0_z_RH.(fieldNames{corI})(1:sampleN); dFF0_z_RH.(fieldNames{corI})(1:end); dFF0_z_RH.(fieldNames{corI})(end-sampleN-1:end)];

            end
            %% save the data
            FiberData.Ach.expCorrected = expCorrected_LH;
            FiberData.NE.expCorrected = expCorrected_RH;

            FiberData.Ach.lowPassData = lowPassData_LH;
            FiberData.NE.lowPassData = lowPassData_RH;

            FiberData.Ach.zScored= zScored_LH;
            FiberData.NE.zScored = zScored_RH;

            FiberData.Ach.motionCorrect= motionCorrect_LH;
            FiberData.NE.motionCorrect = motionCorrect_RH;

            FiberData.Ach.CBVCorrection = CBVCorrection_LH;
            FiberData.NE.CBVCorrection = CBVCorrection_RH;

            FiberData.Ach.dFF0_z = CBVCorrection_LH;
            FiberData.NE.dFF0_z = CBVCorrection_RH;

            FiberData.params = Params;

            clearvars -except FiberData fiberDataFilespath trialNum Fiberpath filepath
            
            save(Fiberpath,"FiberData",'-v7.3')

            clearvars -except fiberDataFilespath trialNum filepath
    end