function TDTFiberPhotometry_bilateralGRABNE_Percent_modified(filepath,rawDataFilespath)

            if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
            end
        %% Check if the raw fiber data exist
        FiberDataFileStruct = dir('*_FiberData.mat');
        FiberDataFiles = {FiberDataFileStruct.name}';
        if size(FiberDataFiles,1) ~= 4 
            Data = TDTbin2mat(filepath);
            %%
            Params.DataFs=round(Data.streams.x560B.fs);
            Params.RawFs = Data.streams.Fitr.fs;
            Params.low_Freq = 1;
            Params.Fit_Freq=0.05;
            Params.plot_Freq=0.1;
            Params.ball_Freq=10;
            %% Extract the onsets 
            Params.DataOnset = Data.epocs.Top_.onset; % get the onset from the top epoc
            Params.DataOffset = Data.epocs.Bot_.offset;% get the offset from the bot epoc
            %% Filter Parameters
            [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
            [Params.sos_Low,Params.g_Low]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.ball_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for locomotion data
            [Params.sos_ball,Params.g_ball]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.Fit_Freq/(0.5*Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
            [Params.sos_Fit,Params.g_Fit]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.plot_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for plot
            [Params.sos_plot,Params.g_plot]=zp2sos(z,p,k);            
        end
%%
    for trialNum = 1:1:length(rawDataFilespath)
        %% generate file save path for each trials
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            Params.savepath =  [saveExt(1:strIdx-1) 'FiberData'];
            fiberFile = [Params.savepath '.mat'];
            %% Fiber signals
            if exist(fiberFile,"file")
                load(fiberFile);
                Raw_LH = FiberData.RawData.Raw_LH;
                Raw_RH = FiberData.RawData.Raw_RH;
                ForceSensor_raw = FiberData.RawData.ForceSensor_raw;
                Params = FiberData.params;
            else
                Raw_LH = [double(Data.streams.x405A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x465A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x560B.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
                Raw_RH = [double(Data.streams.x405C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x465C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x560D.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
                
                ForceSensor_raw = double(Data.streams.Forc.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))))';
            end
            %% remove data from the front to fix the exponential decay
            timeN = 15; % remove first 15 seconds;
            sampleStart = floor(timeN*Params.DataFs);
            RawData_LH = Raw_LH;
            RawData_RH = Raw_RH;
            ForceSensor = ForceSensor_raw;
            % adding dummy data to match index
            RawData_LH(1:sampleStart,:) = Raw_LH(sampleStart+1:sampleStart+sampleStart,:);
            RawData_RH(1:sampleStart,:) = Raw_RH(sampleStart+1:sampleStart+sampleStart,:);
            ForceSensor(1:sampleStart,:) = ForceSensor_raw(sampleStart+1:sampleStart+sampleStart,:);

            RawData_LH(sampleStart+1:end,:) = Raw_LH(sampleStart+1:end,:);
            RawData_RH(sampleStart+1:end,:) = Raw_RH(sampleStart+1:end,:);
            ForceSensor(sampleStart+1:end,:) = ForceSensor_raw(sampleStart+1:end,:);
            % add 15s back to first 
            ForceSensor_filt = filtfilt(Params.sos_ball,Params.g_ball,ForceSensor);
            %% detect any sudden spikes in the data
             time_N = (1:length(RawData_LH(:,2)))/Params.DataFs;
             % Outliers are defined as elements more than three 
             % local scaled MAD from the local median over a window length specified by window. This method is also known as a Hampel filter.
             [OR_RawData_LH,~] = filloutliers(RawData_LH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N); % slide using a 0.25 second window
             time_N = (1:length(RawData_RH(:,2)))/Params.DataFs;
             [OR_RawData_RH,~] = filloutliers(RawData_RH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N);
            %% remove outliers using median filter
            [DenoisedData_LH] =  Correct_motionartifact(OR_RawData_LH,Params);  
            [DenoisedData_RH] =  Correct_motionartifact(OR_RawData_RH,Params);
            %% Correct exponential decay
            time_Trim_Start = 210; %Data removed from the start
            time_Trim_End = 30; %Data removed from the end
            % [expCorrected_LH, ~] = Correct_ExponentialDecay_HighPass_Percent(OR_RawData_LH,ForceSensor,Params,'ACh',time_Trim,'y');
            % [expCorrected_RH, ~] = Correct_ExponentialDecay_HighPass_Percent(OR_RawData_RH,ForceSensor,Params,'NE',time_Trim,'y');
            % [expCorrected_LH, ~] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_LH,ForceSensor,Params,'ACh',time_Trim,'y');
            % [expCorrected_RH, ~] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_RH,ForceSensor,Params,'NE',time_Trim,'y');
            [expCorrectedData.Percentage.LH,expCorrectedData.ZScored.LH,Params] = Correct_Decay_Percentage_Change_HighPass_Lowpass(DenoisedData_LH,ForceSensor_filt,Params,'ACh',time_Trim_Start,time_Trim_End);
            [expCorrectedData.Percentage.RH,expCorrectedData.ZScored.RH,Params] = Correct_Decay_Percentage_Change_HighPass_Lowpass(DenoisedData_RH,ForceSensor_filt,Params,'NE',time_Trim_Start,time_Trim_End);
            %% Low pass the data  
            
            % % % if Data.info.Subject == "NEACh004"
            % % % expCorrectedData.ZScored.LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F405);
            % % % expCorrectedData.ZScored.LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F465);
            % % % expCorrectedData.ZScored.LH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.LH.F560);
            % % % expCorrectedData.ZScored.RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F405);
            % % % expCorrectedData.ZScored.RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F465);
            % % % expCorrectedData.ZScored.RH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.RH.F560);
            % % % elseif Data.info.Subject == "NEACh005"
            % % % expCorrectedData.ZScored.LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F405);
            % % % expCorrectedData.ZScored.LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F465);
            % % % expCorrectedData.ZScored.LH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.LH.F560);
            % % % expCorrectedData.ZScored.RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F405);
            % % % expCorrectedData.ZScored.RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F465);
            % % % expCorrectedData.ZScored.RH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.RH.F560);
            % % % elseif Data.info.Subject == "NEACh006"
            % % % expCorrectedData.ZScored.LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F405);
            % % % expCorrectedData.ZScored.LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F465);
            % % % expCorrectedData.ZScored.LH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.LH.F560);
            % % % expCorrectedData.ZScored.RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F405);
            % % % expCorrectedData.ZScored.RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F465);
            % % % expCorrectedData.ZScored.RH.F560=filtfilt(Params.sos_Low,Params.g_Low,-expCorrected.Percentage.RH.F560);
            % % % else
            % % % expCorrectedData.ZScored.LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F405);
            % % % expCorrectedData.ZScored.LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F465);
            % % % expCorrectedData.ZScored.LH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.LH.F560);
            % % % expCorrectedData.ZScored.RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F405);
            % % % expCorrectedData.ZScored.RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F465);
            % % % expCorrectedData.ZScored.RH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected.Percentage.RH.F560);
            % % % end
            %% hemodynamic Corrections
            CBVCorrectionSavePath = [Params.savepath 'CBVCorrectionParameters' '.mat'];
            [AnimalID,~,~] = GetFileInfo_FP(saveExt);
            AnimalIDStr = strfind(AnimalID,'NE');

            if AnimalID == "JF027" %%%|| AnimalID ==  "NEAChM003"
                if ~(exist(CBVCorrectionSavePath,'file') == 2)
                [CBVCorrection.Percent.LH.Slope,CBVCorrection.Percent.LH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.Percentage.LH,Params,'Percent','LH');
                [CBVCorrection.Percent.RH.Slope,CBVCorrection.Percent.RH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.Percentage.RH,Params,'Percent','RH');
                [CBVCorrection.ZScored.LH.Slope,CBVCorrection.ZScored.LH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.ZScored.LH,Params,'ZScored','LH');
                [CBVCorrection.ZScored.RH.Slope,CBVCorrection.ZScored.RH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.ZScored.RH,Params,'ZScored','RH');
                save(CBVCorrectionSavePath,"CBVCorrection",'-mat');
                else
                    load(CBVCorrectionSavePath); 
                end
            else % import the parameters used for correction
                if contains(AnimalID,'NE') % these parameters are calculated using GRAB mutants and should only be used for GRAB Sensors
                load('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Science_Research\SleepExperiment\CBV_Correction_Data\CBVMle_Hemodynamic_Correction_Parameters.mat')
                end
            end

            %% import the parameters used for correction
            % if exist("CBVCorrection","var")
            %     CorrectionParameters = CBVCorrection;
            % 
            % else
            % 
            % end
            %% perform the correction
            CBVCorrectionData.Percentage.LH = expCorrectedData.Percentage.LH;
            CBVCorrectionData.Percentage.RH = expCorrectedData.Percentage.RH;

            CBVCorrectionData.Percentage.LH.F465 = CBVCorrectionData.Percentage.LH.F465 - ((CBVCorrectionData.Percentage.LH.F560.*CBVCorrection.Percent.LH.Slope) + CBVCorrection.Percent.LH.Inter);
            CBVCorrectionData.Percentage.RH.F465 = CBVCorrectionData.Percentage.RH.F465 - ((CBVCorrectionData.Percentage.RH.F560.*CBVCorrection.Percent.RH.Slope) + CBVCorrection.Percent.RH.Inter);

            CBVCorrectionData.ZScored.LH = expCorrectedData.ZScored.LH;
            CBVCorrectionData.ZScored.RH = expCorrectedData.ZScored.RH;

            CBVCorrectionData.ZScored.LH.F465 = CBVCorrectionData.ZScored.LH.F465 - ((CBVCorrectionData.ZScored.LH.F560.*CBVCorrection.Percent.LH.Slope) + CBVCorrection.Percent.LH.Inter);
            CBVCorrectionData.ZScored.RH.F465 = CBVCorrectionData.ZScored.RH.F465 - ((CBVCorrectionData.ZScored.RH.F560.*CBVCorrection.Percent.RH.Slope) + CBVCorrection.Percent.RH.Inter);
            %plot the corrected signal
            plot_CBVCorrectedSignal(expCorrectedData.Percentage,CBVCorrectionData.Percentage,Params,'Percent')    
            plot_CBVCorrectedSignal(expCorrectedData.ZScored,CBVCorrectionData.ZScored,Params,'ZScored')
            %% house keeping
            sampleStart = floor(time_Trim_Start*Params.DataFs);
            sampleEnd= floor(time_Trim_End*Params.DataFs);
          
           % Match the data size
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    expCorrectedData.Percentage.LH.(fieldNames{corI}) = [expCorrectedData.Percentage.LH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.Percentage.LH.(fieldNames{corI})(1:end); expCorrectedData.Percentage.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    expCorrectedData.Percentage.RH.(fieldNames{corI}) = [expCorrectedData.Percentage.RH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.Percentage.RH.(fieldNames{corI})(1:end); expCorrectedData.Percentage.RH.(fieldNames{corI})(end-sampleEnd-2:end)];

                    expCorrectedData.ZScored.LH.(fieldNames{corI}) = [expCorrectedData.ZScored.LH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.ZScored.LH.(fieldNames{corI})(1:end); expCorrectedData.ZScored.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    expCorrectedData.ZScored.RH.(fieldNames{corI}) = [expCorrectedData.ZScored.RH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.ZScored.RH.(fieldNames{corI})(1:end); expCorrectedData.ZScored.RH.(fieldNames{corI})(end-sampleEnd-2:end)];
            end
            
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    CBVCorrectionData.ZScored.LH.(fieldNames{corI}) = [CBVCorrectionData.ZScored.LH.(fieldNames{corI})(1:sampleStart-1); CBVCorrectionData.ZScored.LH.(fieldNames{corI})(1:end); CBVCorrectionData.ZScored.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    CBVCorrectionData.ZScored.RH.(fieldNames{corI}) = [CBVCorrectionData.ZScored.RH.(fieldNames{corI})(1:sampleStart-1); CBVCorrectionData.ZScored.RH.(fieldNames{corI})(1:end); CBVCorrectionData.ZScored.RH.(fieldNames{corI})(end-sampleEnd-2:end)];

                    CBVCorrectionData.Percentage.LH.(fieldNames{corI}) = [CBVCorrectionData.Percentage.LH.(fieldNames{corI})(1:sampleStart); CBVCorrectionData.Percentage.LH.(fieldNames{corI})(1:end); CBVCorrectionData.Percentage.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    CBVCorrectionData.Percentage.RH.(fieldNames{corI}) = [CBVCorrectionData.Percentage.RH.(fieldNames{corI})(1:sampleStart); CBVCorrectionData.Percentage.RH.(fieldNames{corI})(1:end); CBVCorrectionData.Percentage.RH.(fieldNames{corI})(end-sampleEnd-2:end)];
            end
            %% save the data
            FiberData.RawData.Raw_LH = Raw_LH;
            FiberData.RawData.Raw_RH = Raw_RH;
            FiberData.RawData.ForceSensor_raw =  ForceSensor_raw;
                      
            FiberData.pressureSensor  = ForceSensor;
            FiberData.pressureSensor_filt  = ForceSensor_filt;

            FiberData.ACh.expCorrected.ZScored = expCorrectedData.ZScored.LH;
            FiberData.NE.expCorrected.ZScored = expCorrectedData.ZScored.RH;

            FiberData.ACh.dFF0_p = CBVCorrectionData.Percentage.LH;
            FiberData.NE.dFF0_p = CBVCorrectionData.Percentage.RH;

            FiberData.ACh.dFF0_z = CBVCorrectionData.ZScored.LH;
            FiberData.NE.dFF0_z = CBVCorrectionData.ZScored.RH;

            FiberData.params = Params;
            
            clearvars -except -regexp Data$ path$ Params$ trialNum$
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            savepathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            save(savepathmat,"FiberData",'-v7.3')
    end