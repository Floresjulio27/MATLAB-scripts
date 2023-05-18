% close all;
% clear all; 
% clc
function TDTFiberPhotometry_bilateralGCaMP7s_CheckFilter(filepath,rawDataFilespath)
Data = TDTbin2mat(filepath);
%%
    Params.DataFs=round(Data.streams.x560B.fs);
    Params.RawFs = Data.streams.Fitr.fs;
    Params.low_Freq = 1;
    Params.Fit_Freq=0.05;
    Params.plot_Freq=0.1;
    Params.ball_Freq=10;
    Params.decay_Freq=1/720; % 0.0014
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

    [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
    [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);

    for trialNum = 1:1:length(Params.DataOnset)
        %% generate file save path for each trials
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            Params.savepath =  [saveExt(1:strIdx-1) 'FiberData'];
        %% Fiber signals
            Raw_LH = [double(Data.streams.x405A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x465A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x560B.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
            Raw_RH = [double(Data.streams.x405C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x465C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x560D.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
            
            ForceSensor_raw = double(Data.streams.Forc.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))))';            
          
            %% remove data from the front to fix the exponential decay
            timeN = 1; % remove first 15 seconds;
            sampleN = floor(timeN*Params.DataFs);
            RawData_LH = Raw_LH;
            RawData_RH = Raw_RH;
            ForceSensor = ForceSensor_raw;
            % adding dummy data to match index
            RawData_LH(1:sampleN,:) = Raw_LH(sampleN+1:sampleN+sampleN,:);
            RawData_RH(1:sampleN,:) = Raw_RH(sampleN+1:sampleN+sampleN,:);
            ForceSensor(1:sampleN,:) = ForceSensor_raw(sampleN+1:sampleN+sampleN,:);

            RawData_LH(sampleN+1:end,:) = Raw_LH(sampleN+1:end,:);
            RawData_RH(sampleN+1:end,:) = Raw_RH(sampleN+1:end,:);
            ForceSensor(sampleN+1:end,:) = ForceSensor_raw(sampleN+1:end,:);
            % add 15s back to first 
            ForceSensor_filt = filtfilt(Params.sos_ball,Params.g_ball,ForceSensor);

            %% detect any sudden spikes in the data
             Det_RawData_RH = RawData_RH;
             Det_RawData_LH = RawData_LH;
             time_N = (1:length(Det_RawData_LH(:,2)))/Params.DataFs;

             %	Outliers are defined as elements more than three 
             % local scaled MAD from the local median over a window length specified by window. This method is also known as a Hampel filter.
             [C_RawData_LH,TF_Outliers_RawData_LH] = filloutliers(Det_RawData_LH,"linear","movmedian",round(0.25*Params.DataFs),"SamplePoints",time_N);

             figure;
             plot(time_N,Det_RawData_LH(:,2));
             hold on
             plot(time_N,C_RawData_LH(:,2),'--')
             hold off

             figure;
             plot(time_N,Det_RawData_LH(:,1));
             hold on
             plot(time_N,C_RawData_LH(:,1),'--')
             hold off

              figure;
             plot(time_N,Det_RawData_LH(:,3));
             hold on
             plot(time_N,C_RawData_LH(:,3),'--')
             hold off

             time_N = (1:length(Det_RawData_RH(:,2)))/Params.DataFs;

             [C_RawData_RH,TF_Outliers_RawData_RH]  = filloutliers(Det_RawData_RH,"linear","movmedian",round(0.25*Params.DataFs),"SamplePoints",time_N);

             figure;
             plot(time_N,Det_RawData_RH(:,2));
             hold on
             plot(time_N,C_RawData_RH(:,2),'--')
             hold off

             figure;
             plot(time_N,Det_RawData_RH(:,1));
             hold on
             plot(time_N,C_RawData_RH(:,1),'--')
             hold off

             figure;
             plot(time_N,Det_RawData_RH(:,3));
             hold on
             plot(time_N,C_RawData_RH(:,3),'--')
             hold off

             close all
    end