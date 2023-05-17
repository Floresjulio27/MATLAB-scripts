% close all;
% clear all; 
% clc
function Plot_Raw_TDTFiberPhotometry_bilateralGCaMP7s(filepath,rawDataFilespath)
Data = TDTbin2mat(filepath);
%%
    Params.DataFs=round(Data.streams.x560B.fs);
    Params.RawFs = Data.streams.Fitr.fs;
    Params.low_Freq = 1;
    Params.Fit_Freq=0.05;
    Params.plot_Freq=0.1;
    %% Extract the onsets 
    Params.DataOnset = Data.epocs.Top_.onset; % get the onset from the top epoc
    Params.DataOffset = Data.epocs.Bot_.offset;% get the offset from the bot epoc
%% Filter Parameters
    [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
    [Params.sos_Low,Params.g_Low]=zp2sos(z,p,k);

    [z,p,k]=butter(3,10/(0.5*Params.DataFs),'low'); %Low pass filter for locomotion data
    [Params.sos_ball,Params.g_ball]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.Fit_Freq/(0.5*Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
    [Params.sos_Fit,Params.g_Fit]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.plot_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for locomotion data
    [Params.sos_plot,Params.g_plot]=zp2sos(z,p,k);

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
            %% Collected signal 
            Mod_Raw_LH_long = [double(Data.streams.Fitr.data(:,round(Params.DataOnset(trialNum)*Data.streams.Fitr.fs)+1:round(Data.streams.Fitr.fs*Params.DataOffset(trialNum))))]';
            Mod_Raw_RH_long = [double(Data.streams.Fipr.data(:,round(Params.DataOnset(trialNum)*Data.streams.Fitr.fs)+1:round(Data.streams.Fitr.fs*Params.DataOffset(trialNum))))]';
            resample_freq = 100;
            Mod_Raw_RH = detrend(resample(Mod_Raw_RH_long(round(180*Data.streams.Fitr.fs):end-round(180*Data.streams.Fitr.fs),:),resample_freq, round(Data.streams.Fitr.fs)),'constant');
            Mod_Raw_LH = detrend(resample(Mod_Raw_LH_long(round(180*Data.streams.Fitr.fs):end-round(180*Data.streams.Fitr.fs),:),resample_freq, round(Data.streams.Fitr.fs)),'constant');
            ForceSensor_raw_N = detrend(resample(ForceSensor_raw(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');

            MplotTime = (1:length(Mod_Raw_LH(:,1)))/(resample_freq*60);
            NPlotTime = (1:length(ForceSensor_raw_N))/(resample_freq*60);
            figure;
            h(1)=subplot(511);
            plot(NPlotTime,ForceSensor_raw_N)
            legend('Force Sensor')
            h(2)=subplot(512);
            plot(MplotTime,Mod_Raw_LH(:,1))
            legend('LH Raw 465')
            title('Raw signals recorded from the fibers')
            h(3)= subplot(513);
            plot(MplotTime,Mod_Raw_LH(:,2))
            legend('LH Raw 560')
            h(4)=subplot(514);
            plot(MplotTime,Mod_Raw_RH(:,1))
            legend('RH Raw 465')
            h(5)=subplot(515);
            plot(MplotTime,Mod_Raw_RH(:,2))
            legend('RH Raw 560')
            ylabel('Raw signals')
            xlabel('Time (min)')
            linkaxes(h);
             xlim([1 54])
%             xlim([1 60])
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawSignal.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawSignal.tiff'],'tiff')
            close 
            
            figure;

            resample_freq = 100;
            Raw_LH_N = detrend(resample(Raw_LH(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');
            Raw_RH_N = detrend(resample(Raw_RH(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');
            NPlotTime = (1:length(ForceSensor_raw_N))/(resample_freq*60);
            h(1)=subplot(711);
            plot(NPlotTime,ForceSensor_raw_N)
            legend('Force Sensor')
            title('Raw signals recorded from the fibers')
            h(2)=subplot(712);
            plot(NPlotTime,Raw_LH_N(:,1))
            legend('LH Raw 405')
            h(3)=subplot(713);
            plot(NPlotTime,Raw_LH_N(:,2))
            legend('LH Raw 465')
            h(4)= subplot(714);
            plot(NPlotTime,Raw_LH_N(:,3))
            legend('LH Raw 560')

            h(5)=subplot(715);
            plot(NPlotTime,Raw_RH_N(:,1))
            legend('RH Raw 405')
            h(6)=subplot(716);
            plot(NPlotTime,Raw_RH_N(:,2))
            legend('RH Raw 405')
            h(7)=subplot(717);
            plot(NPlotTime,Raw_RH_N(:,3))
            legend('RH Raw 560')
            ylabel('Raw signals')
            xlabel('Time (min)')
            linkaxes(h,'x');
            xlim([1 55])
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawData.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawData.tiff'],'tiff')
            close 
            %% save the data
            loadpathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            load(loadpathmat);

            FiberData.LH.ModrawData = Mod_Raw_LH_long;
            FiberData.RH.ModrawData = Mod_Raw_RH_long;

          
 
            clearvars -except -regexp Data$ path$ trialNum$ Params$
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            savepathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            save(savepathmat,"FiberData",'-v7.3')
    end
%% Find movement and separate those from run
%     speedaccBinary = speedProcess_Doric(AnalogWheel,Params.DataFs);
%     
%     Speed_Out = findSpeedSeg(speedaccBinary, Params.DataFs,3, 5);
%     
%     [events_ch1]=eventevokedchanges(dF_ch1,Speed_Out);
%     [events_ch2]=eventevokedchanges(dF_ch2,Speed_Out);
    %%
%    save(['N:\TDTData\GFP-021-220719-152016\FiberprocessedData.mat'],'-v7.3')