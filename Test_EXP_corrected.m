close all;
% clear all; 
clc

    [z,p,k]=butter(3,0.001/(0.5*FiberData.params.DataFs),'high'); %Low pass filter for locomotion data
    [FiberData.params.sos_test_high,FiberData.params.g_test_high]=zp2sos(z,p,k);


            timeN = 15; % remove first 15 seconds;
            sampleN = floor(timeN*FiberData.params.DataFs);

        NRawData = FiberData.LH.rawData.F560;
RawData_LH =filtfilt(FiberData.params.sos_test_high,FiberData.params.g_test_high,NRawData(sampleN:end)); 
% RawData = NRawData(sampleN:end);
%% Filter Parameters
    [z,p,k]=butter(3,0.001/(0.5*FiberData.params.DataFs),'low'); %Low pass filter for locomotion data
    [FiberData.params.sos_test,FiberData.params.g_test]=zp2sos(z,p,k);
%% perform exponential correction
        FitData=filtfilt(FiberData.params.sos_Fit,FiberData.params.g_Fit,RawData_LH); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(FiberData.params.sos_Low,FiberData.params.g_Low,RawData_LH); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
% FiltData =   RawData;      
TESTFItData = filtfilt(FiberData.params.sos_test,FiberData.params.g_test,RawData_LH);
        %
        Spacing = 1:length(FiltData(:,1)); % sample index
        figTime=(1:length(FiltData))/(FiberData.params.DataFs*60);
    %% Correct TRITC blood volume Exp
        [fitVals]=fit(Spacing',FiltData,'exp2');%,'Exclude',FiberData.params.ExcludeVals);
        coeffVals=coeffvalues(fitVals);
        predicted.F560=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
    %% Correct TRITC blood volume Exp
        predicted.F560 = TESTFItData;
        % Plot the exponential fit
        figure; 
        h(1) = subplot(211);
        plot(figTime,detrend(RawData_LH,'constant')); hold on; plot(figTime,detrend(predicted.F560,'constant')); plot(figTime,detrend(FitData,'constant')); %plot(figTime,WheelData);
        title([ 'Detreding of TRITC metabolism ']); xlabel('Time (min)'); legend({'Raw TRITC brightness','Fit Line','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');

        expCorrected.F560_LH=FiltData-predicted.F560;
     
%             expCorrected.F560_LH=FiltData-predicted.F560';

        h(2)=subplot(212);
        plot(figTime,detrend(RawData_LH,'constant')); hold on; plot(figTime,detrend(expCorrected.F560_LH,'constant')); 
        title(['Detreding of TRITC metabolism ']); xlabel('Time (min)'); legend({'Raw TRITC brightness','expCorrected TRITC'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');
        linkaxes(h,'x');
 %%       
 NRawData = FiberData.RH.rawData.F560;
RawData_RH =filtfilt(FiberData.params.sos_test_high,FiberData.params.g_test_high,NRawData(sampleN:end)); 
% RawData = NRawData(sampleN:end);
%% Filter Parameters
    [z,p,k]=butter(3,0.001/(0.5*FiberData.params.DataFs),'low'); %Low pass filter for locomotion data
    [FiberData.params.sos_test,FiberData.params.g_test]=zp2sos(z,p,k);
%% perform exponential correction
        FitData=filtfilt(FiberData.params.sos_Fit,FiberData.params.g_Fit,RawData_RH); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(FiberData.params.sos_Low,FiberData.params.g_Low,RawData_RH); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
% FiltData =   RawData;   
TESTFItData = filtfilt(FiberData.params.sos_test,FiberData.params.g_test,RawData_RH);
        %
        Spacing = 1:length(FiltData(:,1)); % sample index
        figTime=(1:length(FiltData))/(FiberData.params.DataFs*60);
        %% Correct TRITC blood volume Exp
        [fitVals]=fit(Spacing',FiltData,'exp2');%,'Exclude',FiberData.params.ExcludeVals);
        coeffVals=coeffvalues(fitVals);

        predicted.F560=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
            %% Correct TRITC blood volume Exp
        predicted.F560 = TESTFItData;
        % Plot the exponential fit
        figure; 
        h(1) = subplot(211);
        plot(figTime,detrend(RawData_RH,'constant')); hold on; plot(figTime,detrend(predicted.F560,'constant')); plot(figTime,detrend(FitData,'constant')); %plot(figTime,WheelData);
        title([ 'Detreding of TRITC metabolism ']); xlabel('Time (min)'); legend({'Raw TRITC brightness','Fit Line','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');

        expCorrected.F560_RH=FiltData-predicted.F560;
     
%             expCorrected.F560_RH=FiltData-predicted.F560';

        h(2)=subplot(212);
        plot(figTime,detrend(RawData_RH,'constant')); hold on; plot(figTime,detrend(expCorrected.F560_RH,'constant')); 
        title(['Detreding of TRITC metabolism ']); xlabel('Time (min)'); legend({'Raw TRITC brightness','expCorrected TRITC'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');
        linkaxes(h,'x');


          figure;
       plot(figTime,detrend(expCorrected.F560_LH,'constant')); hold on
       plot(figTime,detrend(expCorrected.F560_RH,'constant')); 
       title(['expCorrected TRITC']); xlabel('Time (min)'); legend({'LH','RH'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');
        

           figure;
       plot(figTime,detrend(RawData_LH,'constant')); hold on
       plot(figTime,detrend(RawData_RH,'constant')); 
       title(['RawData']); xlabel('Time (min)'); legend({'LH','RH'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');