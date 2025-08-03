function [] = AddSleepParameters_FP_GRABNE_SingleFiber(procDataFileIDs,RestingBaselines,baselineType,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Organize data into appropriate bins for sleep scoring characterization
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']);
    [~, fileDate, ~] = GetFileInfo_FP(procDataFileID);
    strDay = ConvertDate_FP(fileDate);
    load(procDataFileID)
    specDataFileID = [procDataFileID(1:end-12) 'SpecDataA.mat'];
    load(specDataFileID)

    %% Neural Data Binning
    LH_Delta = ProcData.data.cortical_LH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.(baselineType).cortical_LH.deltaBandPower.(strDay).mean;
    LH_DeltaNeuro = (LH_Delta - LH_baselineDelta) / LH_baselineDelta;
    
    LH_Theta = ProcData.data.cortical_LH.thetaBandPower;
    LH_baselineTheta = RestingBaselines.(baselineType).cortical_LH.thetaBandPower.(strDay).mean;
    LH_ThetaNeuro = (LH_Theta - LH_baselineTheta) / LH_baselineTheta;
    
    LH_Alpha = ProcData.data.cortical_LH.alphaBandPower;
    LH_baselineAlpha = RestingBaselines.(baselineType).cortical_LH.alphaBandPower.(strDay).mean;
    LH_AlphaNeuro = (LH_Alpha - LH_baselineAlpha) / LH_baselineAlpha;
    
    LH_Beta = ProcData.data.cortical_LH.betaBandPower;
    LH_baselineBeta = RestingBaselines.(baselineType).cortical_LH.betaBandPower.(strDay).mean;
    LH_BetaNeuro = (LH_Beta - LH_baselineBeta) / LH_baselineBeta;
    
    LH_Gamma = ProcData.data.cortical_LH.gammaBandPower;
    LH_baselineGamma = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay).mean;
    LH_GammaNeuro = (LH_Gamma - LH_baselineGamma) / LH_baselineGamma;
    
    LH_Cortical = ProcData.data.cortical_LH.corticalPower;
    LH_baselineCortical = RestingBaselines.(baselineType).cortical_LH.corticalPower.(strDay).mean;
    LH_CorticalNeuro = (LH_Cortical - LH_baselineCortical) / LH_baselineCortical;

    LH_tempDeltaStruct = cell(NBins,1);
    LH_tempThetaStruct = cell(NBins,1);
    LH_tempAlphaStruct = cell(NBins,1);
    LH_tempBetaStruct = cell(NBins,1);
    LH_tempGammaStruct = cell(NBins,1);
    LH_tempCorticalStruct = cell(NBins,1);
    
    % Binning neural signals into 5-second bins
    for b = 1:NBins
        startIndex = (150 * (b - 1)) + 1;
        endIndex = min(150 * b, length(LH_DeltaNeuro));
        
        LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro(startIndex:endIndex)};
        LH_tempThetaStruct(b,1) = {LH_ThetaNeuro(startIndex:endIndex)};
        LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro(startIndex:endIndex)};
        LH_tempBetaStruct(b,1) = {LH_BetaNeuro(startIndex:endIndex)};
        LH_tempGammaStruct(b,1) = {LH_GammaNeuro(startIndex:endIndex)};
        LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro(startIndex:endIndex)};
    end

    % Save neural data under ProcData file
    ProcData.sleep.parameters.cortical_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_LH.thetaBandPower = LH_tempThetaStruct;
    ProcData.sleep.parameters.cortical_LH.alphaBandPower = LH_tempAlphaStruct;
    ProcData.sleep.parameters.cortical_LH.betaBandPower = LH_tempBetaStruct;
    ProcData.sleep.parameters.cortical_LH.gammaBandPower = LH_tempGammaStruct;
    ProcData.sleep.parameters.cortical_LH.corticalPower = LH_tempCorticalStruct;
    
    %% Spectrogram Data Binning
    trialDuration_sec = ProcData.notes.trialDuration_sec;   % sec
    offset = 2.5;   % sec
    binWidth = 5;   % sec
    T = round(SpecData.cortical_LH.T,1);
    F = SpecData.cortical_LH.F;
    specLH = SpecData.cortical_LH.normS;
    freqFloor = floor(F);
    
    % Delta spectrogram band
    deltaLow = freqFloor == 1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow,1,'first');
    deltaLowEnd = find(deltaHigh,1,'last');
    deltaSpecLH = specLH(deltaLowStart:deltaLowEnd,:);
    meanDeltaSpecLH = mean(deltaSpecLH,1);
    
    % Theta spectrogram band
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 9;
    thetaLowStart = find(thetaLow,1,'first');
    thetaLowEnd = find(thetaHigh,1,'last');
    thetaSpecLH = specLH(thetaLowStart:thetaLowEnd,:);
    meanThetaSpecLH = mean(thetaSpecLH,1);
    
    % Alpha spectrogram band
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow,1,'first');
    alphaLowEnd = find(alphaHigh,1,'last');
    alphaSpecLH = specLH(alphaLowStart:alphaLowEnd,:);
    meanAlphaSpecLH = mean(alphaSpecLH,1);
    
    % Beta spectrogram band
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow,1,'first');
    betaLowEnd = find(betaHigh,1,'last');
    betaSpecLH = specLH(betaLowStart:betaLowEnd,:);
    meanBetaSpecLH = mean(betaSpecLH,1);
    
    % Gamma spectrogram band
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow,1,'first');
    gammaLowEnd = find(gammaHigh,1,'last');
    gammaSpecLH = specLH(gammaLowStart:gammaLowEnd,:);
    meanGammaSpecLH = mean(gammaSpecLH,1);
    
    % Binning spectrogram data
    LH_tempDeltaSpecStruct = cell(NBins,1);
    LH_tempThetaSpecStruct = cell(NBins,1);
    LH_tempAlphaSpecStruct = cell(NBins,1);
    LH_tempBetaSpecStruct = cell(NBins,1);
    LH_tempGammaSpecStruct = cell(NBins,1);
    
    for c = 1:NBins
        startTime = binWidth * (c - 1);
        endTime = binWidth * c;
        [~, startTime_index] = min(abs(T - startTime));
        [~, endTime_index] = min(abs(T - endTime));
        
        startTime_index = max(1, startTime_index);
        endTime_index = min(endTime_index, length(meanDeltaSpecLH));
        
        LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
        LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
        LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
        LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
        LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
    end
    
    % Save spectrogram data
    ProcData.sleep.parameters.cortical_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specThetaBandPower = LH_tempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specAlphaBandPower = LH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specGammaBandPower = LH_tempGammaSpecStruct;
  %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
binWhiskerAngle = ProcData.data.binWhiskerAngle;
binForceSensor = ProcData.data.binForceSensor;
ForceSensor = ProcData.data.forceSensor;
whiskerAngle = ProcData.data.whiskerAngle;
whiskerAcceleration = diff(whiskerAngle, 2);  % Acceleration from angle differences

% Find the number of whisker bins based on the NBins
whiskerBinNumber = NBins;

% Divide the signals into five-second bins and store them in cell arrays
tempWhiskerStruct = cell(whiskerBinNumber, 1);
tempWhiskerAccelStruct = cell(whiskerBinNumber, 1);
tempBinWhiskerStruct = cell(whiskerBinNumber, 1);
tempForceStruct = cell(whiskerBinNumber, 1);
tempForceRStruct = cell(whiskerBinNumber, 1);

for whiskerBins = 1:whiskerBinNumber
    startIndex = (150 * (whiskerBins - 1)) + 1;
    endIndex = min(150 * whiskerBins, length(whiskerAngle));  % Ensure bounds are respected

    % Ensure that the indices do not exceed the length of each array
    tempWhiskerStruct{whiskerBins, 1} = whiskerAngle(startIndex:endIndex);
    tempWhiskerAccelStruct{whiskerBins, 1} = whiskerAcceleration(startIndex:endIndex);
    tempBinWhiskerStruct{whiskerBins, 1} = binWhiskerAngle(startIndex:endIndex);
    tempForceStruct{whiskerBins, 1} = binForceSensor(startIndex:endIndex);
    tempForceRStruct{whiskerBins, 1} = ForceSensor(startIndex:endIndex);
end

% Store whisker and force sensor data in the ProcData file
ProcData.sleep.parameters.whiskerAngle = tempWhiskerStruct;
ProcData.sleep.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
ProcData.sleep.parameters.binWhiskerAngle = tempBinWhiskerStruct;
ProcData.sleep.parameters.binForceSensor = tempForceStruct;
ProcData.sleep.parameters.ForceSensor = tempForceRStruct;


%% Add Pupil Parameters
ProcData.data.Pupil.mmPerPixel = 0.018;  % Conversion from pixels to mm

% Create fields for Pupil data binning
dataTypes = {'zDiameter'};  % Other data types can be added as needed
samplingRate = ProcData.notes.dsFs;
[z, p, k] = butter(4, 1 / (samplingRate / 2), 'low');
[sos, g] = zp2sos(z, p, k);

for aa = 1:length(dataTypes)
    dataType = dataTypes{aa};
    pupilData = ProcData.data.Pupil.(dataType);
    
    % Filter the pupil data
    try
        filteredData = filtfilt(sos, g, pupilData);
    catch
        filteredData = pupilData;  % If filtering fails, use raw data
    end
    
    tempPupilStruct = cell(NBins, 1);
    for pupilBins = 1:NBins
        startIndex = (150 * (pupilBins - 1)) + 1;
        endIndex = min(150 * pupilBins, length(filteredData));
        tempPupilStruct{pupilBins, 1} = filteredData(startIndex:endIndex);
    end
    
    % Save pupil data
    ProcData.sleep.parameters.Pupil.(dataType) = tempPupilStruct;
end

%% Create folder for the EMG
EMG = ProcData.data.EMG.emg;
normEMG = EMG - RestingBaselines.(baselineType).EMG.emg.(strDay).mean;
tempEMGStruct = cell(NBins, 1);

for EMGBins = 1:NBins
    startIndex = (150 * (EMGBins - 1)) + 1;
    endIndex = min(150 * EMGBins, length(normEMG));
    tempEMGStruct{EMGBins, 1} = normEMG(startIndex:endIndex);
end

ProcData.sleep.parameters.EMG.emg = tempEMGStruct;

%% Raw EMG Signal
rawEMGSignal = ProcData.data.EMG.emgSignal;
normEMGSignal = rawEMGSignal - RestingBaselines.(baselineType).EMG.emgSignal.(strDay).mean;
tempRawEMGStruct = cell(NBins, 1);

for EMGBins = 1:NBins
    startIndex = (150 * (EMGBins - 1)) + 1;
    endIndex = min(150 * EMGBins, length(normEMGSignal));
    tempRawEMGStruct{EMGBins, 1} = normEMGSignal(startIndex:endIndex);
end

ProcData.sleep.parameters.EMG.emgSignal = tempRawEMGStruct;

%% Create folder for right hemisphere fiber data (z-scored)
Z_NE_GFP = ProcData.data.GFP.Z_NE;
Z_NE_CBV = ProcData.data.CBV.Z_NE;

Z_NE_NormCBV = (Z_NE_CBV - RestingBaselines.(baselineType).CBV.Z_NE.(strDay).mean);
Z_NE_NormGFP = (Z_NE_GFP - RestingBaselines.(baselineType).GFP.Z_NE.(strDay).mean);

Z_NE_tempGFPStruct = cell(NBins, 1);
CBVZ_NE_tempCBVStruct = cell(NBins, 1);

for CBVBins = 1:NBins
    startIndex = (150 * (CBVBins - 1)) + 1;
    endIndex = min(150 * CBVBins, length(Z_NE_NormGFP));
    Z_NE_tempGFPStruct{CBVBins, 1} = Z_NE_NormGFP(startIndex:endIndex);
    CBVZ_NE_tempCBVStruct{CBVBins, 1} = Z_NE_NormCBV(startIndex:endIndex);
end

% Save fiber data
ProcData.sleep.parameters.GFP.Z_NE = Z_NE_tempGFPStruct;
ProcData.sleep.parameters.CBV.Z_NE = CBVZ_NE_tempCBVStruct;

%% Create folder for the right hemisphere fiber data (Percentage)
P_NE_GFP = ProcData.data.GFP.P_NE;
P_NE_CBV = ProcData.data.CBV.P_NE;

P_NE_NormCBV = (P_NE_CBV - RestingBaselines.(baselineType).CBV.P_NE.(strDay).mean);
P_NE_NormGFP = (P_NE_GFP - RestingBaselines.(baselineType).GFP.P_NE.(strDay).mean);

P_NE_tempGFPStruct = cell(NBins, 1);
CBVP_NE_tempCBVStruct = cell(NBins, 1);

for CBVBins = 1:NBins
    startIndex = (150 * (CBVBins - 1)) + 1;
    endIndex = min(150 * CBVBins, length(P_NE_NormGFP));
    P_NE_tempGFPStruct{CBVBins, 1} = P_NE_NormGFP(startIndex:endIndex);
    CBVP_NE_tempCBVStruct{CBVBins, 1} = P_NE_NormCBV(startIndex:endIndex);
end

% Save fiber data
ProcData.sleep.parameters.GFP.P_NE = P_NE_tempGFPStruct;
ProcData.sleep.parameters.CBV.P_NE = CBVP_NE_tempCBVStruct;
%% BLOCK PURPOSE: Create folder for the left and right fiber data (z-scored)
if isfield(ProcData.data.GFP, 'Z_ACh')
    Z_ACh_GFP = ProcData.data.GFP.Z_ACh; 
    Z_ACh_CBV = ProcData.data.CBV.Z_ACh;

    % Normalize the CBV and GFP data
    Z_ACh_NormCBV = (Z_ACh_CBV - RestingBaselines.(baselineType).CBV.Z_ACh.(strDay).mean) / RestingBaselines.(baselineType).CBV.Z_ACh.(strDay).std;
    Z_ACh_NormGFP = (Z_ACh_GFP - RestingBaselines.(baselineType).GFP.Z_ACh.(strDay).mean) / RestingBaselines.(baselineType).GFP.Z_ACh.(strDay).std;

    Z_ACh_tempGFPStruct = cell(NBins, 1);
    CBVZ_ACh_tempCBVStruct = cell(NBins, 1);

    for CBVBins = 1:NBins
        startIndex = (150 * (CBVBins - 1)) + 1;
        endIndex = min(150 * CBVBins, length(Z_ACh_NormGFP));  % Ensure bounds are respected

        Z_ACh_tempGFPStruct{CBVBins, 1} = Z_ACh_NormGFP(startIndex:endIndex);
        CBVZ_ACh_tempCBVStruct{CBVBins, 1} = Z_ACh_NormCBV(startIndex:endIndex);
    end

    % Save the normalized fiber data under ProcData file
    ProcData.sleep.parameters.GFP.Z_ACh = Z_ACh_tempGFPStruct;
    ProcData.sleep.parameters.CBV.Z_ACh = CBVZ_ACh_tempCBVStruct;
end

%% BLOCK PURPOSE: Create folder for the left hemisphere fiber data (z-scored)
if isfield(ProcData.data.GFP, 'P_ACh')
    P_ACh_GFP = ProcData.data.GFP.P_ACh; 
    P_ACh_CBV = ProcData.data.CBV.P_ACh;

    % Normalize the CBV and GFP data
    P_ACh_NormCBV = (P_ACh_CBV - RestingBaselines.(baselineType).CBV.P_ACh.(strDay).mean) / RestingBaselines.(baselineType).CBV.P_ACh.(strDay).std;
    P_ACh_NormGFP = (P_ACh_GFP - RestingBaselines.(baselineType).GFP.P_ACh.(strDay).mean) / RestingBaselines.(baselineType).GFP.P_ACh.(strDay).std;

    P_ACh_tempGFPStruct = cell(NBins, 1);
    CBVP_ACh_tempCBVStruct = cell(NBins, 1);

    for CBVBins = 1:NBins
        startIndex = (150 * (CBVBins - 1)) + 1;
        endIndex = min(150 * CBVBins, length(P_ACh_NormGFP));  % Ensure bounds are respected

        P_ACh_tempGFPStruct{CBVBins, 1} = P_ACh_NormGFP(startIndex:endIndex);
        CBVP_ACh_tempCBVStruct{CBVBins, 1} = P_ACh_NormCBV(startIndex:endIndex);
    end

    % Save the normalized fiber data under ProcData file
    ProcData.sleep.parameters.GFP.P_ACh = P_ACh_tempGFPStruct;
    ProcData.sleep.parameters.CBV.P_ACh = CBVP_ACh_tempCBVStruct;
end

% Save the updated ProcData file
save(procDataFileID, 'ProcData');
end 