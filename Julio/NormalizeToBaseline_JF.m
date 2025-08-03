clc, clear; 
load ("AnalysisResults_DID.mat")

%% Define variables 

mice = {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF050'};
%week = {'Baseline', 'DID1', 'DID2', 'DID3', 'DID4', 'Baseline2'};
week = {'Baseline', 'DID1', 'DID2'};
regions = {'PFC', 'S1BF'};
transitions = {'AWAKEtoNREM', 'NREMtoAWAKE', 'NREMtoREM', 'REMtoAWAKE'};

%% Organize data

% Initialize structure to store data
data = struct();

for iRegion = 1:length(regions)
    region = regions{iRegion};

    for iWeek = 1:length(week)
        weeks = week{iWeek};

        for iTransition = 1:length(transitions)
            transition = transitions{iTransition};

            % Initialize an empty matrix for this transition if not already created
            %if ~isfield(data, weeks) || ~isfield(data.(weeks), region) || ~isfield(data.(weeks).(region), transition)
                data.(weeks).(region).(transition) = [];
            %end

            for iMouse = 1:length(mice)
                mouse = mice{iMouse};
                try
                    % Append the value for the current mouse to the next column
                    value = AnalysisResults_DID.(mouse).(weeks).(region).(transition);
                    data.(weeks).(region).(transition) = [data.(weeks).(region).(transition), value]; %concatenating values 
                catch
                    % If data is missing, skip this combination.
                    fprintf('Data missing for %s, %s, %s, %s\n', mouse, region, weeks, transition);
                    data.(weeks).(region).(transition) = [data.(weeks).(region).(transition), NaN]; % Insert NaN for missing data
                end
            end
        end
    end 
end

%% Normalize data

%create structure 
NormalizedData = struct ();

%PFC 
NormAWAKEtoNREMbaseline = data.Baseline.PFC.AWAKEtoNREM ./ data.Baseline.PFC.AWAKEtoNREM;
NormAWAKEtoNREMDID1 = data.Baseline.PFC.AWAKEtoNREM ./ data.DID1.PFC.AWAKEtoNREM;
NormAWAKEtoNREMDID2 = data.Baseline.PFC.AWAKEtoNREM ./ data.DID2.PFC.AWAKEtoNREM;

NormNREMtoAWAKEbaseline = data.Baseline.PFC.NREMtoAWAKE ./ data.Baseline.PFC.NREMtoAWAKE;
NormNREMtoAWAKEDID1 = data.Baseline.PFC.NREMtoAWAKE ./ data.DID1.PFC.NREMtoAWAKE;
NormNREMtoAWAKEDID2 = data.Baseline.PFC.NREMtoAWAKE ./ data.DID2.PFC.NREMtoAWAKE;

NormNREMtoREMbaseline = data.Baseline.PFC.NREMtoREM ./ data.Baseline.PFC.NREMtoREM;
NormNREMtoREMDID1 = data.Baseline.PFC.NREMtoREM ./ data.DID1.PFC.NREMtoREM;
NormNREMtoREMDID2 = data.Baseline.PFC.NREMtoREM ./ data.DID2.PFC.NREMtoREM;

NormREMtoAWAKEbaseline = data.Baseline.PFC.REMtoAWAKE ./ data.Baseline.PFC.REMtoAWAKE;
NormREMtoAWAKEDID1 = data.Baseline.PFC.REMtoAWAKE ./ data.DID1.PFC.REMtoAWAKE;
NormREMtoAWAKEDID2 = data.Baseline.PFC.REMtoAWAKE ./ data.DID2.PFC.REMtoAWAKE;

%S1BF
NormAWAKEtoNREMbaseline = data.Baseline.S1BF.AWAKEtoNREM ./ data.Baseline.S1BF.AWAKEtoNREM;
NormAWAKEtoNREMDID1 = data.Baseline.S1BF.AWAKEtoNREM ./ data.DID1.S1BF.AWAKEtoNREM;
NormAWAKEtoNREMDID2 = data.Baseline.S1BF.AWAKEtoNREM ./ data.DID2.S1BF.AWAKEtoNREM;

NormNREMtoAWAKEbaseline = data.Baseline.S1BF.NREMtoAWAKE ./ data.Baseline.S1BF.NREMtoAWAKE;
NormNREMtoAWAKEDID1 = data.Baseline.S1BF.NREMtoAWAKE ./ data.DID1.S1BF.NREMtoAWAKE;
NormNREMtoAWAKEDID2 = data.Baseline.S1BF.NREMtoAWAKE ./ data.DID2.S1BF.NREMtoAWAKE;

NormNREMtoREMbaseline = data.Baseline.S1BF.NREMtoREM ./ data.Baseline.S1BF.NREMtoREM;
NormNREMtoREMDID1 = data.Baseline.S1BF.NREMtoREM ./ data.DID1.S1BF.NREMtoREM;
NormNREMtoREMDID2 = data.Baseline.S1BF.NREMtoREM ./ data.DID2.S1BF.NREMtoREM;

NormREMtoAWAKEbaseline = data.Baseline.S1BF.REMtoAWAKE ./ data.Baseline.S1BF.REMtoAWAKE;
NormREMtoAWAKEDID1 = data.Baseline.S1BF.REMtoAWAKE ./ data.DID1.S1BF.REMtoAWAKE;
NormREMtoAWAKEDID2 = data.Baseline.S1BF.REMtoAWAKE ./ data.DID2.S1BF.REMtoAWAKE;

%Save the results 
NormalizedData.Baseline.Transitions.AWAKEtoNREM = NormAWAKEtoNREMbaseline;
NormalizedData.Baseline.Transitions.NREMtoAWAKE = NormNREMtoAWAKEbaseline;
NormalizedData.Baseline.Transitions.NREMtoREM = NormNREMtoREMbaseline;
NormalizedData.Baseline.Transitions.REMtoAWAKE = NormREMtoAWAKEbaseline;

NormalizedData.DID1.Transitions.AWAKEtoNREM = NormAWAKEtoNREMDID1;
NormalizedData.DID1.Transitions.NREMtoAWAKE = NormNREMtoAWAKEDID1;
NormalizedData.DID1.Transitions.NREMtoREM = NormNREMtoREMDID1;
NormalizedData.DID1.Transitions.REMtoAWAKE = NormREMtoAWAKEDID1;

NormalizedData.DID2.Transitions.AWAKEtoNREM = NormAWAKEtoNREMDID2;
NormalizedData.DID2.Transitions.NREMtoAWAKE = NormNREMtoAWAKEDID2;
NormalizedData.DID2.Transitions.NREMtoREM = NormNREMtoREMDID2;
NormalizedData.DID2.Transitions.REMtoAWAKE = NormREMtoAWAKEDID2;

%% Plot the data 




                
