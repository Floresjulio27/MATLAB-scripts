%% Set folder
dataLocation = [rootFolder delim 'Data']; % here is the key how to organize it
cd(dataLocation)


% File path
FileName = 'CompleteTable_AlcoholConsumptionData.xlsx';

% Read the data starting from row 2 to skip sub-header
data = readtable(FileName);

% Preallocate structure array
numAnimals = height(data);
AlcoholDrinking = struct();

% Loop through each row to fill the structure
for i = 1:numAnimals
    AlcoholDrinking(i).AnimalID = string(data.AnimalID{i}); %string for not storing the data in cells 
    AlcoholDrinking(i).Sex = string(data.Sex{i});

    % Baseline data
    AlcoholDrinking(i).Baseline.Weight = data.Weight_Baseline(i);
    AlcoholDrinking(i).Baseline.Day1 = data.Day1_Baseline(i);
    AlcoholDrinking(i).Baseline.Day2 = data.Day2_Baseline(i);
    AlcoholDrinking(i).Baseline.Day3 = data.Day3_Baseline(i);
    AlcoholDrinking(i).Baseline.BingeDay = data.BingeDay_Baseline(i);

    % DID1 data
    AlcoholDrinking(i).DID1.Weight = data.Weight_DID1(i);
    AlcoholDrinking(i).DID1.Day1 = data.Day1_DID1(i);
    AlcoholDrinking(i).DID1.Day2 = data.Day2_DID1(i);
    AlcoholDrinking(i).DID1.Day3 = data.Day3_DID1(i);
    AlcoholDrinking(i).DID1.BingeDay = data.BingeDay_DID1(i);

    % DID2 data
    AlcoholDrinking(i).DID2.Weight = data.Weight_DID2(i);
    AlcoholDrinking(i).DID2.Day1 = data.Day1_DID2(i);
    AlcoholDrinking(i).DID2.Day2 = data.Day2_DID2(i);
    AlcoholDrinking(i).DID2.Day3 = data.Day3_DID2(i);
    AlcoholDrinking(i).DID2.BingeDay = data.BingeDay_DID2(i);

    % DID3 data
    AlcoholDrinking(i).DID3.Weight = data.Weight_DID3(i);
    AlcoholDrinking(i).DID3.Day1 = data.Day1_DID3(i);
    AlcoholDrinking(i).DID3.Day2 = data.Day2_DID3(i);
    AlcoholDrinking(i).DID3.Day3 = data.Day3_DID3(i);
    AlcoholDrinking(i).DID3.BingeDay = data.BingeDay_DID3(i);

    % DID4 data
    AlcoholDrinking(i).DID4.Weight = data.Weight_DID4(i);
    AlcoholDrinking(i).DID4.Day1 = data.Day1_DID4(i);
    AlcoholDrinking(i).DID4.Day2 = data.Day2_DID4(i);
    AlcoholDrinking(i).DID4.Day3 = data.Day3_DID4(i);
    AlcoholDrinking(i).DID4.BingeDay = data.BingeDay_DID4(i);
end

%% Save Data
cd([rootFolder delim 'Results_SST_project']);
save('Results_AlcoholDrinking.mat', 'AlcoholDrinking');
cd([rootFolder delim 'Data']);
