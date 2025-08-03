%% Set folder
dataLocation = [rootFolder delim 'Data']; % here is the key how to organize it
cd(dataLocation)

% File path
FileName = 'CompleteTable_AlcoholConsumptionData.xlsx';

% Read the data
data = readtable(FileName);

% Loop through each row to fill the structure
numAnimals = height(data);
for i = 1:numAnimals
    %AnimalID = matlab.lang.makeValidName(data.AnimalID{i});  % ensure valid fieldname
    AnimalID = data.AnimalID{i}; 
    AnimalSex = data.Sex{i}; 

    % Baseline data
    AlcoholDrinking.Baseline.(AnimalID).Sex = AnimalSex;
    AlcoholDrinking.Baseline.(AnimalID).Weight = data.Weight_Baseline(i);
    AlcoholDrinking.Baseline.(AnimalID).Day1 = data.Day1_Baseline(i);
    AlcoholDrinking.Baseline.(AnimalID).Day2 = data.Day2_Baseline(i);
    AlcoholDrinking.Baseline.(AnimalID).Day3 = data.Day3_Baseline(i);
    AlcoholDrinking.Baseline.(AnimalID).BingeDay = data.BingeDay_Baseline(i);

    % DID1 data
    AlcoholDrinking.DID1.(AnimalID).Sex = AnimalSex;
    AlcoholDrinking.DID1.(AnimalID).Weight = data.Weight_DID1(i);
    AlcoholDrinking.DID1.(AnimalID).Day1 = data.Day1_DID1(i);
    AlcoholDrinking.DID1.(AnimalID).Day2 = data.Day2_DID1(i);
    AlcoholDrinking.DID1.(AnimalID).Day3 = data.Day3_DID1(i);
    AlcoholDrinking.DID1.(AnimalID).BingeDay = data.BingeDay_DID1(i);

    % DID2 data
    AlcoholDrinking.DID2.(AnimalID).Sex = AnimalSex;
    AlcoholDrinking.DID2.(AnimalID).Weight = data.Weight_DID2(i);
    AlcoholDrinking.DID2.(AnimalID).Day1 = data.Day1_DID2(i);
    AlcoholDrinking.DID2.(AnimalID).Day2 = data.Day2_DID2(i);
    AlcoholDrinking.DID2.(AnimalID).Day3 = data.Day3_DID2(i);
    AlcoholDrinking.DID2.(AnimalID).BingeDay = data.BingeDay_DID2(i);

    % DID3 data
    AlcoholDrinking.DID3.(AnimalID).Sex = AnimalSex;
    AlcoholDrinking.DID3.(AnimalID).Weight = data.Weight_DID3(i);
    AlcoholDrinking.DID3.(AnimalID).Day1 = data.Day1_DID3(i);
    AlcoholDrinking.DID3.(AnimalID).Day2 = data.Day2_DID3(i);
    AlcoholDrinking.DID3.(AnimalID).Day3 = data.Day3_DID3(i);
    AlcoholDrinking.DID3.(AnimalID).BingeDay = data.BingeDay_DID3(i);

    % DID4 data
    AlcoholDrinking.DID4.(AnimalID).Sex = AnimalSex;
    AlcoholDrinking.DID4.(AnimalID).Weight = data.Weight_DID4(i);
    AlcoholDrinking.DID4.(AnimalID).Day1 = data.Day1_DID4(i);
    AlcoholDrinking.DID4.(AnimalID).Day2 = data.Day2_DID4(i);
    AlcoholDrinking.DID4.(AnimalID).Day3 = data.Day3_DID4(i);
    AlcoholDrinking.DID4.(AnimalID).BingeDay = data.BingeDay_DID4(i);
end

% Read TotalEtohConsumption sheet
totalEtohData = readtable(FileName, 'Sheet', 'TotalEtohConsumption');

% Loop to store total ethanol consumption per animal
for j = 1:height(totalEtohData)
    AnimalID = totalEtohData.AnimalID{j};
    AlcoholDrinking.TotalEtoh.(AnimalID).TotalEtohConsumed = totalEtohData.TotalEtohConsumed(j);
end


%% Save Data
cd([rootFolder delim 'Results_SST_project']);
save('Results_AlcoholDrinking.mat', 'AlcoholDrinking');
cd([rootFolder delim 'Data']);
