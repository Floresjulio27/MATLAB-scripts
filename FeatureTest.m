function FeatureTest
%% load in all the data to create a table of values
startingDirectory = cd;
trainingDirectory = [startingDirectory];% '\Training Data\'];
cd(trainingDirectory)
% character list of all training files
trainingDataFileStruct = dir('*_ModelData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
% Load each updated training set and concatenate the data into table
for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    dataLength = size(ModelData.TrainLabels,1);
    for lp = 1:1:dataLength
        AllFeatures{((bb-1)*dataLength)+lp } = ModelData.Features{lp}';
    end
    
    AllLabels(((bb-1)*dataLength)+1 : (bb*dataLength)) = ModelData.TrainLabels;
end
for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    dataLength = size(ModelData.TrainLabels,1);
    for lp = 1:1:dataLength
       Test_EMG(:,((bb-1)*dataLength)+lp) = ModelData.Features{lp}(:,12);
    end
    
end
for jl=1:1:length(AllLabels)
    if AllLabels(jl) == "Not Sleep"
    AllLabelsN(jl) = 1;
    elseif AllLabels(jl) == "NREM Sleep"
        AllLabelsN(jl) = 2;
    elseif AllLabels(jl) == "REM Sleep"
        AllLabelsN(jl) = 3;
    end
end
g_table = table(AllLabelsN',Test_EMG');
fitglme(g_table,'Var2~1+Var1')
end