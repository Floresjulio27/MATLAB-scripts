function [] = UpdateTrainingDataSets_FP_JFTrim(procDataFileIDs)
    %________________________________________________________________________________________________________________________
    % Written by Kevin L. Turner
    % The Pennsylvania State University, Dept. of Biomedical Engineering
    % https://github.com/KL-Turner
    %________________________________________________________________________________________________________________________
    %
    %   Purpose: Add existing labels to (potentially) re-processed model tables
    %________________________________________________________________________________________________________________________


   % I made this changes because in some data sets I trim more data than
   % usual 
    % Define how much data was removed
    timeRemovedStart = 800; % seconds removed at the start
    timeRemovedEnd = 1;     % seconds removed at the end
    timeResolution = 5;     % sampling resolution in seconds per data point

    for a = 1:size(procDataFileIDs, 1)
        procDataFileID = procDataFileIDs(a, :);
        modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
        trainingDataSetID = [procDataFileID(1:end-12) 'TrainingData.mat'];
        load(modelDataSetID)
        load(trainingDataSetID)
        disp(['Updating training data set for ' trainingDataSetID '...' ]); disp(' ')
        
        % Original behavioral state labels
        behaveState = trainingTable.behavState;

        % Calculate the number of data points removed
        removedStartIdx = ceil(timeRemovedStart / timeResolution);
        removedEndIdx = ceil(timeRemovedEnd / timeResolution);

        % Adjust `behaveState` to match the new size of the data
        if length(behaveState) > (removedStartIdx + removedEndIdx)
            behaveState = behaveState(removedStartIdx + 1:end - removedEndIdx);
        else
            error('Removed too many data points; `behaveState` is too short.');
        end

        % Ensure the adjusted `behaveState` matches the height of `paramsTable`
        paramsTable = trainingTable;
        if length(behaveState) ~= size(paramsTable, 1)
            error('Mismatch between `behaveState` and `paramsTable` after adjustment.');
        end
        paramsTable.behavState = behaveState;

        % Save the updated table
        trainingTable = paramsTable;
        save(trainingDataSetID, 'trainingTable')
    end
end
