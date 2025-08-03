clc; clear; close all;

% Open file selection dialog
[filename, pathname] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, 'Select a Video File');
if filename == 0
    disp('No file selected. Exiting...');
    return; % Exit if no file is selected
end
videoFile = fullfile(pathname, filename);
video = VideoReader(videoFile);

% Read the first frame
frame1 = readFrame(video);
grayFrame1 = rgb2gray(frame1); % Convert to grayscale

% Initialize a cell array to store ROI coordinates
roiCoordinatesList = {};

%% Loop for selecting multiple ROIs
roiConfirmed = false;
while ~roiConfirmed
    % Display the first frame for ROI selection
    figure, imshow(grayFrame1), title('Draw ROI and double-click inside to confirm');
    
    % User selects a rectangular ROI
    h = imrect; 
    position = wait(h); % Wait for user to finish selection
    
    % Create a binary mask from the selected rectangle
    mask = false(size(grayFrame1));
    mask(round(position(2)):(round(position(2)) + round(position(4))), ...
         round(position(1)):(round(position(1)) + round(position(3)))) = true;
    
    % Show the selected ROI
    maskedROI = grayFrame1;
    maskedROI(~mask) = 0; % Set everything outside ROI to black for visualization
    
    figure, imshow(maskedROI), title('Selected ROI (Confirm or Redo)');
    
    % Save the ROI coordinates for future use
    roiCoordinatesList{end+1} = position; % Append new ROI
    
    % Ask user for confirmation or to select another ROI
    choice = questdlg('Would you like to select another ROI?', 'Select ROI', 'Yes', 'No', 'No');
    
    if strcmp(choice, 'No')
        roiConfirmed = true;
    else
        close all; % Close previous figures and reselect
    end
end

% Close ROI selection figure
close all;
%% Initialize result matrix
frameCount = video.NumFrames;
% Define the batch size (1/10 of total frames)
batchSize = floor(frameCount / 10);
remainingFrames = mod(frameCount, 10);

% Process each ROI separately
for roiIdx = 1:length(roiCoordinatesList)
    roiCoordinates = roiCoordinatesList{roiIdx};
    
    % Convert the ROI to NaN-based mask
    nanMask = false(size(grayFrame1));
    x1 = round(roiCoordinates(1));
    y1 = round(roiCoordinates(2));
    w = round(roiCoordinates(3));
    h = round(roiCoordinates(4));
    nanMask(y1:y1+h, x1:x1+w) = true; % Apply mask to selected region
    nanMask = double(nanMask);
    nanMask(nanMask == 0) = NaN; % Outside the selected region is NaN
    
    % only do this if you want to save the global arrays. I am not going to
    % do it because of memory.

    % % Preallocate matrices for results specific to this ROI
    % maskedFrames = NaN(size(grayFrame1, 1), size(grayFrame1, 2), frameCount);
    % differences = NaN(size(grayFrame1, 1), size(grayFrame1, 2), frameCount - 1);

    % Initialize video reader
    video = VideoReader(videoFile);
    
    % Preallocate mean intensity difference array
    meanIntensityDifferences = NaN(1, frameCount - 1);

    %% Process frames in batches
    for batchIdx = 1:10
        % Determine the batch of frames to process
        startFrame = (batchIdx - 1) * batchSize + 1;
        endFrame = min(batchIdx * batchSize, frameCount);
        
        % Initialize temporary arrays for this batch
        batchMaskedFrames = NaN(size(grayFrame1, 1), size(grayFrame1, 2), endFrame - startFrame + 1);
        batchDifferences = NaN(size(grayFrame1, 1), size(grayFrame1, 2), endFrame - startFrame);
        
        % Set video reader to the start frame for the batch
        video.CurrentTime = (startFrame - 1) / video.FrameRate; % Set the start frame
        
        % Process frames in the current batch
        for i = startFrame:endFrame
            if hasFrame(video)
                frame = readFrame(video);
                grayFrame = rgb2gray(frame);
                
                % Apply the mask for the current ROI
                maskedFrame = double(grayFrame) .* nanMask;
                batchMaskedFrames(:, :, i - startFrame + 1) = maskedFrame; % Store in batch array
                
                % Compute absolute difference between consecutive frames
                if i > startFrame
                    batchDifferences(:, :, i - startFrame) = abs(batchMaskedFrames(:, :, i - startFrame + 1) - batchMaskedFrames(:, :, i - startFrame));
                    
                    % Calculate mean pixel intensity difference for the ROI
                    meanIntensityDifferences(i - 1) = mean(batchDifferences(:,:,i - startFrame), 'all', 'omitnan');
                end
            end
        end
        % we will not store the global arrays to reduce memory use. But you
        % can do that using the following code

        % % Store results for the current batch in the global arrays
        % maskedFrames(:, :, startFrame:endFrame) = batchMaskedFrames;
        % differences(:, :, startFrame:endFrame-1) = batchDifferences;
    end
    
   
    % Interpolate NaN values in meanIntensityDifferences
    validIndices = find(~isnan(meanIntensityDifferences)); % Indices with valid values
    if length(validIndices) > 1
        % Interpolate the NaN values linearly
        meanIntensityDifferences = interp1(validIndices, meanIntensityDifferences(validIndices), 1:(frameCount-1), 'linear', 'extrap');
    end
    %% Define output file
    filenameStr = strfind(filename,'.');
    Nfilename = filename(1:filenameStr(1)-1);
    outputFilePath = fullfile(pathname, Nfilename);
    outputIntensityPlotName = [outputFilePath 'ROI_' num2str(roiIdx) '_Meanintensity_Change.png'];
    outputDiffrencePlotName = [outputFilePath 'ROI_' num2str(roiIdx) '_ExampleCalculation.png'];
    outputMatFileName = [outputFilePath 'ROI_' num2str(roiIdx) '_Meanintensity_calculation.mat'];
    %% Example of plotting the mean pixel intensity difference
    figure;
    FrameRate = 15;
    plotTime = (1:length(meanIntensityDifferences))/FrameRate; 
    plot(plotTime, meanIntensityDifferences, 'LineWidth', 1);
    title(['Mean Pixel Intensity Difference for ROI ' num2str(roiIdx)]);
    xlabel('Time (seconds)');
    ylabel('Mean Intensity Difference');    
    saveas(gcf,outputIntensityPlotName)
    close all
    %% Example of plotting the first few frames for this ROI
    frames_to_plot = 1;

    % Extract the region inside the ROI for the masked frame
    croppedFrame = batchMaskedFrames(y1:y1+h, x1:x1+w, frames_to_plot); 
    cropped2ndFrame = batchMaskedFrames(y1:y1+h, x1:x1+w, frames_to_plot+1); 
    
    % Extract the region inside the ROI for the difference frame
    croppedDifference = batchDifferences(y1:y1+h, x1:x1+w, frames_to_plot); 
    
    figure;
    subplot(1, 3, 1);
    imshow(croppedFrame, []);
    title(['Frame ' num2str(frames_to_plot) ' (ROI)' num2str(roiIdx) ]);

    subplot(1, 3, 2);
    imshow(croppedFrame, []);
    title(['Frame ' num2str(frames_to_plot+1) ' (ROI)' num2str(roiIdx)]);
    
    subplot(1, 3, 3);
    imshow(croppedDifference, []);
    title(['Difference ' num2str(frames_to_plot) ' (ROI)' num2str(roiIdx)]);
    saveas(gcf,outputDiffrencePlotName)
    close all
    %% Save the results
    save(outputMatFileName,'meanIntensityDifferences');
end
