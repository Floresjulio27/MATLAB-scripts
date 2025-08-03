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

% Loop until user confirms ROI selection
roiConfirmed = false;
while ~roiConfirmed
    % Display the first frame and allow user to draw a polygon ROI
    figure, imshow(grayFrame1), title('Click to draw ROI, Press Enter to finish');
    % mask = drawrectangle; %roipoly; % Let the user draw a polygon

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
    
    % Ask user for confirmation
    choice = questdlg('Are you okay with the selected ROI?', 'Confirm ROI', 'Yes', 'No', 'Yes');
    
    if strcmp(choice, 'Yes')
        roiConfirmed = true;
    else
        close all; % Close previous figures and reselect
    end
end

close all; % Close ROI confirmation figure

% Save the ROI coordinates for future use
roiCoordinates = position; % [x, y, width, height]
% save('roi_coordinates.mat', 'roiCoordinates');
%% Convert to NaN-based mask
nanMask = double(mask);
nanMask(nanMask == 0) = NaN; % Outside the selected region is NaN

%%
% Initialize result matrix
frameCount = video.NumFrames;
% Preallocate variables to store results for each batch
maskedFrames = NaN(size(grayFrame1, 1), size(grayFrame1, 2), frameCount); % Preallocate 3D matrix for masked frames
differences = NaN(size(grayFrame1, 1), size(grayFrame1, 2), frameCount - 1); % Store frame differences

% Define the batch size (1/10 of total frames)
batchSize = floor(frameCount / 10);
remainingFrames = mod(frameCount, 10);

% Open video file
video = VideoReader(videoFile);

% Process frames in batches
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
            
            % Apply the mask
            maskedFrame = double(grayFrame) .* nanMask;
            batchMaskedFrames(:, :, i - startFrame + 1) = maskedFrame; % Store in batch array
            
            % Compute absolute difference between consecutive frames
            if i > startFrame
                batchDifferences(:, :, i - startFrame) = abs(batchMaskedFrames(:, :, i - startFrame + 1) - batchMaskedFrames(:, :, i - startFrame));
            end
        end
    end
    
    % Store results for the current batch in the global arrays
    maskedFrames(:, :, startFrame:endFrame) = batchMaskedFrames;
    differences(:, :, startFrame:endFrame-1) = batchDifferences;
end
%% Display a comparison of the first two frames
% figure;
% subplot(1, 3, 1), imshow(maskedFrames(:, :, 11), []), title('Frame 1 (Masked)');
% subplot(1, 3, 2), imshow(maskedFrames(:, :, 12), []), title('Frame 2 (Masked)');
% subplot(1, 3, 3), imshow(differences(:, :, 11), []), title('Difference (Frame 2 - Frame 1)');

%% Plot masked frames and differences (showing only inside the ROI)
    frames_to_plot = 1;
    % Get the ROI boundaries
    x1 = round(roiCoordinates(1));
    y1 = round(roiCoordinates(2));
    w = round(roiCoordinates(3));
    h = round(roiCoordinates(4));
    
    % Extract the region inside the ROI for the masked frame
    croppedFirstFrame = maskedFrames(y1:y1+h, x1:x1+w, frames_to_plot); 
    croppedSecondFrame = maskedFrames(y1:y1+h, x1:x1+w, frames_to_plot + 1); 
    
    % Extract the region inside the ROI for the difference frame
    croppedDifference = differences(y1:y1+h, x1:x1+w, frames_to_plot); 
    
    figure;
    subplot(1, 3, 1);
    imshow(croppedFirstFrame, []);
    title(['Frame ' num2str(frames_to_plot) ' (Zoomed ROI)']);

    subplot(1, 3, 2);
    imshow(croppedSecondFrame, []);
    title(['Frame ' num2str(frames_to_plot+1) ' (Zoomed ROI)']);
    
    subplot(1, 3, 3);
    imshow(croppedDifference, []);
    title(['Difference ' num2str(frames_to_plot) ' (Zoomed ROI)']);
    saveas(gcf,'ExampleCalculation.png')
    close all;
%%

% Compute mean intensity difference for motion detection
meanIntensityDifference = nanmean(differences, [1 2]); % Mean over height & width
meanIntensityDifference = meanIntensityDifference(:);

% Display motion detection results
figure;
frameRate = 15;
plotTime = (1:length(meanIntensityDifference))/frameRate; 
plot(plotTime,meanIntensityDifference, 'LineWidth', 1);
xlabel('Time(s)');
ylabel('Mean Intensity Difference');
title('Mean Intensity Difference');
saveas(gcf,'Meanintensity_Change.png')
close all;
% Define output file
filenameStr = strfind(filename,'.');
Nfilename = filename(1:filenameStr(1)-1);
outputFilePath = fullfile(pathname, Nfilename);
outputMatFile = [outputFilePath '_processed_data.mat'];

% Save the results
save(outputMatFile, 'differences', 'meanIntensityDifference', 'roiCoordinates');

disp(['Processed data saved to: ', outputMatFile]);
%% Define output video file
outputVideoFile = [outputFilePath '_processed_frame_differences.avi'];
outputVideo = VideoWriter(outputVideoFile, 'Uncompressed AVI'); % or 'Grayscale AVI'
outputVideo.FrameRate = video.FrameRate; % Set frame rate to match original video
open(outputVideo);

% Process each frame difference and write to video
for i = 1:length(differences(1,1,:))
    if ~isempty(differences(:,:,i))
        % Compute absolute change in pixel intensity
        absDifference = abs(differences(:,:,i)); 

        % Normalize for better visualization (optional)
        absDifference = uint8(255 * mat2gray(absDifference)); 

        % Write the frame to video
        writeVideo(outputVideo, absDifference);
    end
end

% Close the video file
close(outputVideo);

disp(['Difference video saved to: ', outputVideoFile]);