%% PROGRAM FOR VEHICLE NUMBER PLATE RECOGNITION USING STANDARD CHARACTER TEMPLATES- Ashwin Vinoo Version

clear;               % Clears the workspace
pause on;            % Enables the pause function
close all force;     % Closes all figures
clc;                 % Clears the console

%% ----------- Hyper Parameters ----------
% Update the working directory when changing laptops
workingDirectory = 'C:\Users\Ashwin Vinoo\Desktop\Work\License Plate Recognition';
% Set True to generate the standard character templates
createNewTemplates = true;
% This is the size of the character template generated (rows & columns)
characterTemplateSize = [42, 21];
% Specify the font used to create the standard character templates
CharacterFont = 'License Plate Regular';
% Set True to generate the standard character curve profiles
createNewCurveProfiles = true;
% Sets the threshold of the RGB color Filter
colorThreshold = 80;
% The sensitivity for the adaptive thresholding process
adaptiveThresholdSensitivity = 0.35;
% The minimum and maximum number of pixels accepted in an character object
objectPixelArea = [200, 3000];
% The minimum accepted correlation of a region with the character templates
minimumCharacterCorr = 0.40;
% The deviation accepted in pixels from the detected plate number line
centralPlateDeviation = 30;

%% ----------- Creating the standard Character Templates ----------
% make a directory to store the standard templates generated
mkdir(strcat(workingDirectory,'\Standard Character Templates')); 
% The list of characters needed for license plate recognition
characterList = strcat('A':'Z','0':'9');
% Creating a cell array to save the standard templates for fast loading into the workspace
standardTemplates = cell(1,36);
% Create the hyperparameters if requested
if(createNewTemplates)
    % Will loop through each element in this list and create the corresponding character template image
    for i = 1:length(characterList)
        % Sets the figure background color to black
        figure('Color',[0 0 0]);
        % Creates the alphanumeric character on a figure with the specified properties
        text(0,0.5,characterList(i),...
        'FontName',CharacterFont,...
        'Color', 'white',...
        'FontUnits', 'pixels',...
        'FontSize', 500,...
        'FontWeight', 'normal',...
        'BackgroundColor', 'black',...
        'Margin',1)
        % Removes the axis from the figure so that it does not appear in the final image
        axis off;
        % Gets the frame from the current image
        im = getframe(gcf);

        % Extracts the RGB image from the frame obtained
        rgb_image = im.cdata;
        % Converts the RGB image into a grayscale
        grayscale_image = rgb2gray(rgb_image);
        % Converts the grayscale image to binary
        binary_image = imbinarize(grayscale_image);

        % Returns all the rows - column combinations having 
        [rows, columns] = find(binary_image);
        % Finds the topmost row containing atleast one white pixel
        topRow = min(rows);
        % Finds the topmost row containing atleast one white pixel
        bottomRow = max(rows);
        % Finds the leftmost containing atleast one white pixel
        leftColumn = min(columns);
        % Finds the rightmost containing atleast one white pixel
        rightColumn = max(columns);
        % Crops the image so that the images borders are trimmed till the character
        cropped_image = binary_image(topRow:bottomRow, leftColumn:rightColumn);
        % Closes the current figure
        close(gcf);
        % resize this image to that needed for the standard templates
        standardTemplate = imresize(cropped_image,characterTemplateSize);
        % In the case of 'I' add padding so that correlation doesn't fail later on
        if(characterList(i) == 'I')
            standardTemplate(:,1) = 0;
            standardTemplate(:,characterTemplateSize(2)) = 0;
        end;
        % Saves each standard template into the cell array
        standardTemplates{i} = standardTemplate;
        % Write the image over to the Standard templates directory
        imwrite(standardTemplate,strcat(workingDirectory,...
            '\Standard Character Templates\',characterList(i),'.bmp'),'bmp');
    end;
end;
% Saves all the standard templates as a mat file along with the images
save(strcat(workingDirectory, '\Standard Character Templates\',...
            'standardTemplates'),'standardTemplates');
    
%% ----------- Performing license plate recognition on each of the sample plates ----------

% Getting data from the the excel which contains the actual values of the license plates
[~,plateData] = xlsread(strcat(workingDirectory,...
    '\Sample License Plates\Sample License Plates.xlsx'),'Samples','B3:C27');

% Loads all the standard character template images from the mat file
load(strcat(workingDirectory, '\Standard Character Templates\','standardTemplates'));

% Variable used to count the number of plates detected correctly
platesDetectedCorrectly = 0;

% Cycling through each license plate using a for loop
for i = 1:length(plateData)
    
    % Reading the license plate as an RGB image
    currentLicensePlate_rgb = imread(strcat(workingDirectory,...
        '\Sample License Plates\',char(plateData(i,1)),'.jpg'),'jpg');  
    
    % Resizing the RGB image to 500 pixels height while maintaining the same aspect ratio
    currentLicensePlate_rgb = imresize(currentLicensePlate_rgb,[500 NaN]);
    
    % Obtaining the new image dimensions
    [imageRows, imageColumns,~]=size(currentLicensePlate_rgb);  
    
    % Obtaining the red, green and blue components of the image
    red=currentLicensePlate_rgb(:,:,1);
    green=currentLicensePlate_rgb(:,:,2);
    blue=currentLicensePlate_rgb(:,:,3);
    
    % Color Filtering - Any pixel which deviates too much from a shade of gray is made white
    colorFilter = false(imageRows,imageColumns);
    for j = 1:imageRows
        for k = 1:imageColumns
            if(abs(red(j,k)-green(j,k)) + abs(green(j,k)-blue(j,k)) + abs(blue(j,k)-red(j,k)) > colorThreshold)
                colorFilter(j,k) = 1;
            end;    
        end;
    end;
    
    % Converting the RGB image to gray scale (shades of black and white)
    currentLicensePlate_gray = rgb2gray(currentLicensePlate_rgb); 
    
    % Blurring the image with a gaussian filter to remove noise
    currentLicensePlate_gray = imgaussfilt(currentLicensePlate_gray,0.5,'FilterSize',5);
    
    % Sharpening the image so that edges weakened by blurring get resolved
    currentLicensePlate_gray = imsharpen(currentLicensePlate_gray,'Radius',1,'Amount',0.8);
    
    % Image adjustment so that contrast is improved. Useful when brightness is inconsistent
    currentLicensePlate_gray = imadjust(currentLicensePlate_gray,[0.05 0.95],[0 1]);
    
    % Creating an edge mask
    edgeMask = edge(currentLicensePlate_gray,'Canny',0.5);
    
    % Creating the structural element for dilation of the edge mask
    structuralElement = strel('disk',1);
    
    % Dilating the edge mask using the structural element
    edgeMask = imdilate(edgeMask, structuralElement);
    
    % Adaptive thresholding ensures that bright plate number sections are retained during binarization
    adaptiveThreshold = adaptthresh(currentLicensePlate_gray, adaptiveThresholdSensitivity);
    
    % Converting the image to black and white format
    currentLicensePlate_bw = imbinarize(currentLicensePlate_gray, adaptiveThreshold);
    
    % Using the color filter primarily helps problems with yellow license plates
    currentLicensePlate_bw = max(currentLicensePlate_bw, colorFilter);
    
    % Using the edge mask to ensure that the dark license plate background
    % gets separated from the text in case thresholding failed to do so
    currentLicensePlate_bw = max(currentLicensePlate_bw, edgeMask);
    
    % Now we are inverting the image so that object identification can take place
    currentLicensePlate_bw = imcomplement(currentLicensePlate_bw);
    
    imwrite(imcomplement(currentLicensePlate_bw),strcat('C:\Users\Ashwin Vinoo\Desktop\img',num2str(i),'.png'),'png');
    
    % Filtering the image for interconnected objects within the specified pixel range
    filteredLicensePlate_bw = bwareafilt(currentLicensePlate_bw ,objectPixelArea);
    
    % Identifying all the different interconnected components from the image (islands of 1's)
    connectedRegions = regionprops(filteredLicensePlate_bw,'BoundingBox');
    
    % Obtaining the same in a matrix format
    connectedRegions = cat(1,connectedRegions.BoundingBox);
    
    % Creating cells to hold data on regions and detected characters
    regionCharacters = cell(length(connectedRegions),3);
    
    % Creating a variable to help find the row number of the plate numbers
    weightedAverageRow = 0;
    
    % Looping through the different interconnected regions within the image
    for j = 1:length(connectedRegions)
        
        % We will process the image regions only if the width-height ratios are within range
        if(connectedRegions(j,3) < 2 * connectedRegions(j,4) &&...
           connectedRegions(j,4) < 10 * connectedRegions(j,3) &&...
           connectedRegions(j,1) > imageColumns/20 &&...
           connectedRegions(j,1) +  connectedRegions(j,3) < imageColumns * 19/20 &&...
           connectedRegions(j,2) > imageRows/20 &&...
           connectedRegions(j,2) +  connectedRegions(j,4) < imageRows * 19/20)
       
            % Cropping the current region out of the image         
            croppedRegion = imcrop(filteredLicensePlate_bw,connectedRegions(j,:));
            
            % Resizing the cropped region to that used by the standard character templates
            croppedRegion = imresize(croppedRegion, characterTemplateSize);
            
            % Create an array to store the correlation with the standard templates 
            correlationArray = zeros(1,length(standardTemplates));
            
            % Looping through the different character templates to find the best match
            for k=1:length(standardTemplates)
                % Find the correlation of the current region with current template 
                correlationArray(k)= corr2(standardTemplates{k}, croppedRegion); 
            end
            
            % This is the maximum correlation detected for any template
            maxCorr = max(correlationArray);
            
            % We only accept that a character has been detected if a minimum correlation is satisfied
            if(maxCorr >= minimumCharacterCorr)
                
                % The index of the character template that has this correlation
                templateIndex = find(correlationArray == maxCorr, 1);
                
                % Storing the index number of the region
                regionCharacters{j,1} = j;
                
                % Storing the detected character in the character list
                regionCharacters{j,2} = characterList(templateIndex);   
                
                % Storing the detected character in the character list
                regionCharacters{j,3} = maxCorr;   

                % This variable will help find the mean row on which the plate numbers lie
                weightedAverageRow = weightedAverageRow + ...
                    maxCorr * (connectedRegions(j,2) + connectedRegions(j,4)/2);       
            end;
        end;
    end;
    
    % Removing all the empty rows from the cell array
    regionCharacters(any(cellfun(@isempty,regionCharacters),2),:) = []; 
    
    % Determines the most probable row on which the plate numbers lie
    plateNumberRow = weightedAverageRow/sum(cellfun(@double,regionCharacters(:,3)));
    
    % Creating an empty string upon which the license plate letters will be added
    detectedPlateNumber = string('');
    
    % Looping through the different accepted regions
    for j = 1:length(regionCharacters)
        
        % This is the region number assigned earlier
        regionNumber = regionCharacters{j,1};
        
        % Obtaining the starting row of the current region
        regionStartingRow = connectedRegions(regionNumber, 2);
        
        % Obtaining the starting row of the current region
        regionEndingRow = connectedRegions(regionNumber, 2) + connectedRegions(regionNumber, 4);
        
        % Removing any region that isnt near the determined plate number row
        if(regionStartingRow < plateNumberRow + centralPlateDeviation &&...
           regionEndingRow > plateNumberRow -  centralPlateDeviation)
       
            % Adding the detected character to the plate number string
            detectedPlateNumber = strcat(detectedPlateNumber, regionCharacters{j,2});
        end;
    end;
    
    imtool(filteredLicensePlate_bw);
    
    % Writes the detected license plate to the excel where sample information is stored
    xlswrite(strcat(workingDirectory, '\Sample License Plates\Sample License Plates.xlsx'),...
             {char(detectedPlateNumber)}, 'Samples', strcat('D',num2str(i+2)));
         
    % Checking if the plate was detected correctly
    if(strcmp(char(plateData(i,2)),detectedPlateNumber))   
        % If the plate was detected correctly add one to the count of plates correctly detected
        platesDetectedCorrectly = platesDetectedCorrectly+1;
    end;
end;
    
% Displays the overall license plate accuracy based on how many plates were correctly identified
disp(strcat('Overall License Plate Detection Accuracy: ',...
    num2str(platesDetectedCorrectly/length(plateData)*100), '%'));

%---------- END OF CODE ----------