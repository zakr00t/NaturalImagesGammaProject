function [correlationsFull,correlationsSelected,numSelectedImages,predictionString] = analyzeData(subjectName,expDate,protocolName,imageFolderName,imageIndices,selectOptions,radiusMatrixDeg,folderSourceString)

if ~exist('radiusMatrixDeg','var');     radiusMatrixDeg=[];             end

patchSizeDeg = 2;
if isempty(radiusMatrixDeg)
    radiusMatrixDeg = 0.3:0.3:patchSizeDeg;
else
    patchSizeDeg = max(patchSizeDeg,max(radiusMatrixDeg));
end

plottingDetails.displayPlotsFlag=0;
    
if ~exist('folderSourceString','var');  folderSourceString = '';        end
if isempty(folderSourceString)
    folderSourceString = fileparts(pwd);
end

rawImageFolder = fullfile(folderSourceString,'data','images',imageFolderName);

% 1. Get actual gamma power
[powerST,~,electrodeList] = getMeanEnergy(subjectName,expDate,protocolName);
powerST = squeeze(powerST(:,:,imageIndices));

% 2. Get stimulus parameters for image patches
disp('Getting stim params...');
numElectrodes = size(powerST,1);
numImages = size(powerST,2);
allStimParams = cell(numElectrodes,numImages);
for i=1:numImages
    % Load image
    imageFileName = fullfile(rawImageFolder,['Image' num2str(imageIndices(i)) '.png']);
    [patchData,imageAxesDeg] = getImagePatches(imageFileName,electrodeList,subjectName,folderSourceString,patchSizeDeg,plottingDetails);
    
    % Get Stim Parameters
    for j=1:numElectrodes
        stimParams = getSingleImageParameters(rgb2hsv(patchData{j}),imageAxesDeg,[0 0],radiusMatrixDeg,selectOptions,0);
        allStimParams{j,i} = stimParams;
    end
end

% 3. Get predicted gamma power and correlations
correlationsFull = zeros(6,numElectrodes);
correlationsSelected = zeros(6,numElectrodes);
numSelectedImages = zeros(1,numElectrodes);
for i=1:numElectrodes
    actualPower = powerST(i,:);
    stimParams = allStimParams(i,:);
    [correlationsFull(:,i),correlationsSelected(:,i),predictionString,~,selectedImageIndices] = getAllCorrelations(subjectName,stimParams,actualPower);
    numSelectedImages(i) = length(selectedImageIndices);
end
end