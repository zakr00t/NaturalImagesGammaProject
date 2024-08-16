function [correlationsFull, predictionString, correlationsSelected, numSelectedImages] = ...
    analyzeData(subjectName,expDate,protocolName,imageFolderName,imageIndices,versionFlag,patchSizeDeg,radiusMatrixDeg,selectOptions,powerOption,folderSourceString)

    % Input Arguments:
    if ~exist('versionFlag','var');         versionFlag = 0;                end
    if ~exist('patchSizeDeg','var');        patchSizeDeg = 2;               end
    if ~exist('radiusMatrixDeg','var');     radiusMatrixDeg=[];             end
    if isempty(radiusMatrixDeg)
        radiusMatrixDeg = 0.3:0.3:patchSizeDeg;
    else
        patchSizeDeg = max(patchSizeDeg, max(radiusMatrixDeg));
    end
    if ~exist('selectOptions','var');       selectOptions=[];               end
    if isempty(selectOptions)
        selectOptions.meanThr = 5e-2*ones(1, 3);
        selectOptions.stdThr = 2*selectOptions.meanThr;
        selectOptions.measure = 'diff';
        selectOptions.method = 'vector';
    end
    if ~exist('powerOption','var');         powerOption=3;                  end
    if ~exist('folderSourceString','var');  folderSourceString = '';        end
    if isempty(folderSourceString)
        folderSourceString = fileparts(pwd);
    end
    
    gridType = 'Microelectrode';
    RFdata = load(fullfile(folderSourceString,'data','rfData',subjectName,[subjectName gridType 'RFData.mat']));
    rawImageFolder = fullfile(folderSourceString,'data','images',imageFolderName);
    
    % 1. Get stimulus parameters from image patches:
    disp('Getting stim params...');
    electrodeList = RFdata.highRMSElectrodes;
    numImages = length(imageIndices);
    numElectrodes = length(electrodeList);
    allStimParams = cell(numImages,numElectrodes);
    plottingDetails.displayPlotsFlag=0; % Necessary argument in getImagePatches
    for i=1:numImages
        % Load image
        imageFileName = fullfile(rawImageFolder,['Image' num2str(imageIndices(i)) '.png']);
        [patchData,imageAxesDeg] = getImagePatches(imageFileName,electrodeList,subjectName,folderSourceString,patchSizeDeg,plottingDetails);
        % Get Stim Parameters
        for j=1:numElectrodes
            stimParams = getSingleImageParameters(rgb2hsv(patchData{j}),imageAxesDeg,[0 0],radiusMatrixDeg,selectOptions,0);
            allStimParams{i,j} = stimParams;
        end
    end
    
    % 2. Get actual gamma power:
    [powerST,powerBL] = getMeanEnergy(subjectName,expDate,protocolName);
    if powerOption==1
        powerST = squeeze(powerST(:,:,imageIndices)); % Only take stimulus power
    elseif powerOption==2
        powerST = squeeze(powerST(:,:,imageIndices)) ./ squeeze(powerBL(:,:,imageIndices)); % Ratio between ST and BL
    elseif powerOption==3
        [powerST2,powerBL2] = getMeanEnergy(subjectName,expDate,protocolName,'',{[80 150]}); % Take power between 80 to 150 Hz
        powerST = squeeze(powerST(:,:,imageIndices)) ./ squeeze(powerBL(:,:,imageIndices)); % Ratio between ST and BL
        powerST2 = squeeze(powerST2(:,:,imageIndices)) ./ squeeze(powerBL2(:,:,imageIndices)); % Ratio between ST and BL in high gamma
        powerST = powerST - powerST2;
    end
    
    % 3. Get predicted gamma power and correlations:
    if ~versionFlag
        correlationsFull = zeros(6,numElectrodes);
        correlationsSelected = zeros(6,numElectrodes);
        numSelectedImages = zeros(1,numElectrodes);
        for j=1:numElectrodes
            stimParams = allStimParams(:,j);
            actualPower = powerST(j,:);
            [correlationsFull(:,j),correlationsSelected(:,j),~,~,selectedImageIndices] = getAllCorrelations(subjectName,stimParams,actualPower);
            numSelectedImages(j) = length(selectedImageIndices);
            predictionString = ["H", "S", "V", "HS", "HSV", "HSVR"];
        end
    else
        correlationsFull = zeros(7, numElectrodes);
        for j = 1:numElectrodes
            stimParams = allStimParams(:, j);
            actualPower = powerST(j, :);
            correlationsFull(:,j) = getAllCorrelations(subjectName, stimParams, actualPower, [], versionFlag, j, RFdata, rawImageFolder, imageIndices);
        end
        predictionString = ["H", "S", "V", "HS", "HSV", "P", "HSV+P"];
        correlationsSelected = [];
        numSelectedImages = [];
    end
end