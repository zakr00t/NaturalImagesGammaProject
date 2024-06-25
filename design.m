function [X_tensor, Y] = design(idx, radiusMatrixDeg, patchSizeDeg, powerOption)
% Returns the HSV design matrix (tensor) and target matrix/vector for an images dataset
    
    % Startup:
    [experimentalDetails, ~] = getExperimentalDetails;
    subjectName = experimentalDetails{idx}{1};
    imageFolderName = experimentalDetails{idx}{4};
    expDate = experimentalDetails{idx}{5};
    protocolName = experimentalDetails{idx}{6};

    folderSourceString = fileparts(pwd);
    rawImageFolder = fullfile(folderSourceString, 'data', 'images', imageFolderName);
    
    if ~exist('radiusMatrixDeg', 'var'); radiusMatrixDeg = 0.3; end
    if ~exist('patchSizeDeg', 'var'); patchSizeDeg = 2;         end
    if ~exist('powerOption', 'var'); powerOption = 3;         end

    plottingDetails.displayPlotsFlag = 0;
    % [Y, ~, electrodeList] = getMeanEnergy(subjectName, expDate, protocolName);
    [powerST, powerBL, electrodeList] = getMeanEnergy(subjectName, expDate, protocolName);
    if powerOption == 1
        Y = powerST; % Only take stimulus power
    elseif powerOption == 2
        Y = powerST./powerBL; % Ratio between ST and BL
    elseif powerOption == 3
        [powerST2, powerBL2, ~] = getMeanEnergy(subjectName, expDate, protocolName, '', {[80 150]}); % Take power between 80 to 150 Hz
        powerST = powerST./powerBL; % Ratio between ST and BL
        powerST2 = powerST2./powerBL2; % Ratio between ST and BL in high gamma
        Y = powerST - powerST2;
    end

    Y = squeeze(Y)'; % Matrix of actual stimulus gamma energies for all images and electrodes.
    numImages = size(Y, 1);
    numElectrodes = size(Y,  2);
    
    % Populating X_tensor (design tensor with stimParams as features):
    X_tensor = zeros(numImages, 6, numElectrodes);
    X_tensor(:, 1, :) = 1; % Ones vector, to be multiplied with c0
    for i = 1:numImages
        imageFileName = fullfile(rawImageFolder,['Image' num2str(i) '.png']);
        [patchData, imageAxesDeg] = getImagePatches(imageFileName, electrodeList, ...
            subjectName, folderSourceString, patchSizeDeg, plottingDetails);
        % Get Stimulus Parameters:
        for j = 1:numElectrodes
            stimParams = getSingleImageParameters(rgb2hsv(patchData{j}), imageAxesDeg, ...
                [0, 0], radiusMatrixDeg, [], 0);
            X_tensor(i, 2, j) = deg2rad(stimParams.hueDeg); % H
            X_tensor(i, 4, j) = stimParams.sat; % S
            X_tensor(i, 5, j) = 0.5*(1 + sin(deg2rad(stimParams.spatialFreqPhaseDeg))...
                *(stimParams.contrastPC/100)); % V
            X_tensor(i, 6, j) = stimParams.radiusDeg; % R
        end
    end
    X_tensor(:, 3, :) = sin(X_tensor(:, 2, :)); % sinH
    X_tensor(:, 2, :) = cos(X_tensor(:, 2, :)); % cosH
    % X_tensor(:, 3, :) = X_tensor(:, 2, :);
    % X_tensor(:, 2, :) = (cos(X_tensor(:, 2, :)) + 1)/2; % cosH with 0-1 normalization
    % X_tensor(:, 3, :) = (sin(X_tensor(:, 3, :)) + 1)/2; % sinH with 0-1 normalization
end