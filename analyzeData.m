
[experimentalDetails,matchIndex] = getExperimentalDetails;
posList = 5; % Index for which data needs to be saved

for i=1:length(posList)
    tmp = experimentalDetails{posList(i)};
    
    % 1. Get actual gamma power
    [powerST,powerBL,electrodeList,rfStats,psdST,psdBL,freqVals] = getMeanEnergy(tmp{1},tmp{5},tmp{6});
    powerST = squeeze(powerST); powerBL = squeeze(powerBL);
    
    % 2. Get Predicted gamma power
    numElectrodes = size(psdST,1);
    numImages = size(psdST,2);
    
    predictedPower = zeros(numElectrodes,numImages);
    for j=1:numImages
        % Load image
        
        for k=1:numElectrodes
            % Use rfStats to get rfCenter and then to find the appropriate patch
            % For this patch, use getSingleImageParameters to get the appropriate stimulus properties 
            % Now use getPredictedGamma to get the predicted Gamma and
            % populate predictedPower
        end
    end
    
    % For each electrode, find the correlation
    corrVals = zeros(1,numElectrodes);
    for j=1:numElectrodes
        % Get correlation between actual and predicted gamma. Take care of
        % the NaNs for images for which no prediction is possible. Can also
        % use a cutoff so that correlation is only done for electrodes for
        % which power can be predicted for sufficient images
    end
end