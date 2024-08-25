function v = mismatchL2(stimParams, rawImageFolder, fValsToUse, RFdata, elecIdx)
% Returns the HSV mismatch (L2 Norm / Euclidean Distance) between the IN and OUT portions of all images for a particular electrode.
% Currently designed to work with <displaySCNI.m>. Assumes the IN HSV values have already been computed and skips their computation.
    circ = @(X, Y, h, k, r) (((X - h).^2 + (Y - k).^2) <= r.^2); % Anon. function for circular boundaries
    [xAxisDeg, yAxisDeg] = getMonitorDetails;
    [X, Y] = meshgrid(xAxisDeg, yAxisDeg(end:-1:1)); % Applying Y axis flip to convert from matrix to image convention
    
    IN = zeros(length(stimParams), 3); % Initializing matrix to store IN HSV values
    OUT = IN; % Initializing matrix to store OUT HSV values
    for i = 1:size(IN, 1)
        IN(i, :) = [stimParams{i}.hueDeg/360, stimParams{i}.sat, ((stimParams{i}.contrastPC*stimParams{i}.spatialFreqPhaseDeg)/1.8e4) + 0.5]; % IN HSV
        img = rgb2hsv(imread(fullfile(rawImageFolder, sprintf("Image%d.png", fValsToUse(i)))));
        RFcenter = [RFdata.rfStats(RFdata.highRMSElectrodes(elecIdx)).meanAzi,RFdata.rfStats(RFdata.highRMSElectrodes(elecIdx)).meanEle]; % # RF center (dva)
        r = stimParams{i}.radiusDeg; % # RF radius (dva)
        maskIN = circ(X, Y, RFcenter(1), RFcenter(2), r);
        H = img(:, :, 1); H = H(~maskIN); % Isolating OUT H values
        S = img(:, :, 2); S = S(~maskIN);
        V = img(:, :, 3); V = V(~maskIN);
        X_ = S.*cos(H); Y_ = S.*sin(H); % Polar to Cartesian coordinate conversion
        x = mean(X_);
        y = mean(Y_);
        OUT(i, :) = [(mod(2*pi + atan2(y, x), 2*pi))*(0.5*pi), (x^2 + y^2)^0.5, mean(V)]; % OUT HSV
    end
    v = sqrt(sum((IN - OUT).^2, 2)); % L2 Norm / Euclidean Distance
end

