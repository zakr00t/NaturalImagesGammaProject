% This function takes an 3-D image matrix (in HSV format) and finds a
% Grating or Hue Patch that best describes a portion of the image centered
% around rfCenter.

% To find the size of the appropriate hue patch, we need to compute mean
% and standard deviation of hue and saturation values of disks of
% progressively larger sizes centered around rfCenter. We fix the radius
% once these values deviate enough.

% Inputs
% Required:
% imageHSV: 3-D image matrix in HSV format.

% Optional
% imageAxesDeg: imageAxesDeg.xAxisDeg and imageAxesDeg.yAxisDeg describe
% the x and y ranges of the image in degrees. If these are not provided, we
% get them from getMonitorDetails. But this is only valid if the image has
% the size specified in monitorSpecifications.

% rfCenterDeg: [aziDeg eleDeg]: 2-D array that describes the position of
% the center of the RF in degrees. Default: [0 0]

% radiusMatrixDeg: 1-D matrix of radii for which image properties are
% calculated. Default: 0.3:0.3:9.6.

% selectOptions: contains the details of the thresholds that are used to
% determine the size.
% selectOptions.mean_thr: The % +/- deviation from the initial H, S, and V 
% value to be tolerated. Default: 5e-2 (5%)
% selectOptions.std_thr: Similar purpose as above, but non-negative.
% Default: 0.1 (10%, 2*mean_thr)
% selectOptions.measure = 'diff' or 'abs'. When set to diff, we compute the
% mean and std of only the pixels in the disk between the current radius
% and the preceding one. When set to 'abs', we take all pixels. Default: 
% 'diff'
% selectOptions.method: Decides whether to use an HSV cylindrical "vector" 
% based averaging or channel separated averaging ("naive") approach for the 
% patch approximation. Default = 'vector'

function stimParams = getSingleImageParameters(imageHSV,imageAxesDeg,rfCenterDeg,radiusMatrixDeg,selectOptions,displayAnalysisFlag)

if ~exist('imageAxesDeg','var');        imageAxesDeg=[];                end
if ~exist('rfCenterDeg','var');         rfCenterDeg = [0 0];            end
if ~exist('radiusMatrixDeg','var');     radiusMatrixDeg = 0.3:0.3:9.6;  end
if ~exist('displayAnalysisFlag','var'); displayAnalysisFlag=1;          end    
if ~exist('selectOptions','var');       selectOptions=[];               end

%%%%%%%%%%%%%%%%%%%%%%%%% Get selection thresholds %%%%%%%%%%%%%%%%%%%%%%%%
if isempty(selectOptions)
    selectOptions.mean_thr = 5e-2;
    selectOptions.std_thr = 2*selectOptions.mean_thr;
    selectOptions.measure = 'diff';
    selectOptions.method = 'vector';
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Get X and Y axes in degrees %%%%%%%%%%%%%%%%%%%%
if isempty(imageAxesDeg)
    [xAxisDeg,yAxisDeg,monitorSpecifications] = getMonitorDetails;
    imageXRes=size(imageHSV,2);
    imageYRes=size(imageHSV,1);
    
    if ~(isequal(monitorSpecifications.xRes,imageXRes) && isequal(monitorSpecifications.yRes,imageYRes))
        error('If image is not a screenshot, please provide imageAxesDeg');
    end
else
    xAxisDeg = imageAxesDeg.xAxisDeg;
    yAxisDeg = imageAxesDeg.yAxisDeg;
end

%%%%%%%%%%%%%%%%%% Check if the image is achromatic %%%%%%%%%%%%%%%%%%%%%%%
achromaticFlag=0;

if achromaticFlag
else
    numStimuli = length(radiusMatrixDeg);
    meanList = zeros(3,numStimuli);
    stdList = zeros(3,numStimuli);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Get the apertures %%%%%%%%%%%%%%%%%%%%%%%%%
    gaborStim.azimuthDeg=rfCenterDeg(1);
    gaborStim.elevationDeg=-rfCenterDeg(2);
    gaborStim.spatialFreqCPD=0; % These do not matter since we are only interested in the aperture
    gaborStim.sigmaDeg=100000;
    gaborStim.orientationDeg=0;
    gaborStim.contrastPC = 100;

    goodPosPreviousRadius=false(size(imageHSV,1),size(imageHSV,2)); % Keeps the goodPos of the previous radius to perform the 'diff' computation
    for i=1:numStimuli
        gaborStim.radiusDeg=radiusMatrixDeg(i);
        [~,aperture] = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg);

        % Get goodPos
        goodPos = (aperture==1);
        if strcmp(selectOptions.measure,'diff')
            goodPosToUse = xor(goodPos,goodPosPreviousRadius);
        else
            goodPosToUse = goodPos;
        end
        goodPosPreviousRadius = goodPos;
        
        for j=1:3
            tmp = squeeze(imageHSV(:,:,j));
            goodVals = tmp(goodPosToUse);
            
            if j==1 % Hue
                meanList(j,i) = circ_mean(goodVals*2*pi); % Change from 0-1 to 0-2*pi 
                stdList(j,i) = circ_std(goodVals*2*pi);
            else
                meanList(j,i) = mean(goodVals); 
                stdList(j,i) = std(goodVals);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Parameter Computation %%%%%%%%%%%%%%%%%%%%%%%% 
% The radius before the first radius at which a threshold is crossed is chosen:
r_list = [];
for j = 1:3
    foo = abs(meanList(j, :) - meanList(j, 1));
    if j == 1 % Dealing with hues in radians
        foo(foo > pi) = 2*pi - foo(foo > pi); % Ensures that reflex (>180 degree) angles are replaced with their smaller counterpart. E.g., if two vectors are separated by 270deg, it is more appropriate to say they are separated by 90deg.
        if ~isempty(min(find(foo > selectOptions.mean_thr*2*pi)))
            r_list = [r_list, min(find(foo > selectOptions.mean_thr*2*pi)) - 1];
        else
            r_list = [r_list, numStimuli]; % If all values are sub-threshold, the maximum possible radius is chosen
        end
        if ~isempty(min(find(stdList(j, :) > selectOptions.std_thr*2*pi)))
            r_list = [r_list, min(find(stdList(j, :) > selectOptions.std_thr*2*pi)) - 1];
        else
            r_list = [r_list, numStimuli];
        end
    else
        if ~isempty(min(find(foo > selectOptions.mean_thr)))
            r_list = [r_list, min(find(foo > selectOptions.mean_thr)) - 1];
        else
            r_list = [r_list, numStimuli];
        end
        if ~isempty(min(find(stdList(j, :) > selectOptions.std_thr)))
            r_list = [r_list, min(find(stdList(j, :) > selectOptions.std_thr)) - 1];
        else
            r_list = [r_list, numStimuli];
        end
    end
end
r_idx = min(r_list); % The index of the maximum approximation radius which doesn't deviate beyond threshold

% Aperture of best-fit patch:
if r_idx == 0
    error("r = 0, no suitable patch approx. possible. Consider increasing std_thr."); % While computing r for meanList, it is impossible to go below 1. This is because deviation is considered w.r.t. to the 1st value, and can only exceed the threshold as early as r = 2. For std, as early as r = 1 we can go above std_thr.
else
    gaborStim.radiusDeg = radiusMatrixDeg(r_idx);
    [~, aperture] = makeGaborStimulus(gaborStim, xAxisDeg, yAxisDeg);
end   

% Isolating HSV channel values (0-1 range) within aperture
H = squeeze(imageHSV(:, :, 1));
H = H(aperture == 1);
S = squeeze(imageHSV(:, :, 2));
S = S(aperture == 1);
V = squeeze(imageHSV(:, :, 3));
V = V(aperture == 1);

if selectOptions.method == "vector" % Vector averaging for H, S as expected in a cylindrical basis  
    N = numel(H);
    X_tot = sum((S.*cos(2*pi*H)));
    Y_tot = sum((S.*sin(2*pi*H)));
    if (X_tot >= 0) && (Y_tot >= 0) % I Quadrant
        h = atan(Y_tot/X_tot)/(2*pi);
    elseif (X_tot < 0) && (Y_tot >= 0) % II Quadrant
        h = 0.5 - atan(abs(Y_tot/X_tot))/(2*pi);
    elseif (X_tot < 0) && (Y_tot < 0) % III Quadrant
        h = 0.5 + atan(Y_tot/X_tot)/(2*pi);
    else  % IV Quadrant
        h = 1 - atan(abs(Y_tot/X_tot))/(2*pi);  
    end
    s = (((X_tot)^2 + (Y_tot)^2)^0.5)/N;
elseif selectOptions.method == "naive"
    if circ_mean(2*pi*H) >= 0
        h = circ_mean(2*pi*H)/(2*pi);
    else
        h = 1 + circ_mean(2*pi*H)/(2*pi);
    end
    s = mean(S);
end
v = mean(V);

stimParams.radiusDeg = radiusMatrixDeg(r_idx);
stimParams.hueDeg = h*360;
stimParams.saturation = s;
stimParams.value = v;

if displayAnalysisFlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    sgtitle("Image HSV Channels")
    % Image:
    subplot(2, 2, 1)
        % Steps necessary to keep the image upright and with an ascending Y axis: (REFER https://in.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab)
        imshow(flipud(hsv2rgb(imageHSV)), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        axis on
        title("RGB Image")
        xlabel("Azimuth (deg)")
        ylabel("Elevation (deg)")
    % Hue Channel:
    subplot(2, 2, 2)
        imshow(flipud(imageHSV(:, :, 1)), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        axis on
        colormap('hsv')
        caxis([0, 1])
        colormap('hsv')
        colorbar
        title("Hue")
        xlabel("Azimuth (deg)")
        ylabel("Elevation (deg)")
    % Saturation Channel:
    subplot(2, 2, 3)
        imshow(flipud(imageHSV(:, :, 2)), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        axis on
        caxis([0, 1])
        % colormap('gray') % Default
        colorbar
        title("Saturation")
        xlabel("Azimuth (deg)")
        ylabel("Elevation (deg)")
    % Value Channel
    subplot(2, 2, 4)
        imshow(flipud(imageHSV(:, :, 3)), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        axis on
        caxis([0, 1])
        colorbar
        title("Value/Brightness")
        xlabel("Azimuth (deg)")
        ylabel("Elevation (deg)")

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    sgtitle("Image Grating/Patch Approximation")
    % Image:
    subplot(2, 2, 1)
        imshow(flipud(hsv2rgb(imageHSV)), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        title(sprintf("Patch Approx. (H = %0.2g deg, S = %0.2g, V = %0.2g)", h*360, s, v))
        xlabel("Azimuth (deg)")
        ylabel("Elevation (deg)")
        % Applying the patch approximation:
        chan1 = imageHSV(:, :, 1);
        chan1(aperture == 1) = h;
        chan2 = imageHSV(:, :, 2);
        chan2(aperture == 1) = s;
        chan3 = imageHSV(:, :, 3);
        chan3(aperture == 1) = v;
        patch_approx = hsv2rgb(cat(3, chan1, chan2, chan3));
        hold on
        imshow(flipud(patch_approx), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        alpha(0.7)
        % Showing concentric circles for all radii:
        params(1) = rfCenterDeg(1);
        params(2) = rfCenterDeg(2);
        params(5) = 0;
        params(6) = 1;
        for i = 1:numStimuli
            params(3) = radiusMatrixDeg(i);
            params(4) = params(3);
            [~, ~, boundaryX, boundaryY] = gauss2D(params, xAxisDeg, yAxisDeg, []);
            if i == r_idx
                plot(boundaryX, boundaryY, '-g', 'DisplayName', sprintf("r = %0.2g", radiusMatrixDeg(r_idx))) % Highlight the patch radius
            else
                plot(boundaryX, boundaryY, '-k', 'HandleVisibility', 'off')
            end
        end
        axis on
        legend
   subplot(2, 2, 2)
        polarplot(meanList(1, :), radiusMatrixDeg, 'DisplayName', 'Mean')
        hold on
        polarplot(meanList(1, 1)*ones(numStimuli), radiusMatrixDeg, '-g', 'HandleVisibility', 'off') % Initial mean value
        polarplot((meanList(1, 1) + selectOptions.mean_thr*2*pi)*ones(numStimuli), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Upper Bound
        polarplot((meanList(1, 1) - selectOptions.mean_thr*2*pi)*ones(numStimuli), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Lower Bound
        polarplot(stdList(1, :), radiusMatrixDeg, 'DisplayName', 'Std')
        polarplot((selectOptions.std_thr*2*pi)*ones(numStimuli), radiusMatrixDeg, '--r', 'HandleVisibility', 'off') % Upper Bound on Std
        legend
        title("Hue (deg)")
        % For ease of viewing
        if numStimuli > 10 
           rticks(linspace(radiusMatrixDeg(1), radiusMatrixDeg(end), 10))
        else
           rticks(radiusMatrixDeg)
        end
        rtickformat('%0.1g')
        % xlabel("Radius (deg)")
    subplot(2, 2, 3)
        hold on
        plot(radiusMatrixDeg, meanList(2, :), 'DisplayName', 'Mean')
        yline(meanList(2, 1), '-g', 'HandleVisibility', 'off') % Initial mean value
        if (meanList(2, 1) + selectOptions.mean_thr) <= 1    
            yline(meanList(2, 1) + selectOptions.mean_thr, '--g', 'HandleVisibility', 'off') % Upper Bound
        else
            yline(1, '--g', 'HandleVisibility', 'off')
        end
        if (meanList(2, 1) - selectOptions.mean_thr) >= 0     
            yline(meanList(2, 1) - selectOptions.mean_thr, '--g',  'HandleVisibility', 'off') % Lower Bound
        else
            yline(0, '--g', 'HandleVisibility', 'off') 
        end
        plot(radiusMatrixDeg, stdList(2, :), 'DisplayName', 'Std')
        yline(selectOptions.std_thr, '--r', 'HandleVisibility', 'off') 
        title("Saturation")
        legend
        xticks(radiusMatrixDeg)
        xlabel("Radius (deg)")
        ylim([0, 1])
    subplot(2, 2, 4)
        hold on
        plot(radiusMatrixDeg, meanList(3, :), 'DisplayName', 'Mean')
        yline(meanList(3, 1), '-g', 'HandleVisibility', 'off')
        if (meanList(3, 1) + selectOptions.mean_thr) <= 1    
            yline(meanList(3, 1) + selectOptions.mean_thr, '--g', 'HandleVisibility', 'off') % Upper Bound
        else
            yline(1, '--g', 'HandleVisibility', 'off')
        end
        if (meanList(3, 1) - selectOptions.mean_thr) >= 0     
            yline(meanList(3, 1) - selectOptions.mean_thr, '--g',  'HandleVisibility', 'off') % Lower Bound
        else
            yline(0, '--g', 'HandleVisibility', 'off') 
        end
        plot(radiusMatrixDeg, stdList(3, :), 'DisplayName', 'Std')
        yline(selectOptions.std_thr, '--r', 'HandleVisibility', 'off') 
        legend
        title("Value/Brightness")
        xticks(radiusMatrixDeg)
        xlabel("Radius (deg)")
        ylim([0, 1])
end