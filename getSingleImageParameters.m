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
% selectOptions.meanThr: The % +/- deviation from the initial H, S, and V 
% values to be tolerated. Default: 5e-2 (5%). 1x3 array of values.
% selectOptions.stdThr: Similar purpose as above, but non-negative.
% Default: 0.1 (10%, 2*meanThr). 1x3 array of values.
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
    selectOptions.meanThr = 5e-2*ones(1, 3);
    selectOptions.stdThr = 2*selectOptions.meanThr;
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

numStimuli = length(radiusMatrixDeg);
meanList = zeros(3,numStimuli);
stdList = zeros(3,numStimuli);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the apertures %%%%%%%%%%%%%%%%%%%%%%%%%%
gaborStim.azimuthDeg=rfCenterDeg(1);
gaborStim.elevationDeg=rfCenterDeg(2);
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
            stdList(j,i) = real(circ_std(goodVals*2*pi)); % circ_std can give a complex number out if all values are the same
        else
            meanList(j,i) = mean(goodVals);
            stdList(j,i) = std(goodVals);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Parameter Computation %%%%%%%%%%%%%%%%%%%%%%%% 
% The radius before the first radius at which a threshold is crossed is chosen:
rList = [];
for j = 1:3
    foo = abs(meanList(j, :) - meanList(j, 1));
    if j == 1 % Dealing with hues in radians
        foo(foo > pi) = 2*pi - foo(foo > pi); % Ensures that reflex (>180 degree) angles are replaced with their smaller counterpart. E.g., if two vectors are separated by 270deg, it is more appropriate to say they are separated by 90deg.
        foo = foo/(2*pi); % bring it back to 0-1 for this computation
    end
    if ~isempty(find(foo > selectOptions.meanThr(j), 1 ))
        rList = cat(2,rList, find(foo > selectOptions.meanThr(j), 1 ));
    else
        rList = cat(2,rList, numStimuli);
    end
    if ~isempty(find(stdList(j, :) > selectOptions.stdThr(j), 1 ))
        rList = cat(2,rList, find(stdList(j, :) > selectOptions.stdThr(j), 1 ));
    else
        rList = cat(2,rList, numStimuli);
    end
end
rIdx = min(rList); % The index of the maximum approximation radius which doesn't deviate beyond threshold

gaborStim.radiusDeg = radiusMatrixDeg(rIdx);
[~, aperture] = makeGaborStimulus(gaborStim, xAxisDeg, yAxisDeg);
   
% Isolating HSV channel values (0-1 range) within aperture
H = squeeze(imageHSV(:, :, 1));
H = H(aperture == 1);
S = squeeze(imageHSV(:, :, 2));
S = S(aperture == 1);
V = squeeze(imageHSV(:, :, 3));
V = V(aperture == 1);

if selectOptions.method == "vector" % Vector averaging for H, S as expected in a cylindrical basis  
    N = numel(H);
    xTot = sum((S.*cos(2*pi*H)));
    yTot = sum((S.*sin(2*pi*H)));
    h = (2*pi*(atan2(yTot, xTot) < 0) + atan2(yTot, xTot))/(2*pi); % atan2 has output in radians [-pi, pi], which we convert to [0, 2pi]. This function takes care of signs so we don't have to
    s = (((xTot)^2 + (yTot)^2)^0.5)/N;
elseif selectOptions.method == "naive"
    if circ_mean(2*pi*H) >= 0
        h = circ_mean(2*pi*H)/(2*pi);
    else
        h = 1 + circ_mean(2*pi*H)/(2*pi);
    end
    s = mean(S);
end
v = mean(V);

stimParams.radiusDeg = radiusMatrixDeg(rIdx);
stimParams.hueDeg = h*360;
stimParams.sat = s;

% Conversion from value (v) to contrastPC
% Since we are approximating a patch by gabor stimulus, contrastPC=0 means
% value of 0.5. Higher values can be achieved by increasing the contrast
% and making the spatial phase = 90 degrees. Importantly, to get values lower than
% 0.5, contrastPC still needs to be increased, but now spatial phase will
% be -90 degrees. val of 0 corresponds to contrastPC=100 and phi = -90
% degrees, while val of 1 corresponds to contrastPC=100 and phi = +90
stimParams.contrastPC = 2*abs(v-0.5)*100;
stimParams.spatialFreqPhaseDeg = 90*sign(v-0.5);

if displayAnalysisFlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    sgtitle("Image HSV Channels");
    % Image:
    subplot(2, 2, 1);
        % Steps necessary to keep the image upright and with an ascending Y axis: (REFER https://in.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab)
        imshow(hsv2rgb(imageHSV), 'XData', xAxisDeg, 'YData', yAxisDeg);
        set(gca,'YDir', 'normal');
        axis on;
        title("RGB Image");
        xlabel("Azimuth (deg)");
        ylabel("Elevation (deg)");
    % Hue Channel:
    subplot(2, 2, 2);
        imshow(imageHSV(:, :, 1), 'XData', xAxisDeg, 'YData', yAxisDeg);
        set(gca,'YDir', 'normal');
        axis on;
        caxis([0, 1]);
        colormap('hsv');
        colorbar;
        title("Hue");
        xlabel("Azimuth (deg)");
        ylabel("Elevation (deg)");
    % Saturation Channel:
    subplot(2, 2, 3);
        imshow(imageHSV(:, :, 2), 'XData', xAxisDeg, 'YData', yAxisDeg);
        set(gca,'YDir', 'normal');
        axis on;
        caxis([0, 1]);
        colorbar;
        title("Saturation");
        xlabel("Azimuth (deg)");
        ylabel("Elevation (deg)");
    % Value Channel
    subplot(2, 2, 4);
        imshow(imageHSV(:, :, 3), 'XData', xAxisDeg, 'YData', yAxisDeg);
        set(gca,'YDir', 'normal');
        axis on;
        caxis([0, 1]);
        colorbar;
        title("Value/Brightness");
        xlabel("Azimuth (deg)");
        ylabel("Elevation (deg)");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    sgtitle("Image Grating/Patch Approximation")
    % Image:
    subplot(2, 2, 1)
        imshow(hsv2rgb(imageHSV), 'XData', xAxisDeg, 'YData', yAxisDeg)
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
        patchApprox = hsv2rgb(cat(3, chan1, chan2, chan3));
        hold on;
        imshow(patchApprox, 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        alpha(0.7)
%         % Showing concentric circles for all radii:
%         params(1) = rfCenterDeg(1);
%         params(2) = rfCenterDeg(2);
%         params(5) = 0;
%         params(6) = 1;
%         for i = 1:numStimuli
%             params(3) = radiusMatrixDeg(i);
%             params(4) = params(3);
%             [~, ~, boundaryX, boundaryY] = gauss2D(params, xAxisDeg, yAxisDeg, []);
%             if i == rIdx
%                 plot(boundaryX, boundaryY, '-g', 'DisplayName', sprintf("r = %0.2g", radiusMatrixDeg(rIdx))) % Highlight the patch radius
%             else
%                 plot(boundaryX, boundaryY, '-k', 'HandleVisibility', 'off')
%             end
%         end
        axis on
%        legend
   subplot(2, 2, 2)
        polarplot(meanList(1, :), radiusMatrixDeg, 'DisplayName', 'Mean')
        hold on
        polarplot(meanList(1, 1)*ones(numStimuli), radiusMatrixDeg, '-g', 'HandleVisibility', 'off') % Initial mean value
        polarplot((meanList(1, 1) + selectOptions.meanThr(1)*2*pi)*ones(numStimuli), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Upper Bound
        polarplot((meanList(1, 1) - selectOptions.meanThr(1)*2*pi)*ones(numStimuli), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Lower Bound
        polarplot(stdList(1, :), radiusMatrixDeg, 'DisplayName', 'Std')
        polarplot((selectOptions.stdThr(1)*2*pi)*ones(numStimuli), radiusMatrixDeg, '--r', 'HandleVisibility', 'off') % Upper Bound on Std
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
        if (meanList(2, 1) + selectOptions.meanThr(2)) <= 1    
            yline(meanList(2, 1) + selectOptions.meanThr(2), '--g', 'HandleVisibility', 'off') % Upper Bound
        else
            yline(1, '--g', 'HandleVisibility', 'off')
        end
        if (meanList(2, 1) - selectOptions.meanThr(2)) >= 0     
            yline(meanList(2, 1) - selectOptions.meanThr(2), '--g',  'HandleVisibility', 'off') % Lower Bound
        else
            yline(0, '--g', 'HandleVisibility', 'off')
        end
        plot(radiusMatrixDeg, stdList(2, :), 'DisplayName', 'Std')
        yline(selectOptions.stdThr(2), '--r', 'HandleVisibility', 'off')
        title("Saturation")
        legend
        xticks(radiusMatrixDeg)
        xlabel("Radius (deg)")
        ylim([0, 1])
    subplot(2, 2, 4)
        hold on
        plot(radiusMatrixDeg, meanList(3, :), 'DisplayName', 'Mean')
        yline(meanList(3, 1), '-g', 'HandleVisibility', 'off')
        if (meanList(3, 1) + selectOptions.meanThr(3)) <= 1    
            yline(meanList(3, 1) + selectOptions.meanThr(3), '--g', 'HandleVisibility', 'off') % Upper Bound
        else
            yline(1, '--g', 'HandleVisibility', 'off')
        end
        if (meanList(3, 1) - selectOptions.meanThr(3)) >= 0     
            yline(meanList(3, 1) - selectOptions.meanThr(3), '--g',  'HandleVisibility', 'off') % Lower Bound
        else
            yline(0, '--g', 'HandleVisibility', 'off')
        end
        plot(radiusMatrixDeg, stdList(3, :), 'DisplayName', 'Std')
        yline(selectOptions.stdThr(3), '--r', 'HandleVisibility', 'off')
        legend
        title("Value/Brightness")
        xticks(radiusMatrixDeg)
        xlabel("Radius (deg)")
        ylim([0, 1])
end