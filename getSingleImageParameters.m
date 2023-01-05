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
% selectOptions.threshold_method 1 or 2 for single or dual thresholding
% explanation of dual thresholding
% we increase the r from the first radius until the first threshold is crossed
%If the first threshold is crossed for any variable, we check if it reverts
%back to the pre-first threshold value in stop_max iterations. if not, OR
%if another of HSV crosses threshold, then we keep the last sub-threshold radius

%displayAnalysisFlag: 0 if no plot is desired, 1 if 1 plot is desired, 2 if
%two plots are desired

% precomputed indices to pass the 1D NUMERICAL indices , where the 2D mask 

% values are 1. This offers a significant speedup since the masks (which 
% are the same for different images) need not be recomputed every time

% typical usage to obtain HSVR:  getSingleImageParameters_2(imageHSV, imageAxesDeg, rfCenterDeg, 0.3:0.3:9.6, [], 0,patch)  ;
% where patch would have been computed as function [patch,~]=getMaskIndices(0.3:0.3:9.6,rfCenterDeg,xAxisDeg,yAxisDeg,'diff')
% or alernatively patch can be left empty, and will be computed in the function

% Usage to Check particular patch approximations.
% getSingleImageParameters_2(imagehsv, [], [2.0,-3]).
function stimParams = getSingleImageParameters_2(imageHSV,imageAxesDeg,rfCenterDeg,radiusMatrixDeg,selectOptions,displayAnalysisFlag,precomputed_indices)

%set default Values%
if ~exist('imageAxesDeg','var');        imageAxesDeg=[];                end
if ~exist('rfCenterDeg','var');         rfCenterDeg = [0 0];            end
if ~exist('radiusMatrixDeg','var');     radiusMatrixDeg = 0.3:0.3:9.6;  end
if ~exist('displayAnalysisFlag','var'); displayAnalysisFlag=1;          end    
if ~exist('selectOptions','var');       selectOptions=[];               end
if ~exist('precomputed_indices','var'); precomputed_indices={};         end
numRadii = length(radiusMatrixDeg);


%%%%%%%%%%%%%%%%%%%%%%%%% Get selection thresholds %%%%%%%%%%%%%%%%%%%%%%%%
if isempty(selectOptions)
    selectOptions.mean_thr = 5e-2;    
    selectOptions.std_thr = selectOptions.mean_thr*2;
    selectOptions.mean_thr2 = selectOptions.mean_thr*3;% for dual thresholding
    selectOptions.std_thr2 = selectOptions.std_thr*2;
    selectOptions.stop_max=10; % maximum allowed iterations for which the the values can be above 1st threshold and below 2nd
    selectOptions.method = 'vector'; 
    selectOptions.measure = 'diff';
    selectOptions.threshold_method=2; % 2 is dual thresholding, 1 is single thresholding
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

%%%%% Get GoodPosition indices if  precomputed indices are not passed %%%%%
if isempty(precomputed_indices)
   [precomputed_indices,largest_aperture]=getMaskIndices(radiusMatrixDeg,rfCenterDeg,xAxisDeg,yAxisDeg,selectOptions.measure);
end

%%%%%%%%%%%%%%%%%% Check if the image is achromatic %%%%%%%%%%%%%%%%%%%%%%%
achromaticFlag=0;
if achromaticFlag
    return
end
%%%%%%%%%%%%%%%%%%%%%%%% Variable  initializations %%%%%%%%%%%%%%%%%%%%%%%%
if displayAnalysisFlag
    meanList = zeros(3,numRadii);
    stdList = zeros(3,numRadii);
end
r_idx=1;
stop_flag=0; % indicates if thresholds have been crossed;should stay within  [0, stop_max]                
crossed_variable=-1; % keeps track of which variable (H/S/V) is beyond the threhold
goodPosToUse=precomputed_indices{1};
[centralMean,~,~]=getPatch_MeanStd(imageHSV,goodPosToUse,zeros(3,1));

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Computation %%%%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:numRadii    
    
    goodPosToUse =precomputed_indices{i}; 
    [mean_val,std_val,HSV_diff]=getPatch_MeanStd(imageHSV,goodPosToUse,centralMean);

    if displayAnalysisFlag % we need to store the values
        meanList(:,i)=mean_val;
        stdList(:,i)=std_val;
    end

    if(stop_flag<0 || stop_flag>=selectOptions.stop_max)
        if displayAnalysisFlag continue    
        else break;
        end
    end

    for j=1:3       
    if( HSV_diff(j) > selectOptions.mean_thr | std_val(j) > selectOptions.std_thr)% first threshold has been exceeded 
        if (selectOptions.threshold_method==1) % single thresholding
            r_idx = i-1; 
            stop_flag=-1; break;
        end
        % dual thresholding
        if(crossed_variable==-1)
            % no variable has previosuly crosssed threshold, so we note
            % down the current radius and variable crossing threshold
            r_idx = i-1; 
            crossed_variable=j; 
            stop_flag=1;
        end
        if(HSV_diff > selectOptions.mean_thr2 | std_val(j) > selectOptions.std_thr2)
            % second threshold has been exceeded      
            stop_flag=-1;break; 

        elseif (crossed_variable~=j && crossed_variable~=-1)
            % If another variable had exceeded the first threshold in the
            % previous iteration, it means that >=2 variables are around
            % threshold levels, so we stop .
            stop_flag=-1;break;   

        elseif (crossed_variable==j) 
            % this variable has crossed threshold previosuly
            stop_flag=stop_flag+1; 

        elseif(crossed_variable==-1)
            % no variable has previosuly crosssed threshold, so we note
            % down the current radius and variable crossing threshold
            r_idx = i-1; 
            crossed_variable=j; 
        end 
    end
    end
    % if some variale had previosuly crossed 1st threshold, but all of them
    % have currently returned back to sub-1st threshold-levels.
    if(crossed_variable~=-1 & HSV_diff <= selectOptions.mean_thr & std_val <= selectOptions.std_thr )
        stop_flag=0; 
        crossed_variable=-1;
    end        
    
    if(i==numRadii)
            r_idx=i;
    end
end


if r_idx<1 % even patches which are highly variable in the beginning to be chosen
    r_idx=1;
end
if r_idx == 0
    error("r = 0, no suitable patch approx. possible. Consider increasing std_thr."); 
    % While computing r for meanList, it is impossible to go below 1. This is because deviation is considered w.r.t. to the 1st value, and can only exceed the threshold as early as r = 2. For std, as early as r = 1 we can go above std_thr.
    stimParams.radiusDeg = nan;    stimParams.hueDeg = nan;    stimParams.saturation = nan;    stimParams.value = nan;
else
    stimParams.radiusDeg = radiusMatrixDeg(r_idx);
    goodPosToUse=getMaskIndices(stimParams.radiusDeg,rfCenterDeg,xAxisDeg,yAxisDeg,'abs');
    [stimParams.hueDeg, stimParams.saturation,stimParams.value ] =computePatchAvgHSV(imageHSV,goodPosToUse{1},selectOptions.method);
    stimParams.hueDeg=stimParams.hueDeg*360;
end

if displayAnalysisFlag==1
    Display1(imageHSV,xAxisDeg,yAxisDeg,largest_aperture,rfCenterDeg,radiusMatrixDeg,r_idx,selectOptions,meanList,stdList,stimParams)
end
if displayAnalysisFlag==2
    Display2(imageHSV,xAxisDeg,yAxisDeg)
    Display1(imageHSV,xAxisDeg,yAxisDeg,largest_aperture,rfCenterDeg,radiusMatrixDeg,r_idx,selectOptions,meanList,stdList,stimParams)
end
end

% Function that return indices of a concentric circular masks
% Inputs the list of radii, the center, and the pixel coordinates in V.A
function [precomputed_indices,aperture]=getMaskIndices(radiusMatrixDeg,rfCenterDeg,xAxisDeg,yAxisDeg,measure)
    numRadii = length(radiusMatrixDeg);
    precomputed_indices=cell(numRadii,1);
    %setting the defaults for calling makeGaborStimulus()
    gaborStim.azimuthDeg=rfCenterDeg(1);
    gaborStim.elevationDeg=-rfCenterDeg(2);
    gaborStim.spatialFreqCPD=0; % These do not matter since we are only interested in the aperture
    gaborStim.sigmaDeg=100000;
    gaborStim.orientationDeg=0;
    gaborStim.contrastPC = 100;

    goodPosPreviousRadius=[]; % Keeps the goodPos of the previous radius to perform the 'diff' computation

    for i=1:numRadii        
        gaborStim.radiusDeg=radiusMatrixDeg(i);
        [~,aperture] = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg); % generate a circular mask      
        goodPos = find(aperture==1);% Get goodPos, the indices where mask==1

        if strcmp(measure,'diff') 
            goodPosToUse = setdiff(goodPos,goodPosPreviousRadius);% find the indices of the annular disk between the previous radius and current one
            goodPosPreviousRadius=goodPos; % update previous radis goodpos

        elseif strcmp(measure,'abs')% if  measure is not diff 
            goodPosToUse = goodPos;              
        end
        precomputed_indices(i)={goodPosToUse};        
    end
end

%Function that retruns means, std and HSV Diff of HSV image.
%Inputs: hsv Image, 1D mask indices, and the mean of the center of patch
function [mean_val,std_val,HSV_diff]=getPatch_MeanStd(imageHSV,goodPosToUse,centralMean)
    mean_val=zeros(3,1);
    std_val=zeros(3,1);
    HSV_diff=zeros(3,1);
    for j=1:3
        tmp = squeeze(imageHSV(:,:,j));
        goodVals = tmp(goodPosToUse);
        % we compute the average and HSV along with Standard deviations
        if j==1 % Hue
            mean_val(j) = circ_mean(goodVals*2*pi)/(2*pi); % lies in [-0.5, 0.5] 
            std_val(j) = circ_std(goodVals*2*pi)/(2*pi); 
            HSV_diff(j)=abs(mean_val(j)-centralMean(j));
            % if HSV diff is >0.5 , it means we've taken the obtuse
            % difference instead of the acute difference, which we correct
            if HSV_diff(j)>0.5
                HSV_diff(j)=1-HSV_diff(j); 
            end
            HSV_diff(j)=HSV_diff(j)*2;

        else % saturation and value
            mean_val(j) = mean(goodVals); 
            std_val(j) = std(goodVals);
            HSV_diff(j)=abs(mean_val(j)-centralMean(j));
            
        end
    end
end

% function that returns the average HSV values in a Patch
% Inputs HSV image, 1D indices of mask, and the method 'vector' or 'naive'
function [h,s,v]=computePatchAvgHSV(imageHSV,goodPos,method)
    H = squeeze(imageHSV(:, :, 1));
    H = H(goodPos);
    S = squeeze(imageHSV(:, :, 2));
    S = S(goodPos);
    V = squeeze(imageHSV(:, :, 3));
    V = V(goodPos);
    
    if method == "vector" % Vector averaging for H, S as expected in a cylindrical basis  
        N = numel(H);
        X_tot = sum((S.*cos(2*pi*H)));
        Y_tot = sum((S.*sin(2*pi*H)));

        h = (2*pi*(atan2(Y_tot, X_tot) < 0) + atan2(Y_tot, X_tot))/(2*pi); % Radians [-pi, pi]. This function takes care of signs so we don't have to
        s = (((X_tot)^2 + (Y_tot)^2)^0.5)/N;

    elseif method == "naive"
        h = abs(circ_mean(2*pi*H)/(2*pi));
     
        s = mean(S);
    end
    v = mean(V);
    
end

% code to display HSV channels indipendently
function Display2(imageHSV,xAxisDeg,yAxisDeg)

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
end

%code to display image, average values,std in HSV and concentric boundaries 
function Display1(imageHSV,xAxisDeg,yAxisDeg,aperture,rfCenterDeg,radiusMatrixDeg,r_idx,selectOptions,meanList,stdList,stimParams)
    numRadii=length(radiusMatrixDeg);
    h=stimParams.hueDeg/360;  s=stimParams.saturation;  v=stimParams.value;

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
        %imshow(flipud(patch_approx), 'XData', xAxisDeg, 'YData', yAxisDeg)
        set(gca,'YDir', 'normal')
        alpha(0.7)
        % Showing concentric circles for all radii:
        params(1) = rfCenterDeg(1);
        params(2) = rfCenterDeg(2);
        params(5) = 0;
        params(6) = 1;
        for i = 1:numRadii
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
        polarplot(meanList(1, 1)*ones(numRadii), radiusMatrixDeg, '-g', 'HandleVisibility', 'off') % Initial mean value
        polarplot((meanList(1, 1) + selectOptions.mean_thr*2*pi)*ones(numRadii), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Upper Bound
        polarplot((meanList(1, 1) - selectOptions.mean_thr*2*pi)*ones(numRadii), radiusMatrixDeg, '--g', 'HandleVisibility', 'off') % Lower Bound
        polarplot(stdList(1, :), radiusMatrixDeg, 'DisplayName', 'Std')
        polarplot((selectOptions.std_thr*2*pi)*ones(numRadii), radiusMatrixDeg, '--r', 'HandleVisibility', 'off') % Upper Bound on Std
        legend
        title("Hue (deg)")
        % For ease of viewing
        if numRadii > 10 
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

