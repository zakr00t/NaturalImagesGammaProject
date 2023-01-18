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
% selectOptions.measure = 'diff' or 'abs'. When set to diff ,we compute the
% mean and std of only the pixels in the disk between the current radius
% and the preceding one. When set to 'abs', we take all pixels. Default: 'diff'

function getSingleImageParameters(imageHSV,imageAxesDeg,rfCenterDeg,radiusMatrixDeg,selectOptions,displayAnalysisFlag)

if ~exist('imageAxesDeg','var');        imageAxesDeg=[];                end
if ~exist('rfCenterDeg','var');         rfCenterDeg = [0 0];            end
if ~exist('radiusMatrixDeg','var');     radiusMatrixDeg = 0.3:0.3:9.6;  end
if ~exist('displayAnalysisFlag','var'); displayAnalysisFlag=1;          end    
if ~exist('selectOptions','var');       selectOptions=[];               end

%%%%%%%%%%%%%%%%%%%%%%%%% Get selection thresholds %%%%%%%%%%%%%%%%%%%%%%%%
if isempty(selectOptions)
    selectOptions.measure = 'diff';
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Get the aperatures %%%%%%%%%%%%%%%%%%%%%%%%
    gaborStim.azimuthDeg=rfCenterDeg(1);
    gaborStim.elevationDeg=rfCenterDeg(2);
    gaborStim.spatialFreqCPD=0; % These do not matter since we are only interested in the aperature
    gaborStim.sigmaDeg=100000;
    gaborStim.orientationDeg=0;
    gaborStim.contrastPC = 100;

    goodPosPreviousRadius=false(size(imageHSV,1),size(imageHSV,2)); % Keeps the goodPos of the previous radius to perform the 'diff' computation
    for i=1:numStimuli
        gaborStim.radiusDeg=radiusMatrixDeg(i);
        [~,aperature] = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg);

        % Get goodPos
        goodPos = (aperature==1);
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
                meanList(j,i) = mean(goodVals); % Change from 0-1 to 0-2*pi 
                stdList(j,i) = std(goodVals);
            end
        end
    end
end

if displayAnalysisFlag
    colormap gray
    
    if achromaticFlag
    else
        
        hPlots = getPlotHandles(4,3,[0.05 0.05 0.9 0.9],0.05,0.05);
        
        %%%%%%%%%%%%%%%%%%%% Plot the full image %%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(hPlots(1,1));
        imagesc(xAxisDeg,yAxisDeg,hsv2rgb(imageHSV),'parent',hPlots(1,1));
        colorbar;
        
        % Plot the image and also show the radii for which values are calculated
        params(1) = rfCenterDeg(1);
        params(2) = rfCenterDeg(2);
        params(5) = 0;
        params(6) = 1;
        
        colorNamesList = jet(numStimuli);
        imagesc(xAxisDeg,yAxisDeg,hsv2rgb(imageHSV),'parent',hPlots(1,2));
        hold(hPlots(1,2),'on');
        for i=1:numStimuli
            params(3) = radiusMatrixDeg(i);
            params(4) = params(3);
            [~,~,boundaryX,boundaryY] = gauss2D(params,xAxisDeg,yAxisDeg,[]);
            plot(hPlots(1,2),boundaryX,boundaryY,'color',colorNamesList(i,:));
        end
        
        for i=1:3 % Plot the results separately for H,S and V
            axes(hPlots(i+1,1)); %#ok<LAXES>
            if i==1
                pcolor(xAxisDeg,yAxisDeg,360*imageHSV(:,:,i));
                plot(hPlots(i+1,2),radiusMatrixDeg,meanList(i,:)*180/pi);
            else
                pcolor(xAxisDeg,yAxisDeg,imageHSV(:,:,i));
                plot(hPlots(i+1,2),radiusMatrixDeg,meanList(i,:));
            end
            colorbar;
            shading(hPlots(i+1,1),'interp');
            
            plot(hPlots(i+1,3),radiusMatrixDeg,stdList(i,:));
        end
    end
end
end