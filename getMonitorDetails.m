% This program returns the monitor x and y coordinates in degrees. It also stores the default monitor specifications used for this project

function [xAxisDeg,yAxisDeg,monitorSpecifications,viewingDistanceCM] = getMonitorDetails(monitorSpecifications,viewingDistanceCM)

% Default values
if ~exist('monitorSpecifications','var')
    monitorSpecifications.height = 11.8; % Monitor height in inches
    monitorSpecifications.width  = 20.9; % Monitor width in inches
    monitorSpecifications.xRes   = 1280; % Resolution (num of pixels) of monitor
    monitorSpecifications.yRes   = 720; % Resolution (num of pixels) of monitor
end
if ~exist('viewingDistanceCM','var');       viewingDistanceCM = 50;     end


viewingDistance = viewingDistanceCM/2.54; % convert to inches
yDeg = (atan((monitorSpecifications.height/2)/viewingDistance))*180/pi; % Half the space, in degrees
xDeg = (atan((monitorSpecifications.width/2)/viewingDistance))*180/pi;

imageXRes = monitorSpecifications.xRes;
imageYRes = monitorSpecifications.yRes;
xAxisDeg = -xDeg: (2*xDeg/imageXRes) :xDeg; xAxisDeg = xAxisDeg(1:imageXRes);
yAxisDeg = -yDeg: (2*yDeg/imageYRes) :yDeg; yAxisDeg = yAxisDeg(1:imageYRes);
    
end