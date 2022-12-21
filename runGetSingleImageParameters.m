% runGetSingleImageParameters

%%%%%%%%%%%%%%% Make some dummy images and see how it works %%%%%%%%%%%%%%%

% Color patch
gaborStim.azimuthDeg=0;
gaborStim.elevationDeg=0;
gaborStim.spatialFreqCPD=1;
gaborStim.sigmaDeg=100000;
gaborStim.orientationDeg=0;
gaborStim.contrastPC = 100;
gaborStim.spatialFreqPhaseDeg=0; % For color patches, SF=0 and SFphase=90

gaborStim.radiusDeg=5;
gaborStim.hueDeg=45;
gaborStim.sat = 0.5;

[xAxisDeg,yAxisDeg] = getMonitorDetails;
imageRGB = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg,0);
imageAxesDeg.xAxisDeg=xAxisDeg;
imageAxesDeg.yAxisDeg=yAxisDeg;
rfCenterDeg=[gaborStim.azimuthDeg gaborStim.elevationDeg];

getSingleImageParameters(rgb2hsv(imageRGB),imageAxesDeg,rfCenterDeg);