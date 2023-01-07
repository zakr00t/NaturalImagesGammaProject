% runGetSingleImageParameters

%%%%%%%%%%%%%%% Make some dummy images and see how it works %%%%%%%%%%%%%%%

% Color patch
gaborStim.azimuthDeg=0;
gaborStim.elevationDeg=0;

gaborStim.sigmaDeg=100000;
gaborStim.orientationDeg=0;
gaborStim.contrastPC = 25;

gaborStim.spatialFreqCPD=0;
gaborStim.spatialFreqPhaseDeg=90; % For color patches, SF=0 and SFphase=90

gaborStim.radiusDeg=5;
gaborStim.hueDeg=45;
gaborStim.sat = 0.5;

[xAxisDeg,yAxisDeg] = getMonitorDetails;
imageRGB = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg,0);
imageAxesDeg.xAxisDeg=xAxisDeg;
imageAxesDeg.yAxisDeg=yAxisDeg;
rfCenterDeg=[gaborStim.azimuthDeg gaborStim.elevationDeg];

params = getSingleImageParameters(rgb2hsv(imageRGB),imageAxesDeg,rfCenterDeg);