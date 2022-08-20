% This program generates grating stimuli and plots predicted gamma for the
% different grating hue/grating combinations used in the paper

% Inputs
% subjectName: 'alpaH' or 'kesariH'
% stimParamsToDisplay = 'SFOR','SZOR' or 'CNOR' for Gratings; 'Hue','SZHue', 'Sat' or 'Val' for Color stimuli

function plotGammaTuningCurves(subjectName,stimParamsToDisplay,displayStimuliFlag)

if ~exist('displayStimuliFlag','var');  displayStimuliFlag=1;           end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make stimParamsList %%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(stimParamsToDisplay(end-1:end),'OR')
    stimType = 'Grating';
else
    stimType = 'HuePatch';
end

% Fixed parameters. Not needed to get predicted gamma but needed to generate the stimuli
gaborStim.azimuthDeg=0;
gaborStim.elevationDeg=0;
gaborStim.sigmaDeg=100000; % The program makeGaborStimulus actually produces Gabors. However, when sigma is extremely large, it is essentially a grating
    
if strcmp(stimType,'Grating')
    var1List = 0:22.5:157.5; % One trace for each orientation, which is the first variable    
else
    if strcmp(stimParamsToDisplay,'Hue')
        var1List = 1; % redundant variable, since only hue trace is shown.
    else
        var1List = 0:60:300; % One trace for each hue, which is the first variable
    end
end
numVar1 = length(var1List);
colorNameList = jet(numVar1);
    
if strcmp(stimParamsToDisplay,'SFOR')
    var2List = [0.5 1 2 4 8]; % Spatial Frequency
    numVar2 = length(var2List);
    gaborStim.radiusDeg=100; % FS
    gaborStim.contrastPC=100;
    
    stimParamsList = cell(numVar1,numVar2);
    predictedGammaList = zeros(numVar1,numVar2);
    for i=1:numVar1
        for j=1:numVar2
            tmpGaborStim = gaborStim;            
            tmpGaborStim.spatialFreqCPD=var2List(j);
            tmpGaborStim.orientationDeg=var1List(i);
            
            stimParamsList{i,j} = tmpGaborStim;
            predictedGammaList(i,j) = getPredictedGamma(subjectName,stimParamsList{i,j});
        end
    end
    
elseif strcmp(stimParamsToDisplay,'SZOR') % Size-Ori
    var2List = [0.3 0.6 1.2 2.4 4.8 9.6]; % Sizes
    numVar2 = length(var2List);
    gaborStim.spatialFreqCPD= 2;
    gaborStim.contrastPC=100;
    
    stimParamsList = cell(numVar1,numVar2);
    predictedGammaList = zeros(numVar1,numVar2);
    for i=1:numVar1
        for j=1:numVar2
            tmpGaborStim = gaborStim;            
            tmpGaborStim.radiusDeg=var2List(j);
            tmpGaborStim.orientationDeg=var1List(i);
            
            stimParamsList{i,j} = tmpGaborStim;
            predictedGammaList(i,j) = getPredictedGamma(subjectName,stimParamsList{i,j});
        end
    end
elseif strcmp(stimParamsToDisplay,'CNOR') % Con-Ori
    var2List = [0 3.1 6.25 12.5 25 50 100]; % Cons
    numVar2 = length(var2List);
    gaborStim.spatialFreqCPD= 2; % 
    gaborStim.radiusDeg=100; % FS
    
    stimParamsList = cell(numVar1,numVar2);
    predictedGammaList = zeros(numVar1,numVar2);
    for i=1:numVar1
        for j=1:numVar2
            tmpGaborStim = gaborStim;            
            tmpGaborStim.contrastPC=var2List(j);
            tmpGaborStim.orientationDeg=var1List(i);
            
            stimParamsList{i,j} = tmpGaborStim;
            predictedGammaList(i,j) = getPredictedGamma(subjectName,stimParamsList{i,j});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGamma = subplot('Position',[0.05 0.05 0.4 0.9]); hold(hGamma,'on');
for i=1:numVar1
    plot(hGamma,predictedGammaList(i,:),'color',colorNameList(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%% Display Stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayStimuliFlag
    hStimuli = getPlotHandles(numVar1,numVar2,[0.5 0.05 0.45 0.9]);
    
    % Each stimulus must be of the size specified by monitorSpecifications.
    [xAxisDeg,yAxisDeg] = getMonitorDetails;
    colormap gray
    for i=1:numVar1
        for j=1:numVar2
            gaborStim = stimParamsList{i,j};
            gaborPatch = makeGaborStimulus(gaborStim,xAxisDeg,yAxisDeg,0);
            pcolor(hStimuli(i,j),xAxisDeg,yAxisDeg,gaborPatch);
            shading(hStimuli(i,j),'interp'); caxis(hStimuli(i,j),[0 100]);
        end
    end
end