% This program returns the predicted gamma response for any achromatic
% grating or hue patch from the tuning functions that are obtained experimentally 

% Inputs
% Required
% subjectName: 'alpaH' or 'kesariH'
% stimParams: description of an achromatic or chromatic grating 

% Optional
% modelParams: model parameters to construct the tuning functions. If not
% defined, we use the default median values that are saved for each monkey
% scalingFactor: the overall gamma is scaled by scalingFactor(1) and then
% an offset of scalingFactor(2) is added. Default: [1 0]

function gammaVal = getPredictedGamma(subjectName,stimParams,modelParams,scalingFactor)

if ~exist('modelParams','var');         modelParams=[];                 end
if ~exist('scalingFactor','var');       scalingFactor=[1 0];            end

if isfield(stimParams,'hueDeg')
    gammaVal = getPredictedGammaHue(subjectName,stimParams,modelParams,scalingFactor);
else
    gammaVal = getPredictedGammaGrating(subjectName,stimParams,modelParams,scalingFactor);
end
end
function val = getPredictedGammaGrating(subjectName,stimParams,modelParams,scalingFactor)

if isempty(modelParams)
    if strcmp(subjectName,'alpaH')
        modelParams.sfCenter = 2.9399;
        modelParams.sfSigma = 1.1339;
        modelParams.oriCenter = 1.9133;
        modelParams.oriSpreadK = 0.51753;
        modelParams.sizeSlope = 0.38268;
        modelParams.sizeAtHalfMax = 1.664;
        modelParams.conSlope = 1.3317;
        modelParams.conAtHalfMax = 2.637;
        
    elseif strcmp(subjectName,'kesariH')
        modelParams.sfCenter = 1.5304;
        modelParams.sfSigma = 1.1500;
        modelParams.oriCenter =  0.5558;   
        modelParams.oriSpreadK = 0.2498;
        modelParams.sizeSlope = 0.0580;
        modelParams.sizeAtHalfMax = 5.0000;
        modelParams.conSlope = 9.6000;
        modelParams.conAtHalfMax = 2.0433;       
    end
end

f = stimParams.spatialFreqCPD;
o = stimParams.orientationDeg;
s = stimParams.radiusDeg;
c = stimParams.contrastPC;

sfVal = normpdf(log2(f),modelParams.sfCenter,modelParams.sfSigma)./normpdf(modelParams.sfCenter,modelParams.sfCenter,modelParams.sfSigma);
oriVal = circ_vmpdf(2*deg2rad(o),2*modelParams.oriCenter,modelParams.oriSpreadK)./circ_vmpdf(2*modelParams.oriCenter,2*modelParams.oriCenter,modelParams.oriSpreadK);
sizeVal = sigmoidalFunction(log2(s),[modelParams.sizeSlope modelParams.sizeAtHalfMax]);
conVal = sigmoidalFunction(log2(c),[modelParams.conSlope modelParams.conAtHalfMax]);

% Predicted gamma
val = scalingFactor(1).*sfVal.*oriVal.*sizeVal.*conVal + scalingFactor(2);
end
function val = getPredictedGammaHue(subjectName,stimParams,modelParams,scalingFactor)

if isempty(modelParams)
    if strcmp(subjectName,'alpaH')
        modelParams.CenterR = 5.8954;
        modelParams.SpreadR = 1.1214;
        modelParams.CenterB = 3.3357;
        modelParams.SpreadB = 5.3724;
        modelParams.CenterG = 2.3736;
        modelParams.SpreadG = 2.6017;
        modelParams.GainBbyR = 0.3831;
        modelParams.GainGbyR = 0.3310;
        modelParams.sizeSlope = 9.5148;
        modelParams.sizeAtHalfMax = 2.5138;
        modelParams.satSlope = 2.6627;
        modelParams.satAtHalfMax = 0.6760;
        
    elseif strcmp(subjectName,'kesariH')
        modelParams.CenterR = 5.9521;
        modelParams.SpreadR = 3.6366;
        modelParams.CenterB = 3.9936;
        modelParams.SpreadB = 8.3353;
        modelParams.CenterG = 1.7514;
        modelParams.SpreadG = 200.00;
        modelParams.GainBbyR = 0.1905;    
        modelParams.GainGbyR = 0.6608;
        modelParams.sizeSlope = 1.3618;
        modelParams.sizeAtHalfMax = 2.1889;
        modelParams.satSlope = 2.8089;
        modelParams.satAtHalfMax = 0.6792;      
    end
end

o = stimParams.hueDeg;
s = stimParams.radiusDeg;
st= stimParams.saturation;

hueVal =                      (circ_vmpdf(deg2rad(o),modelParams.CenterR,modelParams.SpreadR)./circ_vmpdf(modelParams.CenterR,modelParams.CenterR,modelParams.SpreadR)...  % sum of 3 VonMises
       + modelParams.GainBbyR.*circ_vmpdf(deg2rad(o),modelParams.CenterB,modelParams.SpreadB)./circ_vmpdf(modelParams.CenterB,modelParams.CenterB,modelParams.SpreadB)...
       + modelParams.GainGbyR.*circ_vmpdf(deg2rad(o),modelParams.CenterG,modelParams.SpreadG)./circ_vmpdf(modelParams.CenterG,modelParams.CenterG,modelParams.SpreadG)); 
sizeVal = sigmoidalFunction(log2(s),[modelParams.sizeSlope modelParams.sizeAtHalfMax]);
satVal  = sigmoidalFunction(st,[modelParams.satSlope  modelParams.satAtHalfMax]);

% get predicted gamma
val = scalingFactor(1).*hueVal.*sizeVal.*satVal + scalingFactor(2);

end
function y = sigmoidalFunction(x,coeff)
slope = coeff(1); offset = coeff(2);
y = 1./(1+10.^(slope*(x-offset)));
end