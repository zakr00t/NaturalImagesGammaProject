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
% if actual data (a vector of values of same length as number of stimuli) is sent in as scaling factor,
% then this dats is used to get overall gain and offset on scaled predictions 

% Modified by STK. Only hue pipeline is modified as of now.

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
        modelParams.centerR = 5.9925;
        modelParams.spreadR = 1.8473;
        modelParams.centerC = 3.0529;
        modelParams.spreadC = 2.5153;
        modelParams.gainCbyR = 0.4259;
        modelParams.sizeSlope = 1.4496;
        modelParams.sizeAtHalfMax = 0.4119;
        modelParams.satSlope = 3.2922;
        modelParams.satAtHalfMax = 0.7477;
        modelParams.valSlope = 0.0041;
        modelParams.valIntercept = 1.4746;
        
    elseif strcmp(subjectName,'kesariH')
        modelParams.centerR = 6.2428;
        modelParams.spreadR = 4.0667;
        modelParams.centerC = 3.8399;
        modelParams.spreadC = 2.0575;    
        modelParams.gainCbyR = 0.3909;
        modelParams.sizeSlope = 0.9333;
        modelParams.sizeAtHalfMax = 0.2258;
        modelParams.satSlope = 4.2867;
        modelParams.satAtHalfMax = 0.6463; 
        modelParams.valSlope = 1.7862e-05;
        modelParams.valIntercept = 2.7213;     
    end
end

o = stimParams.hueDeg;
s = stimParams.radiusDeg;
st= stimParams.sat;
v = stimParams.contrastPC/100;

hueCoeffs  = [modelParams.centerR modelParams.spreadR modelParams.centerC modelParams.spreadC  modelParams.gainCbyR];
sizeCoeffs = [modelParams.sizeSlope modelParams.sizeAtHalfMax];
satCoeffs  = [modelParams.satSlope modelParams.satAtHalfMax];
valCoeffs  = [modelParams.valSlope modelParams.valIntercept];

hueVal  = sum2VmFunctionsOfRadx(o,hueCoeffs);
sizeVal = sigmoidalFunctionOfLogx(s,sizeCoeffs);
satVal  = sigmoidalFunctionOfx(st,satCoeffs);
valVal  = linearFunctionOfx(v,valCoeffs);

% get predicted gamma - scaled without gain
val = 1.*hueVal.*sizeVal.*satVal.*valVal + 0;  % 

% get gain and intercept if data available
if length(scalingFactor)>2 && length(scalingFactor(:))==length(val) 
    actualData = scalingFactor(:);
    scalingFactor = fitGandO(actualData,val(:)); % get actual scaling factor
    val = scalingFactor(1).*hueVal.*sizeVal.*satVal.*valVal + scalingFactor(2);
end

end
function y = sigmoidalFunction(x,coeff)
slope = coeff(1); offset = coeff(2);
y = 1./(1+10.^(slope*(x-offset)));
end
function [y,sigmoidFunctionStr] = sigmoidalFunctionOfLogx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end  % slope = coeff(1); offset = coeff(2);
sigmoidFunction = @(coeff,x) (1./(1+10.^(coeff(1)*(log2(coeff(2))-log2(x)))));
y = sigmoidFunction(coeff,x);                % y = 1./(1+10.^(slope*(x-offset)));
sigmoidFunctionStr = func2str(sigmoidFunction);
end
function [y,sigmoidFunctionStr] = sigmoidalFunctionOfx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end  % slope = coeff(1); offset = coeff(2);
sigmoidFunction = @(coeff,x) (1./(1+10.^(coeff(1)*(coeff(2)-x))));
y = sigmoidFunction(coeff,x);                % y = 1./(1+10.^-(slope*(x-offset)));
sigmoidFunctionStr = func2str(sigmoidFunction);
end
function [y,sum2VmFunctionStr] = sum2VmFunctionsOfRadx(x,coeff)
if ~exist('x','var'), x = 0; end
if ~exist('coeff','var'), coeff = [0 1 0 1 0]; end
sum2VmFunction =  @(coeff,x)( circ_vmpdf(deg2rad(x),coeff(1),coeff(2))./circ_vmpdf(coeff(1),coeff(1),coeff(2)) + ...  % keep max at 1.
                    coeff(5).*circ_vmpdf(deg2rad(x),coeff(3),coeff(4))./circ_vmpdf(coeff(3),coeff(3),coeff(4)));
y = sum2VmFunction(coeff,x);
sum2VmFunctionStr = func2str(sum2VmFunction);
end
function [y,linearFunctionStr] = linearFunctionOfx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end
linearFunction = @(coeff,x)(coeff(1).*x + coeff(2));
y = linearFunction(coeff,x);
linearFunctionStr = func2str(linearFunction);
end
% optimise gain & offset
function [fitResult] = fitGandO(yData,xData)

tols    = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
X0  =  [max(yData)-min(yData)  min(yData)];  % overall gain & offset
Low =  [0   0 ];
Upp =  [Inf Inf];
funis = @(p,dd) p(1).*dd + p(2);

[fitResult,~,~] = lsqcurvefit(funis,X0,xData,yData,Low,Upp,options);
end