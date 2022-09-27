% getMeanEnergy – code returns PSDs and power in specified bands for all
% electrodes that are useful (highRMSElectrodes) and also saves them

% Inputs: subjectName, expDate, protocolName: together these specify the dataset 
%         folderSourceString: specifies the location: default - parent directory
%         freqRangeList: cell array of frequency ranges. Default: {[30 80]}
%         tBL,tST - time range for baseline, and stimulus

% outputs: powerST,powerBL: power in each specified frequency range, for all electrodes and stimulus conditions.
%          electrodeList: list of electrodes
%          rfStats: rf Centers for the electrodes  
%          psdST,psdBL: PSDs
%          freqVals: frequency axis

% STK 240822
% Modified by SR 260922 - simplified so that it only works with the images
% dataset now, and only saves mean energy

function [powerST,powerBL,electrodeList,rfStats,psdST,psdBL,freqVals] = getMeanEnergy(subjectName,expDate,protocolName,folderSourceString,freqRangeList,tBL,tST)

gridType = 'Microelectrode';
if ~exist('folderSourceString','var');  folderSourceString = '';        end
if ~exist('freqRangeList','var'),       freqRangeList = [];             end
if ~exist('tBL','var'),                 tBL = [-0.25 0];                end
if ~exist('tST','var'),                 tST = [0.25 0.5];               end

if isempty(folderSourceString)
    folderSourceString = fileparts(pwd);
end
if isempty(freqRangeList)
    freqRangeList{1} = [30 80];
end

numTapers = 3; % used for Multitaper by Chronux
savedDataDir = fullfile('savedData','spectra');  % local folder where saved power response files are kept
if ~exist(savedDataDir,'dir'),          mkdir(savedDataDir);            end

% get spectral data for all electrodes
rfData = load(fullfile(folderSourceString,'data','rfData',subjectName,[subjectName gridType 'RFData.mat']));
electrodeList = rfData.highRMSElectrodes;
rfStats = rfData.rfStats(electrodeList);

[psdST,psdBL,freqVals] = getPSDMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,tBL,tST,numTapers,1,savedDataDir);

% For each frequency range, get appropriate frequency indices
goodFreqPosList = getGoodFreqPos(freqRangeList,freqVals);
numFreqRanges = length(freqRangeList);
numElectrodes = length(electrodeList);
numImages = size(psdST,2);

% Compute power in frequency bands
powerST = zeros(numFreqRanges,numElectrodes,numImages);
powerBL = zeros(numFreqRanges,numElectrodes,numImages);

for i=1:numFreqRanges
    powerST(i,:,:) = sum(psdST(:,:,goodFreqPosList{i}),3);
    powerBL(i,:,:) = sum(psdBL(:,:,goodFreqPosList{i}),3);
end
end

function badFreqPos = getBadFreqPos(freqVals,deltaF)
badFreqs = 50:50:max(freqVals);
if nargin<2; deltaF = 2; end

badFreqPos = [];
for i=1:length(badFreqs)
    badFreqPos = cat(2,badFreqPos,intersect(find(freqVals>=badFreqs(i)-deltaF),find(freqVals<=badFreqs(i)+deltaF)));
end
end
function goodFreqPosList = getGoodFreqPos(freqRangeList,freqVals)
numFreqRanges = length(freqRangeList);
badFreqPos = getBadFreqPos(freqVals);
goodFreqPosList = cell(1,numFreqRanges);
for i=1:numFreqRanges
    goodFreqPosList{i} = setdiff(find((freqVals>=freqRangeList{i}(1) & freqVals<=freqRangeList{i}(2))),badFreqPos);
end
end
function [psdST,psdBL,freqVals] = getPSDMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,timeRangeBL,timeRangeST,numTapers,saveFlag,folderSave)

numElectrodes = length(electrodeList);
if numTapers>1
    tapers = [ceil((numTapers+1)/2) numTapers];
else
    tapers = [1 1];
end

fileSave = fullfile(folderSave,[subjectName expDate protocolName gridType '_N' num2str(numElectrodes) '_st_' num2str(1000*timeRangeST(1)) '_' ...
    num2str(1000*timeRangeST(2)) '_bl_' num2str(1000*timeRangeBL(1)) '_' num2str(1000*timeRangeBL(2)) '_tapers_' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']);

if exist(fileSave,'file')
    disp(['Loading ' fileSave]);
    x=load(fileSave);
    if isequal(x.electrodeList,electrodeList)
        psdST = x.psdST;
        psdBL = x.psdBL;
        freqVals= x.freqVals;
    else
        error('electrodeList does not match');
    end
else
    disp(['Generating ' fileSave]);
    
    % Get bad trials
    tmp = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
    badTrials = tmp.badTrials;
    
    % Get timeVals
    tmp = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat'));
    timeVals = tmp.timeVals;
    
    % Multi-taper
    Fs  =  round(1/(timeVals(2)-timeVals(1)));
    fMax= 250; 
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 fMax];
    params.trialave = 1;
    
    stDur = diff(timeRangeST);
    blDur = diff(timeRangeBL);
    
    if stDur ~= blDur
        error('baseline and stimulus durations must be equal');
    else
        durS = stDur;
    end
    stPos  = find(timeVals>=timeRangeST(1),1) + (0:round(Fs*durS)-1);
    blPos  = find(timeVals>=timeRangeBL(1),1) + (0:round(Fs*durS)-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Good StimPos %%%%%%%%%%%%%%%%%%%%%%%%%%
    clear goodStimPos
    parameterCombinations = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat')); % ParameterCombinations
    
    numImages   = length(parameterCombinations.fValsUnique);    % numImages
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numFreqs = ceil(fMax*durS);
    psdST = zeros(numElectrodes,numImages,numFreqs);
    psdBL = zeros(numElectrodes,numImages,numFreqs);
    
    for i=1:numElectrodes
        electrodeNum = electrodeList(i);
        
        % Get LFP Data
        tmp = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(electrodeNum) '.mat']));
        analogData = tmp.analogData;
        
        for j=1:numImages    
            clear goodPos
            goodPos = setdiff(parameterCombinations.parameterCombinations{1,1,1,j},badTrials);
            psdST(i,j,:)            = mtspectrumc(analogData(goodPos,stPos)',params);
            [psdBL(i,j,:),freqVals] = mtspectrumc(analogData(goodPos,blPos)',params);           
        end
    end
    if saveFlag
        save(fileSave,'psdST','psdBL','freqVals','electrodeList');
    end
end
end