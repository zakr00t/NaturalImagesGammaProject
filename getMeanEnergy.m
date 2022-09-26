% getEnergy – code returns delta power in a given band and power spectra.
% also saves spectral responses of individual protocols

% Inputs: subjectName ,
%         expDate,
%         protocolName,
%         freqRange
%         tBL,tST - time range for baseline, and stimulus 
%         electrodeNum - scalar electrode number or vector of elecs
% outputs: logDelPwrBand: the gamma power change from baseline for that electrode. log10 
%          meanPwrST:    power in time range of ST
%          meanPwrBLall: power in time range of BL
%          freqVals:     frequency axis
%          parameterCombinations: stimulus conditions for the protocol         
% STK 240822 

function [logDelPwrBand,meanPwrST,meanPwrBLall,bandPwrST,bandPwrBL,parameterCombinations]= getMeanEnergy(subjectName,expDate,protocolName,fBandGamma,fBad,tBL,tST,electrodeNum)

gridType            = 'Microelectrode';
folderSourceString  = '/Volumes/Seagate Backup Plus Drive/SIDRAT/';% 'W:\'; % folder where 'data' directory is kept with extracted protocols
if ~exist('fBandGamma','var'),             fBandGamma = [30 80]; end
if ~exist('fBad','var'),                   fBad = [];        end
if ~exist('tBL','var'),                    tBL = [-0.25 0];  end
if ~exist('tST','var'),                    tST = [0.25 0.5]; end
if ~exist('electrodeNum','var'),           electrodeNum = 2; end

numTapers = 3;                                             % used for Multitaper by Chronux
modelDataDir = fullfile('modelData','derivatives','spectra');  % where saved power response files are kept
if ~exist(modelDataDir,'dir'), mkdir(modelDataDir);   end

% get the responses for this particular experiment:
[powerST,powerBL,freqVals,electrodes,parameterCombinations] = getResponses(gridType,folderSourceString,subjectName,expDate,protocolName,tBL,tST,numTapers,modelDataDir);

% select the electrodes:
elPowerST = cell(length(electrodeNum),size(powerST,2),size(powerST,3),size(powerST,4));
elPowerBL = cell(length(electrodeNum),size(powerBL,2),size(powerBL,3),size(powerBL,4));
[~,eind]  = intersect(electrodes,electrodeNum);
elPowerST(:,:,:,:) = powerST(eind,:,:,:,:); 
elPowerBL(:,:,:,:) = powerBL(eind,:,:,:,:);

% get band power and trial averaged spectra
[logDelPwrBand,meanPwrST,meanPwrBLall,bandPwrST,bandPwrBL] = getPowerChanges(elPowerST,elPowerBL,freqVals,{fBandGamma},fBad);
  
end

% % 
function [pwrST,pwrBL,freqVals,electrodeList,parameterCombinations]= getResponses(gridType,folderSourceString,subjectName,expDate,protocolName,tBL,tST,numTapers,saveDir)

if ~exist('numTapers','var')||isempty(numTapers), numTapers = 3; end

% 1. load good electrodes for this monkey
[~,~,LFPElectrodeList,EcogElectrodeList,~] = getRFdetails(subjectName,'modelData'); % modelData dir has 'RFDetails.mat' file
electrodeList=[LFPElectrodeList{1}; EcogElectrodeList{1}]; % unwrap them from cell.

% 2. get responses for protocol. returns in terms of elec x size x sf x con x ori x tf
if ~isempty(electrodeList)
      [pwrST,pwrBL,freqVals,parameterCombinations]  = getPowerMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,tBL,tST,numTapers,1,saveDir);
else   pwrST= nan(5,5,5,5,5); pwrBL= nan(5,5,5,5,5);  % dummy
end
pwrST = squeeze(pwrST); pwrBL = squeeze(pwrBL);

% 3. some updates to parameter combinations
parameterCombinations.sValsUnique(parameterCombinations.sValsUnique<0 | parameterCombinations.sValsUnique>11) = 11.5; % for FullScreen, use same val
[~,indor] = intersect(parameterCombinations.oValsUnique,[22 67 112 157]); % for oris at gap of 22.5 
parameterCombinations.oValsUnique(indor) = parameterCombinations.oValsUnique(indor)+0.5;

end

function [powerST,powerBL,freqVals,parameterCombinations] = getPowerMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,timeRangeBL,timeRangeST,numTapers,saveFlag,folderSave)

numElectrodes = length(electrodeList);
if numTapers>1,  tapers = [ceil((numTapers+1)/2) numTapers]; 
else             tapers = [1 1];
end 

fileSave = fullfile(folderSave,[subjectName expDate protocolName gridType '_N' num2str(numElectrodes) '_st_' num2str(1000*timeRangeST(1)) '_' ...
    num2str(1000*timeRangeST(2)) '_bl_' num2str(1000*timeRangeBL(1)) '_' num2str(1000*timeRangeBL(2)) '_tapers_' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']);

if exist(fileSave,'file')
    disp(['Loading ' fileSave]);
    x=load(fileSave);
    if isequal(x.electrodeList,electrodeList)
        powerST = x.powerST;
        powerBL = x.powerBL;
        freqVals= x.freqVals;
        parameterCombinations= x.parameterCombinations;
    end
else
    disp(['Generating ' fileSave]);
    
    % Get bad trials
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
    % Get timeVals
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat')); 

    % Multi-taper
    Fs  =  round(1/(timeVals(2)-timeVals(1)));
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 250];
    params.trialave = 0;

    stPos  = find(timeVals>=timeRangeST(1),1) + (0:round(Fs*diff(timeRangeST))-1);
    blPos  = find(timeVals>=timeRangeBL(1),1) + (0:round(Fs*diff(timeRangeBL))-1);

     %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Good StimPos %%%%%%%%%%%%%%%%%%%%%%%%%%
    clear goodStimPos
    parameterCombinations = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat')); % ParameterCombinations
    
    aPos=1;ePos=1;          % azimuth, elevation
    numSizes = length(parameterCombinations.sValsUnique);
    numSFs   = length(parameterCombinations.fValsUnique);             % spatial freqs
    numCons  = length(parameterCombinations.cValsUnique);
    numOris  = length(parameterCombinations.oValsUnique);    % orientations
    numTFs   = length(parameterCombinations.tValsUnique);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:numElectrodes
        electrodeNum = electrodeList(i,1);
             % Get LFP Data
        load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(electrodeNum) '.mat']));

        for sPos = 1:numSizes
            for fPos=1:numSFs 
                for cPos = 1:numCons
                    for oPos = 1:numOris
                        for tPos = 1:numTFs
                            clear goodPos
                            goodPos = setdiff(parameterCombinations.parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos},badTrials);
                            [powerST{i,sPos,fPos,cPos,oPos,tPos}]         = mtspectrumc(analogData(goodPos,stPos)',params);
                            [powerBL{i,sPos,fPos,cPos,oPos,tPos},freqVals]= mtspectrumc(analogData(goodPos,blPos)',params);
                            
                        end
                    end 
                end   
            end 
        end 
    end
    if saveFlag
        save(fileSave,'powerST','powerBL','freqVals','electrodeList','parameterCombinations');
    end  
end 
end 

function [logPwrChangeOut,meanPwrSTOut,meanPwrBLall,stBandPwrOut,blBandPwr] = getPowerChanges(powerST,powerBL,freqVals,fBands,fBad)
% function to get change in power in given ranges. fBands is {[range1] [range2]}

[numElecs,numStimuli1,numStimuli2,numStimuli3] = size(powerBL);

meanPwrBL    = zeros([numElecs,numStimuli1,numStimuli2,numStimuli3,length(freqVals)]);
meanPwrST    = zeros([numElecs,numStimuli1,numStimuli2,numStimuli3,length(freqVals)]);
logPwrChange = zeros(numElecs,numStimuli1,numStimuli2,numStimuli3,length(fBands));
stBandPwr    = zeros(numElecs,numStimuli1,numStimuli2,numStimuli3,length(fBands));
blBandPwr    = zeros(numElecs,length(fBands));
[~,fBadInds] = intersect(freqVals,fBad);
    
clear meanPwrBLall   
for elec = 1:numElecs
    for pos1 = 1:numStimuli1  
        for pos2 = 1:numStimuli2       
            for pos3 = 1:numStimuli3    
                allLFPDataBL = (powerBL{elec,pos1,pos2,pos3});
                allLFPDataST = (powerST{elec,pos1,pos2,pos3});
                meanPwrBL(elec,pos1,pos2,pos3,:) = mean(allLFPDataBL,2);          % mean across trials
                meanPwrST(elec,pos1,pos2,pos3,:) = mean(allLFPDataST,2);  
                for ff = 1:length(fBands)
                    f_inds = freqVals>=fBands{ff}(1) & freqVals<=fBands{ff}(2);
                    f_inds(fBadInds) = false;
                    stBandPwr(elec,pos1,pos2,pos3,ff) = sum(squeeze(meanPwrST(elec,pos1,pos2,pos3,f_inds))); 
                end 
            end 
        end 
    end % end stimulus
    meanPwrBLall(elec,:) = squeeze(nanmean(nanmean(nanmean(meanPwrBL(elec,:,:,:,:),4),3),2));        % get average baseline across all stimuli
end
for ff = 1:length(fBands)
    f_inds = freqVals>=fBands{ff}(1) & freqVals<=fBands{ff}(2);
    f_inds(fBadInds)= false;
    blBandPwr(:,ff) = squeeze(sum(meanPwrBLall(:,f_inds),2));
    logPwrChange(:,:,:,:,ff) = log10(stBandPwr(:,:,:,:,ff)) - log10(repmat(blBandPwr(:,ff),[1 numStimuli1 numStimuli2 numStimuli3]));
end

logPwrChangeOut= squeeze(logPwrChange);
meanPwrSTOut   = squeeze(meanPwrST);
meanPwrBLall   = squeeze(meanPwrBLall);
stBandPwrOut   = squeeze(stBandPwr);
blBandPwr      = squeeze(blBandPwr);

end

