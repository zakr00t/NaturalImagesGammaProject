% This program contains the database of all experiments. Optionally, this
% program also returns the index
function [experimentalDetails,matchIndex] = getExperimentalDetails(expDate)

if ~exist('expDate','var');         expDate=[];                         end

experimentalDetails{1} = [{'alpaH'}     {'Animals'}     {'Flora'}       {'ImagesAF'}        {'010817'}	{'GRF_001'}];
experimentalDetails{2} = [{'alpaH'}     {'Humans'}      {''}            {'ImagesH'}         {'080817'}	{'GRF_001'}];
experimentalDetails{3} = [{'alpaH'}     {'HumansBW'}    {''}            {'ImagesH_gray'}    {'110817'}	{'GRF_003'}];
experimentalDetails{4} = [{'alpaH'}     {'AnimalsBW'}   {'FloraBW'}     {'ImagesAF_gray'}	{'230817'}	{'GRF_001'}];
experimentalDetails{5} = [{'alpaH'}     {'Texture'}     {'Landscape'}	{'ImagesTL'}        {'240817'}	{'GRF_002'}];
experimentalDetails{6} = [{'alpaH'}     {'TextureBW'}	{'LandscapeBW'}	{'ImagesTL_gray'}	{'250817'}	{'GRF_002'}];
experimentalDetails{7} = [{'alpaH'}     {'Animals'}     {'AnimalsS'}	{'ImagesAS'}        {'300817'}	{'GRF_002'}];
experimentalDetails{8} = [{'alpaH'}     {'Texture'}     {'TextureS'}	{'ImagesTS'}        {'310817'}	{'GRF_003'}];
				
experimentalDetails{9} = [{'kesariH'}	{'Animals'}     {'AnimalsS'}	{'ImagesAS'}        {'020118'}	{'GRF_001'}];
experimentalDetails{10}= [{'kesariH'}	{'Texture'}     {'Landscape'}	{'ImagesTL'}        {'030318'}	{'GRF_003'}];
experimentalDetails{11}= [{'kesariH'}	{'Texture'}     {'TextureS'}	{'ImagesTS'}        {'040118'}	{'GRF_001'}];
experimentalDetails{12}= [{'kesariH'}	{'HumansBW'}    {''}            {'ImagesH_gray'}    {'050118'}	{'GRF_003'}];
experimentalDetails{13}= [{'kesariH'}	{'Humans'}		{''}            {'ImagesH'}         {'050318'}	{'GRF_003'}];
experimentalDetails{14}= [{'kesariH'}	{'Animals'}     {'Flora'}       {'ImagesAF'}        {'060318'}	{'GRF_002'}];
experimentalDetails{15}= [{'kesariH'}	{'AnimalsBW'}	{'FloraBW'}     {'ImagesAF_gray'}	{'150118'}	{'GRF_001'}];
experimentalDetails{16}= [{'kesariH'}	{'TextureBW'}	{'LandscapeBW'}	{'ImagesTL_gray'}	{'170218'}	{'GRF_001'}];

if ~isempty(expDate)
    numExperiments = length(experimentalDetails);
    allMatchIndex = zeros(1,numExperiments);
    for i=1:numExperiments
        allMatchIndex(i) = isequal(expDate,experimentalDetails{i}{5});
    end
    matchIndex = find(allMatchIndex==1);
    if isempty(matchIndex)
        disp('No matching expDate');
    elseif length(matchIndex)>1
        disp('Multiple expDates found');
    end
else
    matchIndex=[];
end