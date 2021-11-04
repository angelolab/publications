function tof = MibiParseMsdfSpectrum(fileMSDF,paramS,massDS)
% function MibiSpectraNoBareRegion(cfileMSDF,paramS,massDS)
% Function displays the spectra for an image, excluding the bare areas

HEADER_SIZE = 3;

massNum = size(massDS,1);

%% parse header
fileID = fopen(fileMSDF);
header = fread(fileID,HEADER_SIZE,'uint32');
fileVariant = header(1);
binsPerSpectrum = header(2);
numberOfPixels = header(3);
if (fileVariant == 0)
    dataType = 'uint32';
elseif (fileVariant == 1)
    dataType = 'uint16';
else
    disp (['Unknown file variant: ',num2str(fileVariant)]);
    return;
end

%% read SAT
SAT= nan(numberOfPixels,2);
for i=1:numberOfPixels
    SAT(i,1)= fread(fileID,1,'int64');
    SAT(i,2)= fread(fileID,1,dataType);
end

%% get spectra
tmp=fread(fileID,Inf,'uint32');
fclose(fileID);
spectralData(:,1)= tmp([1:2:end]);
spectralData(:,2)= tmp([2:2:end]);
clear('tmp');
% duplicate values with more than one count in spectra
spectralDataDup = repelem(spectralData(:,1),spectralData(:,2));
tof = spectralDataDup * paramS.TimeResolution;