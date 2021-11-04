function MibiSpectraNoBareRegion(fileMSDF,fileNameXML,fileNameMass,countFile,specFile,pointNumber,figureBaselineFileName,figureSpectraFileName,figureBaselineFiltFileName,figureSpectraFiltFileName,figureSpectraRegularAndFiltFileName,figureRatioFileName,removeChannels,calibrateSpectra,spectraVec,figureTofFileName,fileMSDFAccum,calibrateSpectraPerDepth,CalibrationLogFileName)
% function MibiSpectraNoBareRegion(countFile,specFile,removeChannels)
% Function displays the spectra for an image, excluding the bare areas
% This function gets baseline by binning to a predefined region instead of
% finding minimal peaks

HEADER_SIZE = 3;

% get mass calibration parameters
paramS = MibiReadRunXml(fileNameXML,pointNumber);

% get mass edges
massDS = MibiReadMassData(fileNameMass);
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
% disp ('reading sat');
% tic
SAT= nan(numberOfPixels,2);
for i=1:numberOfPixels
    SAT(i,1)= fread(fileID,1,'int64');
    SAT(i,2)= fread(fileID,1,dataType);
end
% toc

%% get spectra
disp ('get spectra');
% tic
tmp=fread(fileID,Inf,'uint16');
fclose(fileID);
spectralData(:,1)= tmp([1:2:end]);
spectralData(:,2)= tmp([2:2:end]);
clear('tmp');
% convert spectra to mass counts
spectralDataDup = repelem(spectralData(:,1),spectralData(:,2));
if (calibrateSpectra==1)
    if calibrateSpectraPerDepth
        spectraVec = MibiGetSpectraCalibrationVecFromMSDF2(fileMSDFAccum,paramS,massDS,spectraVec);
    end
    sol= MibiCalibrateSpectrum (spectraVec);
    paramS.MassGain=eval(sol.a);
    paramS.MassOffset=eval(sol.b);
    if length(paramS.MassGain)>1
        paramS.MassGain=paramS.MassGain(1);
        paramS.MassOffset=paramS.MassOffset(1);
    end
    % write values to file
    fileID = fopen(CalibrationLogFileName,'a');
    fprintf(fileID,[fileMSDFAccum,'\t']);
    fprintf(fileID,'%f\t%f\t%f\t%f\t%f\t%f\n',[spectraVec,paramS.MassGain,paramS.MassOffset]);
    fclose(fileID);
end
mass = (( spectralDataDup * paramS.TimeResolution / 1000 - paramS.MassOffset ) / paramS.MassGain).^2;

% plot mass
massEdges= [0:0.05:250];
h=histogram(mass, massEdges);
f=figure;
plot(h.BinEdges([1:end-1]),h.Values);
xlabel('m/z');
saveas(f,figureSpectraFileName);

clear('spectralDataDup');
% toc

%% assign masses to pixels

% according to the mass file, create a 3d matrix XxYxZ
% with X,Y being the size of the image, and each Z is a different label
disp ('assign mass to pixel');
% tic
edges=[massDS.BaselineStart,massDS.BaselineEnd,massDS.Start,massDS.Stop]';
edges=edges(:);

counts=zeros(numberOfPixels,massNum*2);
spectralPos=[0;cumsum(SAT(:,2))];
for i=1:length(SAT)
    currSpectralData=spectralData([spectralPos(i)+1:spectralPos(i+1)],:);
    currSpectralDataDup = repelem(currSpectralData(:,1),currSpectralData(:,2));
    currMass = (( currSpectralDataDup * paramS.TimeResolution / 1000 - paramS.MassOffset ) / paramS.MassGain).^2;
    currCounts = histcounts(currMass,edges);
    counts(i,:) = currCounts([1:2:end]);
end
disp ('reshape');
baseline= counts(:,[1:2:end]);
counts= counts(:,[2:2:end]);
counts=reshape(counts,[paramS.XSize,paramS.YSize,massNum]);
baseline=reshape(baseline,[paramS.XSize,paramS.YSize,massNum]);
counts=permute(counts,[2,1,3]);
baseline=permute(baseline,[2,1,3]);

%% remove signal from pixels with background signal
% generate counts with removal of background channels
[tf,channelsLoc] = ismember(removeChannels,massDS.Label);
% add pixels and get connected components
pixelsToRemove = sum(counts(:,:,channelsLoc),3);
maskToRemove = pixelsToRemove >0;
maskToRemove3d = repmat(maskToRemove,1,1,length(massDS));
countsFilt = counts;
%countsFilt(maskToRemove3d) = 0.5*countsFilt(maskToRemove3d);
countsFilt(maskToRemove3d) = NaN;
baselineFilt = baseline;
%baselineFilt(maskToRemove3d) = 0.5*baselineFilt(maskToRemove3d);
baselineFilt(maskToRemove3d) = NaN;

%% calculate totalIon, Peaks and Baseline for entire image for filtered and non-filtered data
% total ion
totalIon = nansum(counts,3);
totalIonFilt = nansum(countsFilt,3);
% peaks
peaks = nansum(nansum(counts));
peaks = reshape(peaks,length(massDS),1,1);
peaksFilt = nansum(nansum(countsFilt));
peaksFilt = reshape(peaksFilt,length(massDS),1,1);
% baseline
baselineVal = nansum(nansum(baseline));
baselineVal = reshape(baselineVal,length(massDS),1,1);
baselineValFilt = nansum(nansum(baselineFilt));
baselineValFilt = reshape(baselineValFilt,length(massDS),1,1);


close all
%% save
save(countFile,'counts','countsFilt','peaks','peaksFilt','baseline','baselineFilt','baselineVal','baselineValFilt','totalIon','totalIonFilt','massDS');
% DEBUG: save spectra
% save(specFile,'mass','massFiltered','massDS');