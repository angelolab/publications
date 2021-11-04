%% A low-level Mibi analysis pipeline
% Stages include:
% A. Plot TOF for all frames to get user input on spectra calibration
% B. Extract MSDF for each depth
% C. Coregister depths and sum

%% default parameters

% Instrument is ran at depth profile mode
if ~exist('depthProfile')
    depthProfile=1;
    chemicalImage=0;
end
% Process msdf
if ~exist('processMsdf')
    processMsdf=1;
end

% don't remove channels
if ~exist('removeChannels')
    removeChannels={};
end

% % parameters for extracting spectra
if ~exist('extractSpectraParams')
    extractSpectraParams=1;
end

if ~exist('makePlots')
    makePlots=1;
end

if ~exist('fixRasterCols')
    fixRasterCols=0;
end

if ~exist('fixRasterRows')
    fixRasterRows=0;
end

if ~exist('calibrateSpectraPerDepth')
    calibrateSpectraPerDepth=1;
end

% parameters for summing depths
if ~exist('sumDepths')
    sumDepths=1;
end

% parameters for coregistration (if necessary)
if ~exist('coregister_planes')
    coregister_planes=1;
end

if ~exist('CRchannel')
    CRchannel='totalIon';
end

% parameters for noise filtering
if ~exist('filterNoise')
    filterNoise=0;
end

% parameters for final saves
if ~exist('plotCorrelations')
    plotCorrelations=0;
end

if ~exist('saveTifsNoNoise')
    saveTifsNoNoise=0;
end

if ~exist('saveTifs')
    saveTifs=1;
end

if ~exist('badDepths')
    badDepths = [];
end

tic
mkdir([processedDataDir,'/Point',num2str(pointNumber),'/Results/']);
CalibrationLogFileName = [processedDataDir,'/Point',num2str(pointNumber),'/Results/SpectraCalibrationLog.txt'];
figureCompareTOFFileName = [processedDataDir,'/Point',num2str(pointNumber),'/Results/TOFPerPlain.fig'];
figureCompareMassFileName = [processedDataDir,'/Point',num2str(pointNumber),'/Results/MassPerPlain.fig'];

if First
    MibiPlotSpectraPerDepthForTimeAndMass3(fileNameXML,fileNameMass,pointNumber,depthStart,depthEnd,figureCompareTOFFileName,figureCompareMassFileName,calibrateSpectra,calibrateSpectraPerDepth,spectraVec);
    %MibiPlotSpectraPerDepth(fileNameXML,fileNameMass,pointNumber,depthStart,depthEnd,figureCompareTOFFileName);
    disp ('Look at spectra to make sure that calibration values are ok and rerun with First=0');
else

%% process msdf for each plain
if processMsdf
    disp('Parsing msdf files');
    for d=depthStart:depthEnd
        disp(['processing depth',num2str(d-1)]);
        if depthProfile
            pathSt = [dataDir,'Point',num2str(pointNumber),'/RowNumber0/Depth_Profile0/Depth',num2str(d-1)];
        elseif chemicalImage
            pathSt = [dataDir,'Point',num2str(pointNumber),'/RowNumber0/Chemical_Image0'];
        else
            error ('If processing msdf, either depth profile or chemical image should be specified');
        end
        fileMSDF=[pathSt,'/Image.msdf'];
        fileMSDFAccum=[pathSt,'/ImageAccum.msdf'];
        countFile=[pathSt,'/counts.mat'];
        specFile=[pathSt,'/spectra.mat'];
        figureTofFileName=[pathSt,'/tof.fig'];
        figureBaselineFileName=[pathSt,'/baseline.fig'];
        figureSpectraFileName=[pathSt,'/spectra.fig'];
        figureBaselineFiltFileName=[pathSt,'/baselineFilt.fig'];
        figureSpectraFiltFileName=[pathSt,'/spectraFilt.fig'];
        figureSpectraRegularAndFiltFileName=[pathSt,'/spectraRegularAndFilt.fig'];
        figureRatioFileName=[pathSt,'/ratioSpectraRegularAndFilt.fig'];
        MibiSpectraNoBareRegion4(fileMSDF,fileNameXML,fileNameMass,countFile,specFile,pointNumber,figureBaselineFileName,figureSpectraFileName,figureBaselineFiltFileName,figureSpectraFiltFileName,figureSpectraRegularAndFiltFileName,figureRatioFileName,removeChannels,calibrateSpectra,spectraVec,figureTofFileName,fileMSDFAccum,calibrateSpectraPerDepth,CalibrationLogFileName);
    end
end

%% load data and find peaks and baseline. Plot if needed
disp('Loading data and getting baseline and peaks info');
if extractSpectraParams
    [massDS,countsAllS,countsAllSFilt,baselineVal,baselineValFilt,peaks,peaksFiltered,ratioAll,ratioFiltered,ratioFilteredToNonFiltered,depthsRemovedPerPixel,totalIon,totalIonFilt]= MibiAnalizeSpectra4(processedDataDir,pointNumber,depthStart,depthEnd,makePlots,depthProfile,chemicalImage);
else
    load([processedDataDir,'/Point',num2str(pointNumber),'/ParamsSpectra.mat']);
end
%% sum planes: coregister, shift and sum
if sumDepths
    goodDepths = [depthStart:depthEnd];
    goodDepths(badDepths)=[];
    countsAllSFilt = countsAllSFilt(goodDepths);
    baselineValFilt = baselineValFilt(goodDepths);
    totalIon = totalIon(:,:,goodDepths);
    totalIonFilt = totalIonFilt(:,:,goodDepths);
    
    % fix raster scanner
    if fixRasterCols==1 | fixRasterRows==1
        [countsAllSFilt,totalIon,totalIonFilt]= MibiFixRaster(countsAllSFilt,totalIon,totalIonFilt, fixRasterCols, fixRasterRows);
    end
  
    
    depthEnd=length(countsAllSFilt);
    if depthProfile
        % 1. Coregister depths
        % Calculate coregistration offset based on the first plane of each
        % depth profile
        if coregister_planes==1
            disp('Coregestering planes');
            if strcmp(CRchannel,'totalIon')
                [yd,xd]=planeCoregTotalIonByPrevDepth(totalIonFilt,[processedDataDir,'/Point',num2str(pointNumber),'/Results/CoregistrationOffsets.fig']);
            else
                [yd,xd]=MibiPlaneCoreg(CRchannel,massDS,countsAllSFilt,[processedDataDir,'/Point',num2str(pointNumber),'/Results/CoregistrationOffsets.fig']);
            end
            % shift coregistration params so that we will get the majority
            % of planes
            medXd = round(median(xd));
            medYd = round(median(yd));
            xsize = size(countsAllSFilt{1},1);
            ysize = size(countsAllSFilt{1},2);
            xd = xd + (xsize-medXd);
            yd = yd + (ysize-medYd);
            save([processedDataDir,'/Point',num2str(pointNumber),'/params.mat'],'xd','yd');
        else
            load([processedDataDir,'/Point',num2str(pointNumber),'/params.mat']);
        end

        % 2. Shift planes, sum and gaussian smooth
        %Sum planes for each channel
        % shift
        disp('shift planes');
        countsAllSFiltCR= cell(size(countsAllSFilt));
        for d=1:depthEnd
            disp(['Adding planes for depth',num2str(d)])
            countsAllSFiltCR{d} = zeros(size(countsAllSFilt{d}));
            for i=1:length(massDS)
                 countsAllSFiltCR{d}(:,:,i)=countsAllSFiltCR{d}(:,:,i)+double(applyOffset(xd,yd,d, ...
                     [size(countsAllSFilt{d},1), size(countsAllSFilt{d},2)],countsAllSFilt{d}(:,:,i)));
            end
        end

        % sum
        disp('sum planes');
        countsAllSFiltCRSum = zeros(size(countsAllSFilt{1}));
        for d=1:depthEnd
            countsAllSFiltCRSum = countsAllSFiltCRSum + countsAllSFiltCR{d};
        end
    elseif chemicalImage
        countsAllSFiltCR = countsAllSFilt;
        countsAllSFiltCRSum = countsAllSFilt{1};
        countsAllSFiltCRSum(isnan(countsAllSFiltCRSum)) = 0;
    end

    save([processedDataDir,'/Point',num2str(pointNumber),'/summedPlanes.mat'],'countsAllSFiltCR','countsAllSFiltCRSum','massDS','-v7.3');
else
    load([processedDataDir,'/Point',num2str(pointNumber),'/summedPlanes.mat']);
end
clear('countsAllSFilt');

% Plots
% all markers side by side
% MibiPlotAllMarkersSideNySide (massDS.Label,countsAllSFiltCRSum,1,1,99,1,1,{'Au'});
% saveas(gcf,['Point',num2str(pointNumber),'/Results/AllMarkersTiled.fig']);

% all markers intensity hist
f=plotMarkerIntensityHist (countsAllSFiltCRSum, massDS.Label);
saveas(f,[processedDataDir,'/Point',num2str(pointNumber),'/Results/AllMarkersIntensityHist.fig']);

% total Ion as a function of pixel
totalIonFiltSum=sum(totalIonFilt,3);
% top to bottom
dataTot=totalIonFiltSum';
dataTot=dataTot(:);
dataTotSm= smooth(dataTot,size(countsAllSFiltCRSum,1),'moving');
f=figure;
hold on;
bar(dataTot);
plot([1:length(dataTot)],dataTotSm,'r-');
title('Ion Gradient Top Bottom');
saveas(f,[processedDataDir,'/Point',num2str(pointNumber),'/Results/IonGradientTopBottom.fig']);
% left to right
dataTot=totalIonFiltSum;
dataTot=dataTot(:);
dataTotSm= smooth(dataTot,size(countsAllSFiltCRSum,1),'moving');
f=figure;
hold on;
bar(dataTot);
plot([1:length(dataTot)],dataTotSm,'r-');
title('Ion Gradient Left Right');
saveas(f,[processedDataDir,'/Point',num2str(pointNumber),'/Results/IonGradientLeftToRight.fig']);


%% save
save ([processedDataDir,'/Point',num2str(pointNumber),'/data.mat'],'totalIonFiltSum','countsAllSFiltCRSum','massDS','pointNumber');
toc
% save TIFs
if saveTifs
    mkdir([processedDataDir,'/Point',num2str(pointNumber),'/TIFs/']);
    for i=1:length(massDS)
        figure; imagesc(countsAllSFiltCRSum(:,:,i)); title(massDS.Label(i));
        data=uint8(countsAllSFiltCRSum(:,:,i));
        imwrite(data,[processedDataDir,'/Point',num2str(pointNumber),'/TIFs/',massDS.Label{i},'.tif']);
    end
    figure;
    imagesc(totalIonFiltSum);
    title('Total Ion');
    data=uint8(totalIonFiltSum);
    imwrite(data,[processedDataDir,'/Point',num2str(pointNumber),'/TIFs/totalIon.tif']);
    %close all;
end

toc
end