% where all of your MIBI data is stored, doesn't change from run to run
root_path = '/Users/noahgreenwald/Documents/MIBI_Data/';

% change this each time
run_name = 'Segmentation/Run1';

folder = [root_path, run_name, '/'];
cd (folder);
addpath(genpath(folder));
fileNameXML=[folder, run_name, '.xml' ];

% change this each time
fileNameMass= [folder, 'SegmentationPanel_180720.csv'];

dataDir=folder;
processedDataDir=folder;


%% default parameters
chemicalImage = 0; 
processMsdf = 1;
extractSpectraParams=1;
makePlots=1;
fixRaster=0;
calibrateSpectraPerDepth=1;
coregister_planes=1;
CRchannel='totalIon';
filterNoise=0;
plotCorrelations=0;
saveTifsNoNoise=0;
saveTifs=1;
removeChannels = {};

depthStart = 1;
depthEnd = 1;
depthProfile = 1;
sumDepths=1;
badDepths = [];

parfor i = 1:29
    pointNumber=i;
    % parameters for spectra calibration (if necessary)
    calibrateSpectra=1;
    First = 0;
    spectraVec=[3358,22.930,9500,196.960];

    MibiAnalysis3par(pointNumber, fileNameXML, fileNameMass, processedDataDir, dataDir, spectraVec, depthStart, depthEnd, depthProfile, chemicalImage, processMsdf, calibrateSpectra, First, removeChannels, extractSpectraParams, makePlots, fixRaster, calibrateSpectraPerDepth, sumDepths, coregister_planes, CRchannel, filterNoise, plotCorrelations, saveTifsNoNoise, saveTifs, badDepths)
end
% shut down parallel pool
delete(gcp('nocreate'))