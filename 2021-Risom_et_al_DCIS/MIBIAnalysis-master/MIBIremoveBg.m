% remove background for several cores according to removal params

corePath = {'SampleData/extracted/Point1/'}; % cores to work on. Can add several paths, separated by commas.

% put in here parameters from MIBIgetBgSubtractionParams.m
bgChannel = ('Background');
gausRad= 1;
t= 0.2;
removeVal= 2;
cap = 10;
coreNum = length(corePath);

for i=1:length(corePath)
    disp(['Working on ' num2str(i)]);
    load([corePath{i},'data.mat']);
    [~,bgChannelInd] = ismember(bgChannel,massDS.Label);
    mask = MibiGetMask(countsAllSFiltCRSum(:,:,bgChannelInd),cap,t,gausRad);
    countsNoBg = MibiRemoveBackgroundByMaskAllChannels(countsAllSFiltCRSum,mask,removeVal);
    save ([corePath{i},'dataNoBg.mat'],'massDS','pointNumber','countsNoBg');
    MibiSaveTifs ([corePath{i},'/TIFsNoBg/'], countsNoBg, massDS.Label)
    close all;
end