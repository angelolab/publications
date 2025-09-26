% remove background for all cores in the study

corePath = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master';


% params for filtering Arginase1
bgChannel = 'HH3'; % Channel causing bg
targetChannels = {'MastChyTry'}; % Channel with bg
cap=5;
t=0.2;
removeVal=20;
gausRad=0.1;
coreNum = 41;


for i=1:coreNum
    point=i;
    disp(['Working on ' num2str(i)]);
    % load in data
    load([corePath,'/Point',num2str(point),'/data.mat']); 
    % get background channel data
    [~,bgChannelInd] = ismember(bgChannel,massDS.Label);
    % produce mask
    mask = MibiGetMask(countsNoBg(:,:,bgChannelInd),cap,t,gausRad);
    % make copy of bg-subtracted data
    countsNoBg2=countsNoBg;
    % remove signal from target channel(s)
    for j=1:length(targetChannels)
        [~,targetChannelInd] = ismember(targetChannels{j},massDS.Label);
        countsNoBg2(:,:,targetChannelInd) = MibiRemoveBackgroundByMaskSingleChannel(countsNoBg(:,:,targetChannelInd),mask,removeVal);
    end
    figure; imagesc(mask); plotbrowser on; 
    for j=1:length(targetChannels)
        [~,targetChannelInd] = ismember(targetChannels{j},massDS.Label);
        figure; imagesc(countsNoBg(:,:,targetChannelInd)); plotbrowser on; 
        figure; imagesc(countsNoBg2(:,:,targetChannelInd)); plotbrowser on;
    end
    countsNoBg = countsNoBg2;
    save ([corePath,'/Point',num2str(point),'/data.mat'],'massDS','pointNumber','countsNoBg','countsAllSFiltCRSum');
    MibiSaveTifs ([corePath,'/Point',num2str(point),'/TIFsNoBg/'], countsNoBg, massDS.Label)
    close all;
end