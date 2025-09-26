% function to remove background signal

corePath = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master';
% params for filtering Arginase1
bgChannel = 'Au'; % Channel causing bg
targetChannels = {'CD4','CD14','CD206','CD16','CD20','CD209','Foxp3','GrzB','IFNg'}; % Channel with bg
cap=50;
t=0.21;
removeVal=200;
gausRad=1;
% coreNum = 41;
coreNum=[3,10,11,21,25];


for i=1:length(coreNum)
    point=coreNum(i);
    disp(['Working on ' num2str(point)]);
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
%     for j=1:length(targetChannels)
%         [~,targetChannelInd] = ismember(targetChannels{j},massDS.Label);
        targetChannelInd = 7;
        figure; imagesc(countsNoBg(:,:,targetChannelInd)); plotbrowser on; 
        figure; imagesc(countsNoBg2(:,:,targetChannelInd)); plotbrowser on;
%     end
    countsNoBg = countsNoBg2;
    save ([corePath,'/Point',num2str(point),'/data.mat'],'massDS','pointNumber','countsNoBg','countsAllSFiltCRSum');
    MibiSaveTifs ([corePath,'/Point',num2str(point),'/TIFsNoBg/'], countsNoBg, massDS.Label)
    close all;
end