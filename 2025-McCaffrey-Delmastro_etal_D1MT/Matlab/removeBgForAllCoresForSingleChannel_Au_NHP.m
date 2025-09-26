
corePath = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master';


% params for filtering Arginase1
bgChannel = 'Au'; % Channel causing bg
targetChannels = {'pSMAD3','pS6','SMA','PD1','MastChyTry','CD206','CD209','GrzB','iNOS','CD40'}; % Channel with bg
cap=50;
t=0.1;
removeVal=50;
gausRad=1;
coreNum = 41;
% coreNum=[11,12,21,22,23,29,31,32,33,34,37,38];


for i=4:coreNum
    point=i;
    disp(['Working on ' num2str(point)]);
    % load in data
    load([corePath,'/Point',num2str(point),'/data.mat']); 
    % get background channel data
    [~,bgChannelInd] = ismember(bgChannel,massDS.Label);
    % produce mask
    mask = MibiGetMask(countsAllSFiltCRSum(:,:,bgChannelInd),cap,t,gausRad);
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
        targetChannelInd = 16;
        figure; imagesc(countsNoBg(:,:,targetChannelInd)); plotbrowser on; 
        figure; imagesc(countsNoBg2(:,:,targetChannelInd)); plotbrowser on;
%     end
    countsNoBg = countsNoBg2;
    save ([corePath,'/Point',num2str(point),'/data.mat'],'massDS','pointNumber','countsNoBg','countsAllSFiltCRSum');
    MibiSaveTifs ([corePath,'/Point',num2str(point),'/TIFsNoBg/'], countsNoBg, massDS.Label)
    close all;
end