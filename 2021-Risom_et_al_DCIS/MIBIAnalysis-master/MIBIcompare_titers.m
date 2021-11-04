% MIBIcompare_titers
% Interactive script for choosing titers

corePath = {'SampleData/extracted/Point1/', ...
    'SampleData/extracted/Point2/'}; % cores to work on. Can add several paths, separated by commas.
Headers = {'High','Low'}; % Headers describing each one of the points. Will be used for visualization.
channel = {'CD8'}; % Channel to work on
cap = 5; % Capping value for plotting.
K=25; %Nearest neighbours to use for density estimation
First =1; % 1- If this is the first time running. If the script is sloa, you can change to 0 to save the loading time after the first run.

%% script
coreNum = length(corePath);
% load all cores
if First == 1
    p=cell(coreNum,1);
    for i=1:coreNum
        p{i} = load([corePath{i},'dataNoBg.mat']);
    end
end

[~,channelInd] = ismember(channel,p{1}.massDS.Label);

% 1. plot titration with same cap
for i=1:coreNum
    MibiPlotDataAndCap(p{i}.countsNoBg(:,:,channelInd),cap,[channel , ' - ' , Headers{i}]); plotbrowser on;
end

% 2. plot intensity histograms
figure;
for i=1:coreNum
    currData = p{i}.countsNoBg(:,:,channelInd);
    currDataLin = currData(:);
    currDataLin(currDataLin == 0) = [];
    hold on;
    histogram(currDataLin,'Normalization','probability','DisplayStyle','stairs');
    xlabel('Intensity');
    ylabel('Counts');
    plotbrowser on;
end

% 3. calculate NN histograms and plot
for i=1:coreNum
    p{i}.IntNormD{channelInd}=MibiGetIntNormDist(p{i}.countsNoBg(:,:,channelInd),p{i}.countsNoBg(:,:,channelInd),K,2,K);
end

figure;
for i=1:coreNum
    hold on;
    histogram(p{i}.IntNormD{channelInd},'DisplayStyle','stairs');
    xlabel('Mean distance to nearest neighbours');
    ylabel('Counts');
    plotbrowser on;
end