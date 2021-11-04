% MIBIgetNNThreshold
% Interactive script for getting thresholds for nn-filtering

% parameters
corePath = {'SampleData/extracted/Point1/','SampleData/extracted/Point1/'}; % path to cores that you want to evaluate for noise rduction. Specify several paths by separating with commas
massPath = 'SampleData/SamplePanel.csv'; % path to panel csv
load_data = 1; % after the first time that you run the script you can change to 0 to save the loading time.
plotChannel = 'CD8'; % channel that you want to denoise.
new_channel = 1; % after the first time that you run the script for a specific channel you can change to 0 to save the calculation time.
t = 3.5; % threshold used for separating signal and noise. Play with this number until you're happy with the results.
capImage = 10; % capping value for plotting. Set to lower to see dynamic range of low-abundant antigens
K = 25; % number of neighbors to use for density calculation. Usually can be kept as 25.

%% script
massDS = MibiReadMassData(massPath);
coreNum= length(corePath);
vec=[1:coreNum];
[~, plotChannelInd] = ismember(plotChannel,massDS.Label);

% load data. Do only for the first run of the script
if load_data
    p=cell(coreNum,1);
    for i=vec
        disp(['Loading core number ', num2str(i)]);
        p{i}=load([corePath{i},'dataNoBg.mat']);
    end

    disp('finished loading');
end

%get the NN values for the channel for all cores
if new_channel
    for i=vec
        p{i}.IntNormD{plotChannelInd}=MibiGetIntNormDist(p{i}.countsNoBg(:,:,plotChannelInd),p{i}.countsNoBg(:,:,plotChannelInd),K,2,K);
    end
end

chanelInAllPoints = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);
chanelInAllPointsCapped = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);

% plot the NN histograms for all points in a single plot, use to find
% noiseT cutoff
f=figure;
hedges = [0:0.25:30];
hline=zeros(coreNum,length(hedges)-1);
for j = vec
    data = p{j}.IntNormD{plotChannelInd};
    h=histogram(data,hedges,'Normalization','probability');
    hline(j,:)=h.Values;
end
clear('h','data');
a = 1:coreNum ;
labels = strread(num2str(a),'%s');
plot(hedges([1:end-1]),hline);
legend(labels);
plotbrowser on;

% test the threshold
for i=1:coreNum
    countsNoNoise{i} = MibiFilterImageByNNThreshold(p{i}.countsNoBg(:,:,plotChannelInd),p{i}.IntNormD{plotChannelInd},t);
    MibiPlotDataAndCap(p{i}.countsNoBg(:,:,plotChannelInd),capImage,['Core number ',num2str(i), ' - Before']); plotbrowser on;
    MibiPlotDataAndCap(countsNoNoise{i},capImage,['Core number ',num2str(i), ' - After. T=',num2str(t)]); plotbrowser on;
end
