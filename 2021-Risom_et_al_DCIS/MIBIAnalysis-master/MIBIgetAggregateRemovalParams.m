% MIBIgetAggregateRemovalParams
% Find thresholds for aggregate removal
% The script works by gaussian-smoothing the data and then removing connected components below a certain size. 

% parameters
corePath = {'SampleData/extracted/cleanData/Point1/','SampleData/extracted/cleanData/Point2/'}; % path to cores that you want to evaluate for aggregate removal. Specify several paths by separating with commas
massPath = 'SampleData/SamplePanel.csv'; % path to panel csv
load_data = 1; % after the first time that you run the script you can change to 0 to save the loading time.
plotChannel = 'CD4'; % channel that you want to work on.
gausFlag = 1; % flag of whether to do gaussian smoothing or not.
gausRad = 1; % gauss radius for smoothing.
capImage = 5; % value for capping the images when plotting. If colors are saturated, increase this number.

t = 100; % Important: threshold used for aggregate removal. Components smaller than this size will be removed. Play with this number until you're happy with the results.

%% script

massDS = MibiReadMassData(massPath);
coreNum= length(corePath);

% load data. Do only for the first run of the script
if load_data
    p=cell(coreNum,1);
    q=cell(coreNum,1);
    for i=1:coreNum
        disp(['Loading core number ', num2str(i)]);
        p{i}=load([corePath{i},'dataDeNoiseCohort.mat']);
    end

    disp('finished loading');
end

% perform aggregate removal:
[~, plotChannelInd] = ismember(plotChannel,massDS.Label);
for i=1:length(corePath)
    q{i}.countsNoNoiseNoAgg(:,:,plotChannelInd) = MibiFilterAggregates(p{i}.countsNoNoise(:,:,plotChannelInd),gausRad,t,gausFlag);

    % plot
    MibiPlotDataAndCap(p{i}.countsNoNoise(:,:,plotChannelInd),capImage,['Point ',num2str(i), ' - Before - ',massDS.Label{plotChannelInd}]); plotbrowser on;
    MibiPlotDataAndCap(q{i}.countsNoNoiseNoAgg(:,:,plotChannelInd),capImage,['Point ',num2str(i), ' - After - ',massDS.Label{plotChannelInd}]); plotbrowser on;

end

