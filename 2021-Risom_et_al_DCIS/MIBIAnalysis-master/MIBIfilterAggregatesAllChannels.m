% Filter aggregates for all cores in the study

% params
corePath = 'SampleData/extracted/cleanData'; % path to points to filter. Script assumes points in folder are named Point1, Point2 etc.
massPath = 'SampleData/SamplePanel.csv'; % path to panel csv. Make sure it has a column 'AggFilter' which has the aggregate threshold for each channel
coreNum = 2; % number of cores to work on
gausRad = 1; % radius for gaussian

%% script
massDS = MibiReadMassData(massPath);
for i=1:coreNum
    disp(['Working on point ',num2str(i)]);
    load([corePath,'/Point',num2str(i),'/dataDeNoiseCohort.mat']);
    countsNoNoiseNoAgg = countsNoNoise;
    for j=1:length(massDS)
        gausFlag = massDS.GausFlag(j);
        t = massDS.AggFilter(j);
        countsNoNoiseNoAgg(:,:,j) = MibiFilterAggregates(countsNoNoise(:,:,j),gausRad,t,gausFlag);
    end
    save([corePath,'/Point',num2str(i),'/dataNoAgg.mat'],'countsNoNoiseNoAgg');
    mkdir([corePath,'/Point',num2str(i)]);
    MibiSaveTifs ([corePath,'/Point',num2str(i),'/TIFsNoAgg/'], countsNoNoiseNoAgg, massDS.Label);
end

