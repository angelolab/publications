% MIBIPlotCellMEGradientsForCohort.m
% Author: Erin McCaffrey (adapted from original script by Leeat Keren:
% MibiPlotImmunePopulationsForCohort170911.m) 
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its ME.
% Keeps the colors consistent with those used for plotting in all other
% figures. 

%% 1. Import data and define points

path = [cd,'/data/']; % path to sc data and colorkey
pathSegment = [cd,'/data/segmentation/']; % path to segmentation masks
pathResults = [path,'/MaxPMs'];
if ~exist(pathResults)
    mkdir(pathResults)
end

points = [33]; %selected point

% import data
sc_data = [path,'allTB-ME_annotated.csv'];
opts = detectImportOptions(sc_data);
dataAll = readtable(sc_data,opts);

%converting to matrix for easy indexing
dataAllMat=table2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(:,[1,2,49:57])); % keep just SampleID, label, and ME info

%% 2. Define ME to plot gradient for 
% (ie. will only give values to cells within a given topic, but will still 
% plot the probability instead of a binary assignment

% topic of interest
topic = 1;

% column containing the topping loadings for the given topic
topic_prob_col = topic + 3;

% column containing the topic assignment
topic_assn_col = 11;

%% 3. Produce ME loading map

for i=1:length(points)
    %load data
    disp(points(i));
    point = points(i);
    newLmod = imread([pathSegment,'segmentationmask_SampleID',num2str(point),'.tif']);
    %get stats for objects in image
    stats = regionprops(newLmod,'PixelIdxList');
    %get just data for current point
    currInds = (dataAllMat(:,1) == point);
    currCellData = dataAllMat(currInds,:);
    %get labels of objects in image
    currLabelIdentity = unique(newLmod);
    labelNum = length(currLabelIdentity);
    %create vector of all pixels for assigning color map
    imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
    for j=1:labelNum 
        %get index of current cell
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % assign a color based on topic assignment
            if(currCellData(cellInd,topic_assn_col) == topic)
                loading = currCellData(cellInd,topic_prob_col);
                imageL(stats(j).PixelIdxList)=loading;
            end
        end
    end
    %shape back into an image
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    % plot
    f1=figure;
    colormap('parula');
    imagesc(imageLReshape);
    imwrite(imageLReshape(31:994,31:994),[pathResults,'/Point',num2str(point),'_ME',num2str(topic),'_gradient.tif'],'tif');
end