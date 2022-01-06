% MIBIPlotCellMEsForCohort.m
% Author: Leeat Keren (original script:
% MibiPlotImmunePopulationsForCohort170911.m) adapted by Erin McCaffrey
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

points = [21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,...
          26,27,40,61,47,48,54,55,92,93]; %TB cohort data

% import data
sc_data = [path,'allTB-ME_annotated.csv'];
opts = detectImportOptions(sc_data);
dataAll = readtable(sc_data,opts);

%converting to matrix for easy indexing
dataAllMat=table2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(:,[1,2,57])); % keep just SampleID, label, and max ME

%% 2. Load ME key, define # of MEs, and ME col
MEKey=readtable([path,'ME_numkey.csv']); 
numMEs = length(MEKey.ME);
MEcol = 3;
%% 3. Define color code for topic

%load color key with hexcodes and RGB values
colorKey=readtable([path,'colorkey_MEs.csv']); 

% go through all cell types and produce a key with numerical codes for
% pheno and color
MECodes=[];
MEColors=[];
for i=1:numMEs
    MECodes(i,1)=i;
    ME = string(MEKey(MEKey.Code == i, :).ME);
    phenoR = colorKey(colorKey.ME == ME, :).R;
    phenoG = colorKey(colorKey.ME == ME, :).G;
    phenoB = colorKey(colorKey.ME == ME, :).B;
    phenoRGB = [phenoR, phenoG, phenoB];
    MEColors(i,1:3) = phenoRGB;
end

%combine into single matrix
phenoRGBkey = [MECodes,MEColors];

% colorcode:
% [ 0 bg - white, cells - rainbow ]
cmap = MEColors;
cmap01 = cmap/255;

%% 4. Recolor cell masks by ME assignment

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
    %go through all labels to determine color assignment
    for j=1:labelNum
        %get index of current cell
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % assign a color based on pheno
            cellVal = currCellData(cellInd,MEcol)+1;
            imageL(stats(j).PixelIdxList)=cellVal;
        end
    end
    %shape back into an image and store in image stack
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    % plot single
    f1=figure;
    imagesc(label2rgb(imageLReshape,cmap01,'w'));
    title(['Point ', num2str(point)])
    imwrite(label2rgb(imageLReshape(31:994,31:994),cmap01,'w'),[pathResults,'/Point',num2str(point),'_MaxPM.tif'],'tif');
end