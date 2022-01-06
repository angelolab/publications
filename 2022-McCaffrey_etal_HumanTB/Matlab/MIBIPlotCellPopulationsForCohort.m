% MIBIPlotCellPopulationsForCohort.m
% Author: Leeat Keren (original script:
% MibiPlotImmunePopulationsForCohort170911.m) adapted by Erin McCaffrey
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its phenotype.
% Keeps the colors consistent with those used for plotting in all other
% figures. 

%% 1. Import data and define points

path = [cd,'/data/']; % path to sc data and colorkey
pathSegment = [cd,'/data/segmentation/']; % path to segmentation masks
pathResults = [path,'/CPMs'];
if ~exist(pathResults)
    mkdir(pathResults)
end

points = [21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,...
          26,27,40,61,47,48,54,55,92,93]; %TB cohort data

% import data
sc_data = [path,'allTB-sarcoid-scdata_matlab.csv'];
opts = detectImportOptions(sc_data);
dataAll = readtable(sc_data,opts);

%converting to matrix for easy indexing
dataAllMat=table2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(:,[1,2,50])); % keep just SampleID, label, and cluster

%% 2. Define some useful variables (FlowSOM ID column, # phenos)

phenoCol=3;
numPheno=length(unique(dataAllMat(:,phenoCol)));

%% 3. Load pheno key and define the num/name of phenos
phenoKey=readtable([path,'cellpheno_numkey.csv']); 

%% 4. Define color code for pheno and major lineage (want both images)

%load color key with hexcodes and RGB values
colorKey=readtable([path,'colorkey_R.csv']);

% go through all cell types and produce a key with numerical codes for
% pheno and color
phenoCodes=[];
phenoColors=[];
for i=1:numPheno
    phenoCodes(i,1)=i;
    pheno = string(phenoKey(phenoKey.Code == i, :).Pheno);
    phenoR = colorKey(colorKey.imm_order == pheno, :).R;
    phenoG = colorKey(colorKey.imm_order == pheno, :).G;
    phenoB = colorKey(colorKey.imm_order == pheno, :).B;
    phenoRGB = [phenoR, phenoG, phenoB];
    phenoColors(i,1:3) = phenoRGB;
end

%combine into single matrix
phenoRGBkey = [phenoCodes,phenoColors];

% colorcode:
% [ 0 bg - white, cells - rainbow ]
cmap = phenoColors;
cmap01 = cmap/255;

% separate color code for lineage map
% [immune, endothelial, epithelial, fibroblast]
cmaplin = [51,102,255; 204,51,51; 255,153,102; 153,255,153];
cmaplin01 = cmaplin/255;

%% 5. Recolor cell masks by cell type

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
    imageLlin = zeros(size(newLmod,1)*size(newLmod,2),1);
    %go through all labels to determine color assignment
    for j=1:labelNum
        %get index of current cell
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % assign a color based on pheno
            cellVal = currCellData(cellInd,phenoCol);
            imageL(stats(j).PixelIdxList)=cellVal;
            % modify assignment for lineage map
            %endothelial
            if cellVal == 1 
                imageLlin(stats(j).PixelIdxList)=2;
            %epithelial 
            elseif cellVal == 18
                imageLlin(stats(j).PixelIdxList)=3;
            %fibroblast
            elseif cellVal == 8
                imageLlin(stats(j).PixelIdxList)=4;
            %immune
            else 
                imageLlin(stats(j).PixelIdxList)=1;
            end
        end
    end
    %shape back into an image and store in image stack
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    imageLlinReshape = reshape(imageLlin,size(newLmod,1),size(newLmod,2));
    % plot single
    f1=figure;
    colormap('parula');
    imagesc(label2rgb(imageLReshape,cmap01,'w'));
    title(['Point ', num2str(point)])
    imwrite(label2rgb(imageLReshape(31:994,31:994),cmap01,'w'),[pathResults,'/Point',num2str(point),'_CPM.tif'],'tif');
    f2=figure;
    imagesc(label2rgb(imageLlinReshape,cmaplin01,'w'));
	title(['Point ', num2str(point)])
    imwrite(label2rgb(imageLlinReshape(31:994,31:994),cmaplin01,'w'),[pathResults,'/Point',num2str(point),'_linCPM.tif'],'tif');
end
