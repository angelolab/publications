% PAHPlotNeighborhoodAssignmentForCohort.m
% Author: Leeat Keren, modified by Erin McCaffrey
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its phenotype.
% Keeps the colors consistent with those used for plotting in all other
% figures. 

%% Initiate paths to data and define important variables

path = '';
pathSegment = '';
resultsDir = '';
mkdir(resultsDir);

% define points for analysis
points = linspace(42,90,49); 
imSize=1024;
dataAll=dataset('File',[path,'PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples

%convert to matrix for easy indexing
dataAllCell=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllCell(2:size(dataAllCell,1),[1,2,16]));

% annotate columns 
pointCol = 1;
cellLabelCol = 2;

% column containing the numerically encoded phenotype
phenoCol=3;
numPheno=length(unique(dataAllMat(:,phenoCol)));

%% Define color code for phenotypes

%load color key with hexcodes and RGB values
colorKey=dataset('File',[path,'neighborhood_colorkey.csv'],'Delimiter',',');

% go through all cell types and produce a key with numerical codes for
% pheno and color
phenoCodes=[];
phenoColors=[];
for i=1:numPheno
    phenoCodes(i,1)=i;
    pheno = colorKey(i, :).neighborhood;
    phenoR = colorKey(pheno, :).R;
    phenoG = colorKey(pheno, :).G;
    phenoB = colorKey(pheno, :).B;
    phenoRGB = [phenoR, phenoG, phenoB];
    phenoColors(i,1:3) = phenoRGB;
end

%combine into single matrix
phenoRGBkey = [phenoCodes,phenoColors];

% colorcode:
% [ 0 bg - white, cells - rainbow ]
cmap = phenoColors;
cmap01 = cmap/255;

%% Recolor cell masks by cell type

for i=1:length(points)

    %load data
    disp(points(i));
    point = points(i);
    seg_mask = double(imread([pathSegment,'/Point',num2str(point),'_labels.tiff']));

    %get stats for objects in image
    stats = regionprops(seg_mask,'PixelIdxList');

    %get just data for current point
    currInds = (dataAllMat(:,pointCol) == point);
    currCellData = dataAllMat(currInds,:);

    %get labels of objects in image
    labels = currCellData(:,cellLabelCol);
    labelNum = length(labels);

    %create vector of all pixels for assigning color map
    px_label = zeros(size(seg_mask,1)*size(seg_mask,2),1);

    %go through all labels to determine color assignment
    for j=1:labelNum
        label = labels(j);
        
        %get index of current cell
        cellInd = find(ismember(currCellData(:,cellLabelCol),label));

        % assign a color based on pheno
        if ~isempty(cellInd)
            cellVal = currCellData(cellInd,phenoCol);
            px_label(stats(label).PixelIdxList)=cellVal;
        end
    end

    %shape back into an image and store in image stack
    px_label_reshape = reshape(px_label,size(seg_mask,1),size(seg_mask,2));

    % plot 
    f1=figure;
    colormap('parula');
    imagesc(label2rgb(px_label_reshape,cmap01,'w'));
    title(['Point ', num2str(point)])

    % export 
    imwrite(label2rgb(px_label_reshape,cmap01,'w'),[resultsDir,'/Point',num2str(point),'_NC4.tif'],'tif');
end
