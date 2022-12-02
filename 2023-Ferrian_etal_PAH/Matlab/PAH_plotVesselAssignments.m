% PAHPlotVesselAssignments.m
% Author: Erin McCaffrey
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its
% association to a particle region of the blood vessel mask.

%% Initiate paths to data and define important variables

path = '';
pathSegment = '';
resultsDir = [pathSegment,'/vessel_region_overlays/'];
mkdir(resultsDir);
points =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,...
    26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,...
    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,...
    70,71,72,73,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,90];

dataAll=dataset('File',[path,'celldata_region_annotated.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%convert to matrix for easy indexing
dataAllCell=dataset2cell(dataAll);
keep_cols = [1,41,55,56];
dataAllMat=cell2mat(dataAllCell(2:size(dataAllCell,1),keep_cols));

% annotate columns 
pointCol = 1;
cellLabelCol = 2;

% column containing the region data
regionCol=3;
regionPeriphCol=4;

%% Define color codes

% colorcode (three assignments, RGB values):
cmap = ([200, 200, 199 ; 198, 112, 81; 101, 30, 56]);
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
    px_labels = zeros(size(seg_mask,1)*size(seg_mask,2),1);

    %go through all labels to determine color assignment
    for j=1:labelNum
        label = labels(j);

        %get index of current cell
        cellInd = find(ismember(currCellData(:,cellLabelCol),label));
        if ~isempty(cellInd)

            % assign a color based on regoion distribution
            regionVal = currCellData(cellInd,regionCol);
            regionPeriphVal = currCellData(cellInd,regionPeriphCol);

            % option 1 = neither in vessel or on border
            if regionVal == 0 && regionPeriphVal == 0
                px_labels(stats(label).PixelIdxList)=1;
            % option 2 = on the vessel border
            elseif regionVal == 0 && regionPeriphVal == 1
                px_labels(stats(label).PixelIdxList)=2;
            % option 3 = in the vessel
            elseif regionVal == 1 
                px_labels(stats(label).PixelIdxList)=3;
            end
        end
    end
    
    %shape back into an image and store in image stack
    px_labels_reshape = reshape(px_labels,size(seg_mask,1),size(seg_mask,2));

    % plot
    f1=figure;
    colormap('parula');
    imagesc(label2rgb(px_labels_reshape,cmap01,'w'));
    title(['Point ', num2str(point)])

    % export as tif
    imwrite(label2rgb(px_labels_reshape,cmap01,'w'),[resultsDir,'/Point',num2str(point),'_region_assignment.tif'],'tif');
end
