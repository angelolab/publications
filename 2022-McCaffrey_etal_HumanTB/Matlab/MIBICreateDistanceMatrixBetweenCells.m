% % MIBICreateDistanceMatrixBetweenCells.m
% Author: Leeat Keren, modified by Erin McCaffrey
% For each point, reads in the massDS and the segmentation path containing
% the relevant segmentation params. It next produces a distance matrix of
% the center of all cells to all other cells for downstream spatial
% analysis. 

%% path to the segmentation data
path = [cd,'/data/segmentation/'];

%% points to create distance matrix for
points = [21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,...
          26,27,40,61,47,48,54,55,92,93];

%% Create distance matrices
for i=1:length(points)
    point=points(i);
    disp(['point',num2str(point)]);
    % load the segmentation mask
    newLmod = imread([path,'segmentationmask_SampleID',num2str(point),'.tif']);
    % get centers of cells
    stats = regionprops(newLmod,'centroid');
    % get x,y coordinates in a mat
    centroidCoord = zeros(length(stats),2);
    for c=1:length(stats)
        centroidCoord(c,1) = stats(c).Centroid(1);
        centroidCoord(c,2) = stats(c).Centroid(2);
    end
    % create a matrix of distances between all cells
    distancesMat = pdist2(centroidCoord,centroidCoord);
    save([path,'/Point',num2str(point),'_cellDistances.mat'],'distancesMat');
end
