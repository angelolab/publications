% % MIBICreateDistanceMatrixBetweenCells
% Author: Leeat Keren, modified by Erin McCaffrey
% For each point, reads in the massDS and the segmentation path containing
% the relevant segmentation params. It next produces a distance matrix of
% the center of all cells to all other cells for downstream spatial
% analysis. 

%% Define path and points to analyze

% path to the segmentation masks
path_segment = '';

% points to create distance matrix for
points = linspace(1,90,90); 

%% Create distance matrix

for i=1:length(points)
    point=points(i);
    disp(['point',num2str(point)]);

    % load the segmentation
    seg_mask = imread([path_segment,'/Point',num2str(point),'_labels.tiff']);

    % get centers of cells
    stats = regionprops(seg_mask,'centroid');

    % get x,y coordinates in a mat
    centroidCoord = zeros(length(stats),2);
    for c=1:length(stats)
        centroidCoord(c,1) = stats(c).Centroid(1);
        centroidCoord(c,2) = stats(c).Centroid(2);
    end

    % create a matrix of distances between all cells
    distancesMat = pdist2(centroidCoord,centroidCoord);

    % save as a Matlab object
    save([path_segment,'/Point',num2str(point),'_cellDistances.mat'],'distancesMat');
end
