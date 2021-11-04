%% MIBIgetCellNeighbors.m
% This script reads in the segmentation mask for a given sample. Next it
% for each cell it determines the neighbors of that cell based on shared
% borders as defined by expanding the cell object (and its border) by 1
% px and determining which other cells it shares a border with. In
% the end it saves a binary matrix of n cells x n cells where a value of 1
% (true) indicates that cell is a neighbor while a value of 0 (false means
% that cell is not a neighbor. It excludes itself (ie self v self = 0).
%% Define paths and point number
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
% points = [64,65,21,84,42,88,28,89,85,13,35,36,14,15,57,58,19,87,6,...
%     7,33,34,26,27,40,61,47,48,54,55,67,68,69,70,71,72,73,74,75,76];
points =  [7,13,21,28,33,42,47,85,88,89];

%% Calculate pairwise shared borders
for p=1:length(points)
    point=points(p);
    disp(['point',num2str(point)]);
    pointNumber=point;
    % load data
    load([path,'/Point',num2str(pointNumber),'/cellData3px.mat']);
    load([path,'/Point',num2str(pointNumber),'/segmentationParams3px.mat']);
    % create cell x cell matrix for border shared stat
    cells = unique(newLmod(newLmod~=0)); %all cell labels
    neighborMatrix = zeros(length(cells),length(cells)); 
    % for each cell determine its neighbors
    for i = 1:length(cells)
        % define cell object including the boundary
        cell = cells(i); 
        obj1 = imdilate((newLmod==cell),ones(3,3)); %include cell + border pixels
        % get the neighbors
        se = ones(3);   % 8-connectivity for neighbors - could be changed
        neighbors = imdilate(obj1, se); % get neighboring objects including self
        neighborLabels = unique(newLmod(neighbors)); %neighbor labels
        neighborLabels(neighborLabels==cell)=[]; %remove self
        neighborLabels(~ismember(neighborLabels,cells)) = []; %remove non-cell objects
        neighborMatrix(cell,neighborLabels)=1; %indicate cell and neighbors
    end 
    save([path,'/Point',num2str(point),'/cellNeighbors.mat'],'neighborMatrix')
end