%% MIBIgetMinDistFromCellCentroidToObjectPerimeter.m
% Author: Erin McCaffrey
% Modified by: (Insert your name here).

% define path and points
pathMask = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets/PAH data/Masks/';
% should be path to segmentation masks for points
path_segment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets/PAH data/segmentation_data/';
% should be path to concatenated csv, needs column with point ID and cell
% label
path_data ='/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets/';
points = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,90]; %cohort data to analyze
dataAll=dataset('File',[path_data,'/celldata_region_annotated.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
% converting to matrix for easy indexing
dataAllCell=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllCell(2:82135,[1,41]));

%define column with patient Idx and cell label for filetering
patientIdx = 1; %column with patient label
cellLabelIdx = 2; %column with cell label corresponding to segmentation mask

%create matrix with just patient IDs and cell labels for points with masks
dataAllPatientAndCells=dataAllMat(:,[patientIdx,cellLabelIdx]);
dataAllPatientAndCells= dataAllPatientAndCells(ismember(dataAllPatientAndCells(:,1), points), :);

%create output matrix
cell_perim_data = zeros(size(dataAllPatientAndCells,1),size(dataAllPatientAndCells,2)+1);
cell_perim_data(:,1)=dataAllPatientAndCells(:,1);
cell_perim_data(:,2)=dataAllPatientAndCells(:,2);

%counter for indexing summary matrix
count=0;

for p=1:length(points)
    point=points(p);
    disp(['point',num2str(point)]);
    
    % load data for mask
    load([pathMask,'/Point',num2str(point),'/regionmasks.mat']);
    
    % load segmentation label mask
    newLmod = imread([path_segment,'/Point',num2str(point),'_labels.tiff']);
    
    % Get mask perimeter
    vessel_perimeter = bwboundaries(filled_mask);
    perimeter_coordinates = cell2mat(vessel_perimeter);
    x = perimeter_coordinates(:, 2);
    y = perimeter_coordinates(:, 1);
    
    % get data just for current point
    patientInds = dataAllPatientAndCells(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllPatientAndCells(patientInds,:); %get just the data for the current point
    stats_cells = regionprops(newLmod,'Centroid'); %get spatial information for the current point
    
    % iterate through all cells and get their location based on pixel idxs
    cells = patientData(:,2); 
    for j = 1:length(cells)
        % define cell object including the boundary
        cell = cells(j);
        cell_x = stats_cells(cell).Centroid(1);
        cell_y = stats_cells(cell).Centroid(2);
        % get distance from centroid to all perimeter points
        distances = sqrt((x-cell_x).^2 + (y-cell_y).^2);
        % get min distance
        min_dist = min(distances);
        % append to output matrix
        cell_perim_data(j+count,3) = min_dist;
    end
    %update counter
    count=count+length(cells);
end

colLabels = {'Point_num','label','MinDist'};
TEXT.PnS = colLabels;
csvwrite_with_headers([pathMask,'/cell_vessel_distance.csv'],cell_perim_data,TEXT.PnS)