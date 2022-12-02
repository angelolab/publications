%% MIBIannotateCellsinGranulomaZones.m
% Author: Erin McCaffrey
% Modified by: (Insert your name here)
% This script takes in a mask of any kind of zones. For
% each point and cell it creates a binary matrix of whether or not that
% cell label is in the mask area or on the border of the mask.

% define path and points
pathMask = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets/PAH data/Masks/';
% should be path to segmentation masks for points
path_segment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets/PAH data/segmentation_data';
% should be path to concatenated csv, needs column with point ID and cell
% label
path_data ='/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets/';
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
n_classes = 2; % number of classifications/regions (ie. mask and periphery)
cell_mask_data = zeros(size(dataAllPatientAndCells,1),size(dataAllPatientAndCells,2)+n_classes);
cell_mask_data(:,1)=dataAllPatientAndCells(:,1);
cell_mask_data(:,2)=dataAllPatientAndCells(:,2);
maskIdx = 3;
periphIdx = 4;

%counter for indexing summary matrix
count=0;

for p=1:length(points)
    point=points(p);
    disp(['point',num2str(point)]);
    
    % load data for mask
    load([pathMask,'/Point',num2str(point),'/regionmasks.mat']);
    
    % load segmentation label mask
    newLmod = imread([path_segment,'/Point',num2str(point),'_labels.tiff']);
    
    % combine masks
    filled_mask(filled_mask>0)=1;
    periph_mask(periph_mask>0)=2;
    combined_mask=filled_mask+periph_mask;
    
    % get region stats for mask
    stats_mask = regionprops(combined_mask,'PixelIdxList');
    core_px = stats_mask(1).PixelIdxList;
    periph_px = stats_mask(2).PixelIdxList;
    
    % get data just for current point
    patientInds = dataAllPatientAndCells(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllPatientAndCells(patientInds,:); %get just the data for the current point
    stats_cells = regionprops(newLmod,'PixelIdxList'); %get spatial information for the current point
    
    % iterate through all cells and get their location based on pixel idxs
    cells = patientData(:,2); 
    for j = 1:length(cells)
        % define cell object including the boundary
        cell = cells(j);
        cell_px = stats_cells(cell).PixelIdxList;
        %check IM
        mask_overlap = intersect(cell_px,core_px);
        if (~isempty(mask_overlap))
            cell_mask_data(j+count,maskIdx) = 1;
        end
        %check peripheral zone
        periph_overlap = intersect(cell_px,periph_px);
        if (~isempty(periph_overlap))
            cell_mask_data(j+count,periphIdx) = 1;
        end
    end
    %update counter
    count=count+length(cells);
end

colLabels = {'SampleID','label','In_Mask','In_Mask_Periph'};
TEXT.PnS = colLabels;
csvwrite_with_headers([pathMask,'/cell_vessel_distance.csv'],cell_mask_data,TEXT.PnS)