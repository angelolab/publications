%% MIBIannotateCellsinGranulomaZones.m
% Author: Erin McCaffrey
% Modified by: (Insert your name here)
% This script takes in a mask of any kind of zones. For
% each point and cell it creates a binary matrix of whether or not that
% cell label is in the mask area or on the border of the mask.

% define path and points
pathMask = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
% should be path to segmentation masks for points
path_segment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
% should be path to concatenated csv, needs column with point ID and cell
% label
path_data ='/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell/';
points = [8,9,13,14,15,16,18,25,27,29,30,31,39,40,41,1,2,10,11,12,19,...
    20,21,23,24,28,32,33,34,35,36,37,38]; %cohort data to analyze
dataAll=dataset('File',[path_data,'/NHP_cohort_data_norm_annotated_matlab.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
% converting to matrix for easy indexing
dataAllCell=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllCell(2:70574,[1,2]));

%define column with patient Idx and cell label for filetering
patientIdx = 1; %column with patient label
cellLabelIdx = 2; %column with cell label corresponding to segmentation mask

%create matrix with just patient IDs and cell labels for points with masks
dataAllPatientAndCells=dataAllMat(:,[patientIdx,cellLabelIdx]);
dataAllPatientAndCells= dataAllPatientAndCells(ismember(dataAllPatientAndCells(:,1), points), :);
% sort the points and matrix to ensure correct assignment
points = sort(points);
dataAllPatientAndCells = sortrows(dataAllPatientAndCells, patientIdx);

%create output matrix
n_classes = 1; % number of classifications/regions (ie. mask and periphery)
cell_mask_data = zeros(size(dataAllPatientAndCells,1),size(dataAllPatientAndCells,2)+n_classes);
cell_mask_data(:,1)=dataAllPatientAndCells(:,1);
cell_mask_data(:,2)=dataAllPatientAndCells(:,2);
maskIdx = 3;

%counter for indexing summary matrix
count=0;

for p=1:length(points)
    point=points(p);
    disp(['point',num2str(point)]);
    
    % load data for mask
    mask = imread([pathMask,'Point',num2str(point),'/myeloid_mask.tif']);
    mask(mask>0) = 1;
    
    % load segmentation label mask
    load([path_segment,'/Point',num2str(point),'/segmentationParams.mat']);
   
    % get region stats for mask
    stats_mask = regionprops(mask,'PixelIdxList');
    if(~isempty(stats_mask))
        core_px = stats_mask.PixelIdxList;
    else
        count=count+length(cells);
        break
    end
    
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
        %check myeloid zone
        mask_overlap = intersect(cell_px,core_px);
        if (~isempty(mask_overlap))
            cell_mask_data(j+count,maskIdx) = 1;
        end
    end
    %update counter
    count=count+length(cells);
end

colLabels = {'SampleID','cellLabelInImage','In_Mask'};
TEXT.PnS = colLabels;
csvwrite_with_headers([pathMask,'/dataPerCell/core-mask-data.csv'],cell_mask_data,TEXT.PnS)