%% MIBIannotateCellsinMaskZones.m
% Author: Erin McCaffrey
% This script takes in a mask of any kind of zones. For
% each point and cell it creates a binary matrix of whether or not that
% cell label is in the mask area, on the border of the mask, or in the area
% surrounding the mask

%% 1. Import data and define points
pathMask = 'path_to_zone_masks';
% should be path to segmentation masks for points
path_segment = 'path_to_segmentation_masks';
% should be path to concatenated csv w/ column for point ID and cell label
path_data = 'path_to_sc_data';
% points to analyze
points = [3108]; 

% import data
sc_data = [path_data,'sc_data.csv'];
opts = detectImportOptions(sc_data);

% convert to a matrix
dataAll = readtable(sc_data,opts);
dataAllMat=table2cell(dataAll);
% subset just the SampleID, cellLabelInImage
dataAllMat=cell2mat(dataAllMat(:,[1:2,53]));

%% 2. define column with patient Idx and cell label for filetering

patientIdx = 1; %column with patient label
cellLabelIdx = 57; %column with cell label corresponding to segmentation mask

%% 3. Create matrix with just patient IDs and cell labels for points with masks

dataAllPatientAndCells=dataAllMat(:,[patientIdx,cellLabelIdx]);
dataAllPatientAndCells = dataAllPatientAndCells(ismember(dataAllPatientAndCells(:,1), points), :);

%% 4. Create output matrix
n_classes = 2; % number of classifications/regions (ie. mask and periphery)
cell_mask_data = zeros(size(dataAllPatientAndCells,1),size(dataAllPatientAndCells,2)+n_classes);
cell_mask_data(:,1)=dataAllPatientAndCells(:,1);
cell_mask_data(:,2)=dataAllPatientAndCells(:,2);
maskIdx = 3;
periphIdx = 4;

%% 5. Run analysis for each point

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
        %check myeloid zone
        mask_overlap = intersect(cell_px,core_px);
        if (~isempty(mask_overlap))
            cell_mask_data(j+count,maskIdx) = 1;
        end
        %check periph zone
        periph_overlap = intersect(cell_px,periph_px);
        if (~isempty(periph_overlap))
            cell_mask_data(j+count,periphIdx) = 1;
        end
    end
    %update counter
    count=count+length(cells);
end

%% 6. Export data
colLabels = {'SampleID','cellLabelInImage','In_Mask','In_Mask_Periph'};
TEXT.PnS = colLabels;
csvwrite_with_headers([pathMask,'/cell_mask_annotations.csv'],cell_mask_data,TEXT.PnS)