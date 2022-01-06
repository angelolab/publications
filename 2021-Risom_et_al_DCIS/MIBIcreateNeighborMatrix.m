%% MIBIcreateNeighborMatrix.m
% Author: Erin McCaffrey
% Date created: 12/01/2019
% Script reads in the annotated data matrix for cell data for the whole
% cohort. For each point it reads in the distance matrix and produces a
% matrix of all cells x all phenotypes where it counts the number of cells
% with a given distance threshold (t) that are each phenotype. It creates
% an additional matrix with the frequency (normalized by all cells within
% the threshold) that belong to each phenotype. It exports the resulting
% matrix as a csv and saves it as a variable environment.

%% 1. Read in cell data

path = 'path_to_data';

% import data
sc_data = [path,'sc_data.csv'];
opts = detectImportOptions(sc_data);

% convert to a matrix
dataAll = readtable(sc_data,opts);
dataAllMat=table2cell(dataAll);
% subset just the SampleID, cellLabelInImage, and FlowSOMID
dataAllMat=cell2mat(dataAllMat(:,[1:2,53]));

%% 2. Read in pheno data

phenoKey=readtable([path,'/pheno_key.csv']); 

% sort in descending order of pheno keys
phenoKey = sortrows(phenoKey, {'Code'}, {'ascend'});

% pull out phenotype names and codes
pTitles = phenoKey(:,1); %names of the cell types
phenoTitles=table2cell(pTitles);
pCodes = phenoKey(:,2); %numerical code for cell types
pCodes=table2cell(pCodes);
pCodes=cell2mat(pCodes);
phenoNum = length(pCodes); %number of markers to compare

%% 3. initiate empty matrices for cell neighborhood data

cell_neighbor_counts = zeros(size(dataAllMat,1), phenoNum+2);
cell_neighbor_freqs = zeros(size(dataAllMat,1), phenoNum+2);
% add SampleID and cell label info to the matrices
cell_neighbor_counts(:,1:2) = dataAllMat(:,1:2);
cell_neighbor_freqs(:,1:2) = dataAllMat(:,1:2);

%% 4. Initiate distance threshold in px and define indices

t = 50;
patientIdx = 1;
cellIdx = 2;
phenoIdx = 3;

%% 5. Build data matrix for each point

points = unique(dataAllMat(:,1)); % points in dataset
cell_count = 1; % for accurate indexing in neighbor matrices

% iterate through all points
for i=1:length(points)
    point=points(i);
    
    % load relevant data (distance matrix)
    disp(['Working on point:',num2str(point)]);
    load([path,'/Point',num2str(point),'/cellDistances.mat']);
    
    % subset the 
    patientInds = dataAllMat(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllMat(patientInds,:); %get just the data for the current point
    
    % iterate through all cells
    for j = 1:size(patientData,1)
        % get cell
        cell = patientData(j, cellIdx);
        % get all cell neighbors within threshold distance
        cellDistMat = distancesMat(cell,:);
        cellDistMatBin = zeros(size(cellDistMat));
        cellDistMatBin(cellDistMat<t) = 1;
        % get indices (labels) of close cells
        neighbor_labels = find(cellDistMatBin);
        % filter out non-cellular objects
        neighbor_labels_cells = intersect(neighbor_labels, patientData(:,cellIdx));
        % count phenotypes in cell neighbors
        count_vec = zeros(1,phenoNum);
        neighborInds = ismember(patientData(:,cellIdx), neighbor_labels_cells(:,1));
        pheno_vec = patientData(neighborInds, phenoIdx);
        for k = 1:phenoNum
            count_vec(1,k) = sum(pheno_vec==k);
        end
        % add to neighborhood matrices
        cell_neighbor_counts(cell_count,3:22) = count_vec(1,:);
        cell_neighbor_freqs(cell_count,3:22) = (count_vec(1,:)/length(neighbor_labels_cells));
        % update cell count
        cell_count = cell_count+1;
    end 
end
        
%% 6. Export as csv

resultsPath = [path, '/results'];
channelLabels = ['SampleID';'cellLabelInImage';pTitles.Pheno];
TEXT.PnS = channelLabels;
csvwrite_with_headers([resultsPath,'/cell_neighbor_counts_50px.csv'],cell_neighbor_counts,TEXT.PnS)
csvwrite_with_headers([resultsPath,'/cell_neighbor_freqs_50px.csv'],cell_neighbor_freqs,TEXT.PnS)
        
