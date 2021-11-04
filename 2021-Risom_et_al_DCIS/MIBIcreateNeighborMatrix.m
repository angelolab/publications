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

%% read in cell data
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
dataAll=dataset('File',[path,'/dataPerCell_3px/granA_cellpheno_CS-asinh-norm_matlab_revised.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%converting to matrix. This is not a clean way to do it but oh well for
%now..
dataAllMat=dataset2cell(dataAll);
% subset just the SampleID, cellLabelInImage, and FlowSOMID
dataAllMat=cell2mat(dataAllMat(2:55713,[1:2,53]));

%% read in pheno data
phenoKey=dataset('File',[path,'/dataPerCell_3px/cellpheno_numkey.csv'],'Delimiter',','); 
phenoKey = sortrows(phenoKey, {'Code'}, {'ascend'});
pTitles = phenoKey(:,1); %names of the cell types
phenoTitles=dataset2cell(pTitles);
phenoTitles=phenoTitles(2:21,:);
pCodes = phenoKey(:,2); %numerical code for cell types
pCodes=dataset2cell(pCodes);
pCodes=pCodes(2:21,:);
pCodes=cell2mat(pCodes);
phenoNum = length(pCodes); %number of markers to compare

%% initiate empty matrices for cell neighborhood data
cell_neighbor_counts = zeros(size(dataAllMat,1), phenoNum+2);
cell_neighbor_freqs = zeros(size(dataAllMat,1), phenoNum+2);
% add SampleID and cell label info to the matrices
cell_neighbor_counts(:,1:2) = dataAllMat(:,1:2);
cell_neighbor_freqs(:,1:2) = dataAllMat(:,1:2);

%% Initiate distance threshold in px and define indices
t = 50;
patientIdx = 1;
cellIdx = 2;
phenoIdx = 3;

%% Build data matrix for each point
points = unique(dataAllMat(:,1)); % points in dataset
cell_count = 1; % for accurate indexing in neighbor matrices
% iterate through all points
for i=1:length(points)
    point=points(i);
    
    % load relevant data
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
        
%% export as csv
resultsPath = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px/';
channelLabels = ['SampleID';'cellLabelInImage';pTitles.Pheno];
TEXT.PnS = channelLabels;
csvwrite_with_headers([resultsPath,'/cell_neighbor_counts_50px.csv'],cell_neighbor_counts,TEXT.PnS)
csvwrite_with_headers([resultsPath,'/cell_neighbor_freqs_50px.csv'],cell_neighbor_freqs,TEXT.PnS)
        
