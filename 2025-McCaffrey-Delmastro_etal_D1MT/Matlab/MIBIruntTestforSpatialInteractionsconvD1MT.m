% MIBIruntTestforSpatialInteractionsconvD1MT.m
% Author: Erin McCaffrey
% Script reads in the z-scores for all pairwise protein interactions and
% cell type interactions. It produces concatenated matrices of each for the
% con and D1MT skewed groups separately. Next for each interaction it runs a
% t-test between the two groups and stores the p-value in a summary matrix
% (one for protein interactions and one for cell type interactions). 

%% Set up paths and define parameters 

% Define paths and point info
pathSpatial = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
resultsDir = [pathSpatial,'/dataPerCell/plots/spatial/convD1MTpooled/'];
mkdir(resultsDir);

D1MT_points = [8,9,13,14,15,16,18,25,27,29,30,31,39,40,41]; %D1MT skewed
con_points = [1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38]; %con skewed

%number of distinct markers and phenotypes
markerNum = 30;
phenoNum = 13;

%create empty matrices to store p-values
pvalues_marker = zeros(markerNum,markerNum,2);
pvalues_pheno = zeros(phenoNum,phenoNum,2);

%% Compare pairwise protein interactions

% create single matrix of z-scores for D1MT skewed

zAll_D1MT_marker = zeros(markerNum,markerNum,length(D1MT_points)); 

for i = 1:length(D1MT_points)
    point = D1MT_points(i);
    load([pathSpatial,'/Point',num2str(point),'/spatialAnalysis_100px.mat']);
    % add z-score
    zAll_D1MT_marker(:,:,i) = z;
end

% create single matrix of z-scores for con skewed

zAll_con_marker = zeros(markerNum,markerNum,length(con_points)); 

for i = 1:length(con_points)
    point = con_points(i);
    load([pathSpatial,'/Point',num2str(point),'/spatialAnalysis_100px.mat']);
    % add z-score
    zAll_con_marker(:,:,i) = z;
end

% go through each row and col of the combined z-score matrices and run
% t-test between the stack of z-scores at each position, store p value in
% the summary matrix and whether or not it was significant (p<0.05)

for i = 1:markerNum
    for j = 1:markerNum
        D1MT_zscores = reshape(zAll_D1MT_marker(i,j,:),1,length(D1MT_points));
        con_zscores = reshape(zAll_con_marker(i,j,:),1,length(con_points));
        [h,p] = ttest2(D1MT_zscores,con_zscores);
        pvalues_marker(i,j,1) = p;
        pvalues_marker(i,j,2) = h;
    end
end

%% Compare pairwise cell type interactions

% create single matrix of z-scores for D1MT skewed

zAll_D1MT_pheno = zeros(phenoNum,phenoNum,length(D1MT_points)); 

for i = 1:length(D1MT_points)
    point = D1MT_points(i);
    load([pathSpatial,'/Point',num2str(point),'/spatialAnalysisCellPheno.mat']);
    % add z-score
    zAll_D1MT_pheno(:,:,i) = z;
end

% create single matrix of z-scores for con skewed

zAll_con_pheno = zeros(phenoNum,phenoNum,length(con_points)); 

for i = 1:length(con_points)
    point = con_points(i);
    load([pathSpatial,'/Point',num2str(point),'/spatialAnalysisCellPheno.mat']);
    % add z-score
    zAll_con_pheno(:,:,i) = z;
end

% go through each row and col of the combined z-score matrices and run
% t-test between the stack of z-scores at each position, store p value in
% the summary matrix and whether or not it was significant (p<0.05)

for i = 1:phenoNum
    for j = 1:phenoNum
        D1MT_zscores = reshape(zAll_D1MT_pheno(i,j,:),1,length(D1MT_points));
        con_zscores = reshape(zAll_con_pheno(i,j,:),1,length(con_points));
        [h,p] = ttest2(D1MT_zscores,con_zscores);
        pvalues_pheno(i,j,1) = p;
        pvalues_pheno(i,j,2) = h;
    end
end

%% Create adjusted p-values with multiple comparison correction (FDR)

[h_marker, crit_p_marker, adj_ci_cvrg_marker, adj_p_marker]=fdr_bh(pvalues_marker(:,:,1),.1,'pdep','yes');
[h_pheno, crit_p_pheno, adj_ci_cvrg_pheno, adj_p_pheno]=fdr_bh(pvalues_pheno(:,:,1),.05,'pdep','yes');
       
%% Export matrices for visualization and downstream analysis

%pairwise protein csv
TEXT.PnS = markerTitles;
csvwrite_with_headers([resultsDir,'/t_test_marker.csv'],pvalues_marker(:,:,1),TEXT.PnS)

%pairwise pheno csv
TEXT.PnS = phenoTitles;
csvwrite_with_headers([resultsDir,'/t_test_pheno.csv'],pvalues_pheno(:,:,1),TEXT.PnS)

%Matlab variables
save([resultsDir,'/ttest_D1MTvcon.mat'],'pvalues_marker','pvalues_pheno',...
    'phenoTitles','markerTitles', 'zAll_con_pheno','zAll_D1MT_pheno',...
    'zAll_con_marker','zAll_D1MT_marker','adjpvalues_marker','adjpvalues_pheno');

