% % MibiQuantifyImmuneTumorMixingForCohort
% Author: Leeat Keren (with modifications made by Erin McCaffrey)
% Quantify using close cells. The myeloid-lymphocyte mixing score is
% defined as the number of myeloid-lymphoid interactions divided by the
% total myeloid-myeloid interactions.

path = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
pathSegment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
points = [1,2,8,9,10,11,12,13,14,15,16,18,19,20,21,22,...
    23,24,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41]; %cohort data
dataAll=dataset('File',[pathSegment,'/dataPerCell/NHP_cohort_data_norm_annotated_matlab.csv'],'Delimiter',','); %data with lin numerically coded
% code myeloid=1, lymph=2, granulocyte=3, imm_other=4, non_immune=5
%converting to matrix. This is not a clean way to do it but oh well for
%now..
dataAllMat=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(2:70574,[1,2,48])); %keep just the sampleID, cell label, and lineage code
sampleIdx=1; %column index of sampleID
cellLabelIdx=2; %column index of cell label
linIdx=3; %column index of lineage code


%vectors to save results
totalInteractions = zeros(length(points),2); %total cell interactions among myeloid cells and lymphocytes
myeloidLymphoInteractions = zeros(length(points),2); %total cell interactions between myeloid and lymphocytes
totalInteractions(:,1)=points;
myeloidLymphoInteractions(:,1)=points;

for i=1:length(points)
    point=points(i);
    disp(point);
    load([pathSegment,'/Point',num2str(point),'/cellNeighbors.mat']);
    % reduce neighbor matrix only to cells classified as myeloid
    currLabel = (dataAllMat(:,sampleIdx) == point); %get indices of data for this point
    currLabelIdentity = dataAllMat(currLabel,:); %create subsetted matrix for just the patient
    currInds1 = currLabelIdentity((currLabelIdentity(:,linIdx)==6 | currLabelIdentity(:,linIdx)==8),cellLabelIdx); %get just myeloid cells
    neighborMatCellsMyeloid = neighborMatrix(currInds1,currInds1); %neighbor matrix of just myeloid cells
    % Get total amount of interactions (neighborMatrix=1)
    totalInteractions(i,2) = sum(sum(neighborMatCellsMyeloid)); %sum of all myeloid-myeloid interactions
    % Get interactions between myeloid and lymphoid cells
    currInds2 = currLabelIdentity((currLabelIdentity(:,linIdx)==2),cellLabelIdx); %get lymphocyte indices
    % get myeloid and lymphoid inds to get mixed interactions
    neighborMatCellsMix = neighborMatrix(currInds1,currInds2); %get the matrix of myeloid vs lymphocytes
    myeloidLymphoInteractions(i,2) = sum(sum(neighborMatCellsMix));
end

percentMix = myeloidLymphoInteractions./totalInteractions;
percentMix(:,1)=points;

% sort patients according to interactions percentage
percentMixS = sortrows(percentMix,2);
figure;
bar(percentMixS(:,2));
set(gca,'XTick',(1:length(points)),'XTickLabel',percentMixS(:,1));

resultsDir = [pathSegment,'/dataPerCell_3px/plots/draft_figs/spatial/mixing/'];
mkdir(resultsDir);
saveas(gcf,([resultsDir,'mixing_score_bar']),'epsc')
csvwrite_with_headers([pathSegment,'/dataPerCell_3px/mixingscore.csv'],percentMix,{'SampleID','Mixing_Score'})
