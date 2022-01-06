%MIBIpooledSpatialOrganizationAnalysis.m
%Author: Erin McCaffrey
%Date created: 06/09/19
%This script performs different pooled analysis from the spatial protein
%enrichment analysis that was performed with Leeat's
%MIBIQuantifySpatialOrganizationByAdjacenyMatandBeostr.m, First it produces
%an average z-score matrix. Seconldy, it determines the percent images
%significant or not for each pairwise enrichment (p-value for # times
%bootstrap >= real) and for each pairwise avoidance (p-value for # times
%bootstrap<= real). Lastly based on the average z-score matrix, it
%determines whether each interaction is on average an avoidance (-z),
%enrichment (+z), or neither (0) and determines the % images + for the
%enriched or avoided interactions.

%% 1. Define paths and point info
pathSpatial = [cd,'/segmentation'];
points = [21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,...
          26,27,40,61,47,48,54,55,92,93]; 
resultsDir = [cd,'/data'];

%% 2. Define important values
numMarkers=36; %number of markers for pairwise comparison (no HH3)
p_cutoff=0.05; %p-value cutoff

%% 3. Create matrices for summary stats
allZscores=zeros(numMarkers,numMarkers,length(points)); % store z scores for all points and interactions
avgZscore=zeros(numMarkers,numMarkers); % matrix for avg z-score for all interactions
pValueAvoid=zeros(numMarkers,numMarkers); % matrix for % significant for avoidance
pValueEnrich=zeros(numMarkers,numMarkers); % matrix for % significant enrichment
pValueCombo=zeros(numMarkers,numMarkers); % matrix for summary interactions based on avg z-score

%% 4. Summarize
for i=1:length(points)
    point=points(i);
    disp(['Working on point:',num2str(point)]);
    load([pathSpatial,'/Point',num2str(point),'/spatialAnalysisv2.mat']);
    %add z-score matrix to average
    allZscores(:,:,i)=z;
    for j=1:numMarkers
        for k=1:numMarkers
            %add significant interactions/avoidances if applicable
            if(p(j,k,1)<p_cutoff && adj_p(j,k)<p_cutoff)
                pValueEnrich(j,k)=pValueEnrich(j,k)+1;
            end
            if(p(j,k,2)<p_cutoff && adj_p(j,k)<p_cutoff)
                pValueAvoid(j,k)=pValueAvoid(j,k)+1;
            end
        end
    end
end

% get average z-scores
avgZscore=nanmean(allZscores,3);

% get percent enrichment and avoided
pValueEnrich=(pValueEnrich/length(points))*100;
pValueAvoid=(pValueAvoid/length(points))*100;

%pooled enrichment + analysis with decision for interaction based on avg z
%score
for i=1:length(points)
    for j=1:numMarkers
        for k=1:numMarkers
            % interaction
            if(avgZscore(j,k)>0)
                %grab percent positive from the enriched p-value percent
                pValueCombo(j,k)=pValueEnrich(j,k);
            elseif(avgZscore(j,k)<0)
                pValueCombo(j,k)=-1*pValueAvoid(j,k);
            end
        end
    end
end

%% 5. Plot and save

% avg Z score
avgZ_plot=avgZscore;
avgZ_plot(isnan(avgZ_plot)) = 0;
avgZ_plot(isinf(avgZ_plot)) = 0;
hmap=clustergram(avgZ_plot([1:34],[1:34]),'RowLabels', markerTitles([1:34]),...
    'ColumnLabels', markerTitles([1:34]), 'Colormap', 'redbluecmap','DisplayRange',...
    5, 'DisplayRatio', 0.1);
addTitle(hmap, 'Average Z-score')

fig1=plot(hmap);

TEXT.PnS = markerTitles;
csvwrite_with_headers([pathSpatial,'/pooled-enrich_TB.csv'],avgZscore, TEXT.PnS)


% percent positive for enrichments
penrich_plot=pValueEnrich;
pe_hmap=clustergram(penrich_plot([1:34,36],[1:34,36]),'RowLabels', markerTitles([1:34,36]),...
    'ColumnLabels', markerTitles([1:34,36]), 'Colormap', 'redbluecmap','DisplayRange',...
    100, 'DisplayRatio', 0.1);
addTitle(pe_hmap, '% of Points Significant for Enrichment')

fig2=plot(pe_hmap);

% percent positive for avoidances
pavoid_plot=pValueAvoid;
pa_hmap=clustergram(pavoid_plot([1:34,36],[1:34,36]),'RowLabels', markerTitles([1:34,36]),...
    'ColumnLabels', markerTitles([1:34,36]), 'Colormap', 'redbluecmap','DisplayRange',...
    100, 'DisplayRatio', 0.1);
addTitle(pa_hmap, '% of Points Significant for Avoidance')

% percent positive for all interactions binarized
pcombo_plot=pValueCombo;
pc_hmap=clustergram(pcombo_plot([1:34,36],[1:34,36]),'RowLabels', markerTitles([1:34,36]),...
    'ColumnLabels', markerTitles([1:34,36]), 'Colormap', 'redbluecmap','DisplayRange',...
    100, 'DisplayRatio', 0.1);
addTitle(pc_hmap, '% of Points Significant for Interaction')


save([path,'/spatialAnalysisPooledTB.mat'],'allZscores','avgZscore','pValueAvoid','pValueEnrich','pValueCombo','markerTitles');