% MIBI_macrophagee_mixing_scores.m
% This script calculates the mixing score between T cells and macrophages
% in: 1) The whole image 2)


%% Set up paths
path = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2';
seg_path = [path, '/final_segmentation'];
results_path = [path, '/spatial_analysis/mixing'];

%% Read in single cell data, subset to only include sample, tiled label, and pheno
sc_data = readtable([path, '/spatial_analysis/cell_cohort_data_metabolic_zones.csv']);
sc_data_sub = sc_data(:, ["sample","tiled_label","pheno_corrected", "glyco_zone", "IDO1_zone"]);

%% Convert sample and phenotype columns to categorical vars
sc_data_sub.sample = categorical(sc_data_sub.sample);
sc_data_sub.pheno_corrected = categorical(sc_data_sub.pheno_corrected);

%% Optionally subset to create zone-specific mixing scores
sc_data_mcore = sc_data_sub(sc_data_sub.glyco_zone == 1 | sc_data_sub.IDO1_zone == 1, :);
sc_data_IDOzone = sc_data_sub(sc_data_sub.IDO1_zone == 1, :);
sc_data_glycozone = sc_data_sub(sc_data_sub.glyco_zone == 1, :);

%% Choose input
input_data = sc_data_sub;

%% Define samples to operate on
samples = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,...
    29,32,33,34,35,37,38,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,...
    58,59];

%% Define cells to generate mixing score for
macrophages = {'CD11c+_Mac', 'CD14+_Mac_Mono', 'CD14+CD11c+_Mac', 'CD163+_Mac',...
    'CD206+_Mac', 'CD68+_Mac', 'FN1+_Mac','giant_cell'};

%% Create threshold qualify as 'associated'
t = 25; % 25 px is approx 10 um

%% Create empty struct for output and counter to index
output = struct;
count = 1;

%% Generate and save mixing score per sample and per cell type
for i=1:length(samples)

    sample_num=samples(i);
    disp(['sample',num2str(sample_num)]);

    % load the segmentation and the cell data
    seg_mask = imread([seg_path,'/sample',num2str(sample_num),'_labels.tiff']);

    % get centers of cells
    stats = regionprops(seg_mask,'centroid');

    % get x,y coordinates in a mat
    centroidCoord = zeros(length(stats),2);
    for c=1:length(stats)
        centroidCoord(c,1) = stats(c).Centroid(1);
        centroidCoord(c,2) = stats(c).Centroid(2);
    end

    % subset cell table for sample
    sample_data = input_data(input_data.sample == ['sample',num2str(sample_num)], :);

    % get unique labels
    sample_labels = unique(sample_data.tiled_label);

    % get T cell labels and coords for the sample
    cd4t_labels = sample_data(sample_data.pheno_corrected == "CD4+Tcell", :).tiled_label;
    cd8t_labels = sample_data(sample_data.pheno_corrected == "CD8+Tcell", :).tiled_label;
    all_t_labels = vertcat(cd4t_labels, cd8t_labels);
%     all_t_labels = cd8t_labels;
    t_coords =  centroidCoord(all_t_labels,:); 

    % get all macrophage labels and coords for sample
    all_mac_labels = sample_data(ismember(sample_data.pheno_corrected, macrophages), :).tiled_label;
    all_mac_coords = centroidCoord(all_mac_labels,:); 

    % generate an overall t-mac mixing score
    allmacDistMat = pdist2(all_mac_coords,all_mac_coords);
    allmacDistMatBin = zeros(size(allmacDistMat));
    allmacDistMatBin(allmacDistMat<=t) = 1;
    allmacDistMatBin_nonself = allmacDistMatBin - diag(diag(allmacDistMatBin));
    all_n_self_interactions = sum(sum(allmacDistMatBin_nonself));

    totaltmacDistMat = pdist2(all_mac_coords,t_coords);
    totaltmacDistMatBin = zeros(size(totaltmacDistMat));
    totaltmacDistMatBin(totaltmacDistMat<=t) = 1;
    n_t_mac_interactions = sum(sum(totaltmacDistMatBin));

    overall_mix_score = n_t_mac_interactions / all_n_self_interactions;


    % get mixing scores per mac type
    for j=1:length(macrophages) 

        mac_pop = macrophages{j};

        % get labels and coords for this population
        mac_labels = sample_data(sample_data.pheno_corrected == mac_pop, :).tiled_label;
        mac_coords = centroidCoord(mac_labels,:);

        % determine n mac-mac interactions
        macDistMat = pdist2(mac_coords,mac_coords);
        macDistMatBin = zeros(size(macDistMat));
        macDistMatBin(macDistMat<=t) = 1;
        macDistMatBin_nonself = macDistMatBin - diag(diag(macDistMatBin));
        n_self_interactions = sum(sum(macDistMatBin_nonself));

        % determine n mac-t cell interactions
        tDistMat = pdist2(mac_coords,t_coords);
        tDistMatBin = zeros(size(tDistMat));
        tDistMatBin(tDistMat<=t) = 1;
        n_t_interactions = sum(sum(tDistMatBin));

        % calculate mixing score
        output(count).sample = strcat('sample', num2str(sample_num));
        output(count).total_mixing_score = overall_mix_score;
        output(count).mac_pop = mac_pop;
        output(count).pop_count = length(mac_labels);
        output(count).t_count = length(all_t_labels);
        output(count).n_mac_int = n_self_interactions;
        output(count).n_t_int = n_t_interactions;
        output(count).pop_mixing_score = n_t_interactions / n_self_interactions;
        
        % update counter before proceeding
        count = count + 1;
  
    end
end

%% convert results to table and export
output_table = struct2table(output);
writetable(output_table, [results_path, '/mac_mixing_score_whole_gran.csv']);