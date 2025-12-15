%% path
addpath('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/scripts')  

%% Read tables
updatedFilePath = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/Ducts_table.csv';
T_ducts =  readtable(updatedFilePath,'VariableNamingRule' , 'preserve');

DCIS_Metadata  =  readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/DCIS_Metadata.csv','VariableNamingRule' , 'preserve');

updatedFilePath = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/singleCellTable.csv';
singleCellTable_new =  readtable(updatedFilePath,'VariableNamingRule' , 'preserve');


%% definition of tresholds
r_deep       = [0 0.85];
r_shallow    = [0.85 1];
d_periductal = 100;

%% Go over cell table and assign compartment
numRows = height(singleCellTable_new); 
cellArray = repmat({'Duct uninvolved'}, numRows, 1);
singleCellTable_new.compartment = cellArray;

ind_deep       = find(singleCellTable_new.DuctNumber>0 & singleCellTable_new.r>=r_deep(1) & singleCellTable_new.r<=r_deep(2));
ind_shallow    = find(singleCellTable_new.DuctNumber>0 & singleCellTable_new.r>=r_shallow(1) & singleCellTable_new.r<=r_shallow(2));
ind_periductal = find(singleCellTable_new.DuctNumber==0 & singleCellTable_new.distance_to_duct_edge<=d_periductal);

singleCellTable_new.compartment(ind_deep) = {'Deep infiltrating'};
singleCellTable_new.compartment(ind_shallow) = {'Shallow infiltrating'};
singleCellTable_new.compartment(ind_periductal) = {'Periductal'};

%% get immune cells
cell_meta_clusters   = unique(singleCellTable_new.cell_meta_cluster);
immune_meta_clusters = cell_meta_clusters([1 10:16 18 20:23]);
subTab_immune = singleCellTable_new(find(ismember(singleCellTable_new.cell_meta_cluster,immune_meta_clusters)),[1 2 5 311 355:358 359 360 size(singleCellTable_new,2)]);

%% sum up per fov per cell type stats for immune cells
tr_too_few = 3;
G = groupsummary(subTab_immune,{'fov','cell_meta_cluster','compartment'},"IncludeEmptyGroups",true);
fov_sum = groupsummary(subTab_immune,{'fov','cell_meta_cluster'});
fov_sum.Properties.VariableNames{3} = 'N cells';

G = innerjoin(G,fov_sum,'Keys',{'fov','cell_meta_cluster'});
G.compartment_freq = G.GroupCount./G.("N cells");
ind_too_few = find(G.("N cells")<tr_too_few);
G.compartment_freq(ind_too_few) = nan;


fov            = unique(G.fov);
clusters       = unique(G.cell_meta_cluster);
compartments   = unique(G.compartment);
immune_loc_tab = table(fov);

for f = 1:numel(fov)
    cur_fov = immune_loc_tab.fov{f};
    for c=1:numel(clusters)
        cur_cluster = clusters{c};
        for i=1:numel(compartments)
            cur_compartment = compartments{i};
            ind_G = find(strcmp(G.fov,cur_fov) & strcmp(G.cell_meta_cluster,cur_cluster)  & strcmp(G.compartment,cur_compartment));
            if ~isempty(ind_G)
                cur_colName = [cur_compartment '.' cur_cluster];
                immune_loc_tab.(cur_colName)(f) = G.compartment_freq(ind_G);
            else
                immune_loc_tab.(cur_colName)(f) = nan;
            end
        end
        
    end
    
end

%% Repeat with less granular clusters
tr_too_few = 3;
% substitute clusters fro less granular ones in subTab immune
clust_levels = table(immune_meta_clusters);
clust_levels.immune_groups =  repmat({'e'}, size(clust_levels,1), 1);
clust_levels.immune_groups{1} = 'APC';
clust_levels.immune_groups{2} = 'Immune_other';
clust_levels.immune_groups{3} = 'Mast_cells';
clust_levels.immune_groups{9} = clust_levels.immune_meta_clusters{9};
clust_levels.immune_groups([4:8]) = {'Mono_Mac_all'};
clust_levels.immune_groups([10:13]) = {'Tcell_all'};
clust_levels.immune_broad = clust_levels.immune_groups;
clust_levels.immune_broad([3:9]) = {'myeloid'};
clust_levels.immune_broad([10:13]) = {'lymphoid'};
clust_levels.immune_bin = clust_levels.immune_groups;
clust_levels.immune_bin([1:end]) = {'Immune'};

for l=1:size(clust_levels,2)
    % create immune cell table with current granularity
    subTab_immune_l = subTab_immune;
    if l>1
        for c=1:size(clust_levels,1)
            ind_cur = find(strcmp(subTab_immune_l.cell_meta_cluster,clust_levels.immune_meta_clusters{c}));
            subTab_immune_l.cell_meta_cluster(ind_cur) = clust_levels{c,l};
        end
    end
    % charecterize localization per cluster
    G = groupsummary(subTab_immune_l,{'fov','cell_meta_cluster','compartment'},"IncludeEmptyGroups",true);
    fov_sum = groupsummary(subTab_immune_l,{'fov','cell_meta_cluster'});
    fov_sum.Properties.VariableNames{3} = 'N cells';
    
    G = innerjoin(G,fov_sum,'Keys',{'fov','cell_meta_cluster'});
    G.compartment_freq = G.GroupCount./G.("N cells");
    ind_too_few = find(G.("N cells")<tr_too_few);
    G.compartment_freq(ind_too_few) = nan;
    
    
    fov            = unique(G.fov);
    clusters       = unique(G.cell_meta_cluster);
    compartments   = unique(G.compartment);
    immune_loc_tab = table(fov);
    
    for f = 1:numel(fov)
        cur_fov = immune_loc_tab.fov{f};
        for c=1:numel(clusters)
            cur_cluster = clusters{c};
            for i=1:numel(compartments)
                cur_compartment = compartments{i};
                ind_G = find(strcmp(G.fov,cur_fov) & strcmp(G.cell_meta_cluster,cur_cluster)  & strcmp(G.compartment,cur_compartment));
                if ~isempty(ind_G)
                    cur_colName = [cur_compartment '.' cur_cluster];
                    immune_loc_tab.(cur_colName)(f) = G.compartment_freq(ind_G);
                else
                    immune_loc_tab.(cur_colName)(f) = nan;
                end
            end
            
        end
        
    end
    eval([ 'immune_loc_tab_' num2str(l) '=immune_loc_tab;'])
end

%% combine for different granularity levels
keys_1  = intersect(immune_loc_tab_1.Properties.VariableNames, immune_loc_tab_2.Properties.VariableNames);
immune_loc_tab  = innerjoin(immune_loc_tab_1, immune_loc_tab_2, 'Keys', keys_1);

keys_2  = intersect(immune_loc_tab.Properties.VariableNames, immune_loc_tab_3.Properties.VariableNames);
immune_loc_tab  = innerjoin(immune_loc_tab, immune_loc_tab_3, 'Keys', keys_2);

keys_3  = intersect(immune_loc_tab.Properties.VariableNames, immune_loc_tab_4.Properties.VariableNames);
immune_loc_tab  = innerjoin(immune_loc_tab, immune_loc_tab_4, 'Keys', keys_3);

%% Save immune loc tab
tbl_name = 'immune_loc_tab';
move_old_table_to_old_folder(tbl_name,'/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/');
updatedFilePath = ['/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/' tbl_name '.csv'];
writetable(immune_loc_tab, updatedFilePath);
