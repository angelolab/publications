%% load glycan data
glycan_DE_limma = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/de_glycans_TA511_PVbinned.csv');
glycan_exp      = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/IFmaldi_TA511_mean_gly_permask_PVbinned.csv','PreserveVariableNames',true);

%% Arrange DEG table 
close_comp_only = 0;
remove_EVTP     = 1;

var_names       = glycan_DE_limma.Properties.VariableNames;
inds            = find(contains(var_names,'P_Value'));
names           = var_names(inds)';

if remove_EVTP==0
    types           = {'PV','EVTA','EVTI','EVTP','EVTE'};
else
    types           = {'PV','EVTA','EVTI','EVTE'};
    ind_rem         = find(strcmp(glycan_exp.cell_anno,'EVT_P'));
    glycan_exp(ind_rem,:) = [];
end

arep     = strrep(names,'group','');

if close_comp_only==1
    bin_keep = zeros(1,numel(names));
    for i=1:numel(names)
        temp = strsplit(arep{i},'_');
        typ1 = temp{1};
        typ2 = temp{2};
        try
            ind_typ1 = find(strcmp(typ1,types));
            ind_typ2 = find(strcmp(typ2,types));
        catch
            continue
        end
        if abs(ind_typ1-ind_typ2)==1
            bin_keep(i) = 1;
        end
    end
    tab_DEG = glycan_DE_limma(:,[1 inds(bin_keep==1) inds(bin_keep==1)+1 inds(bin_keep==1)-1]);
else
    ind_no_p           = find(~contains(var_names,'EVTP'));
    if remove_EVTP==1
        tab_DEG = glycan_DE_limma(:,ind_no_p );
    else
        tab_DEG = glycan_DE_limma;
    end
end

T = tab_DEG;
sortedNames = sort(T.Properties.VariableNames(2:end));
T1 = T(:,sortedNames);
T2 = [T(:,1) T(:,sortedNames)];
tab_DEG = T2;

%% set TRs
p_tr  = 0.05;
fc_tr = 0.25;
q_tr  = 0.25;

%% get DEG
var_names       = tab_DEG.Properties.VariableNames;
inds_p          = find(contains(var_names,'P_Value'));
names           = var_names(inds_p)';
comps           = strrep(names,'_P_Value','');
DEG_all         = {};

for i=1:numel(comps)
    cur_comp = comps{i};
    p  = tab_DEG.([cur_comp '_P_Value']);
    q  = tab_DEG.([cur_comp '_adj_P_Val']);
    fc = tab_DEG.([cur_comp '_logFC']);
    
    DEG_cur = tab_DEG.Glycan(find(p<=p_tr & q<=q_tr & abs(fc)>=fc_tr));
    DEG_all = [DEG_all;DEG_cur];

end
DEG_all = unique(DEG_all);

%% Calc median by EVT type
[G,id]          = findgroups(glycan_exp.cell_anno_update);
med_by_type     = splitapply(@median,glycan_exp{:,6:end},G);
temp            = array2table(med_by_type);
typ             = table(id);
med_by_type_tab = [typ temp ];
med_by_type_tab.Properties.VariableNames(2:end) = glycan_exp.Properties.VariableNames(6:end);

T2 = rows2vars(med_by_type_tab);
T2(1,:) = [];
T2.Properties.VariableNames = ['TargetName' cellstr(id)'];

med_by_type_tab = T2(:,[1 5 2 4 3]);

%% Get DEG expression data
indin      = get_gene_inds_of_list1_in_list2(DEG_all,med_by_type_tab.TargetName);
exp_dat    = med_by_type_tab(indin,:);

inter_traj = cell2mat(exp_dat{:,2:end});%[exp_dat.FV ,exp_dat.VCT, exp_dat.EVT_A, exp_dat.EVT_I, exp_dat.EVT_E];
gene_name  = exp_dat.TargetName;

%% calculate COM for trending genes
max_vec   = max(abs(inter_traj),[],2);
trend_mat = inter_traj./repmat(max_vec,1,size(inter_traj,2));
com       = sum(trend_mat.*repmat(1:size(trend_mat,2),size(trend_mat,1),1),2)./sum(trend_mat,2);

%% normalize for comparability between genes
inter_traj(inter_traj<0)=0;
norm_flag = 2; 

if norm_flag==2
    mean_by_stage_norm = normalize(inter_traj,2);
elseif norm_flag==1
    % max normalize
    max_vec   = max(abs(inter_traj),[],2);
    mean_by_stage_norm = inter_traj./repmat(max_vec,1,size(inter_traj,2));
end
mean_by_stage_table_norm = array2table(mean_by_stage_norm);
mean_by_stage_table_norm = [gene_name mean_by_stage_table_norm];
mean_by_stage_table_norm.Properties.VariableNames{1} = 'gene_name';
mean_by_stage_table_norm.com = com;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% graph based approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load glycan type grouping 
% gly_anno = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/glycan_peaklist_typeAnnotations_051524.csv');

gly_anno = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/glycan_peaklist_typeAnnotations_nopartialgal_091525.csv');

%% read gene by gly type annotation
% read annotation file

% file_path = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/glycan_types_enzymes_051524.xlsx';
% sheet_name = '051524_gene_types';

file_path = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/glycan_types_enzymes_nopartialgal_091525.xlsx';
sheet_name = 'glycan_types_enzymes_nopartialg';

gene_by_gly_type = readtable(file_path,'Sheet' ,sheet_name);

% read DE Genes table
EVT_DEG_by_type = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/EVT_DEG_heatmap_DESeq2_updated.csv');
ind_DE_enzymes = get_gene_inds_of_list1_in_list2(gene_by_gly_type.Var1,EVT_DEG_by_type.gene_name);
DE_enzymes     = EVT_DEG_by_type(ind_DE_enzymes(ind_DE_enzymes>0),:);

%% K clustering
% concat DE genes and glycans for k-clustering
% DE_gly_enz = [DE_enzymes]; % do only enzymes
DE_gly_enz = [DE_enzymes(:,1:end-1);mean_by_stage_table_norm];

%re-do k clustering 
N_Ks = 7;
mean_by_stage_table_norm_gly_enz = k_heatmap(N_Ks,DE_gly_enz,types,1,1);

%% create graph nodes
DE_gly_nodes   = mean_by_stage_table_norm.gene_name;
k_nodes        = arrayfun(@(x) num2str(x), (1:N_Ks)', 'UniformOutput', false);
enz_nodes      = DE_enzymes.gene_name;
gly_type_nodes = gene_by_gly_type.Properties.VariableNames(2:end)'; 
node_names     = [DE_gly_nodes;k_nodes ;enz_nodes;gly_type_nodes];

%% create adj mat
adj_mat = zeros(numel(node_names),numel(node_names)); 

% do gly->k edges 
for i=1:numel(DE_gly_nodes)
    cur_name         = DE_gly_nodes{i};
    cur_k            = mean_by_stage_table_norm_gly_enz.k_cluster(find(strcmpi(cur_name,mean_by_stage_table_norm_gly_enz.gene_name)));
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_k_node_ind   = find(strcmpi(num2str(cur_k),node_names));
    adj_mat(cur_name_node_ind,cur_k_node_ind)        = 1;
    adj_mat(cur_k_node_ind,cur_name_node_ind)        = 1;
end

% do enz->k edges 
for i=1:numel(enz_nodes)
    cur_name          = enz_nodes{i};
    cur_k             = mean_by_stage_table_norm_gly_enz.k_cluster(find(strcmpi(cur_name,mean_by_stage_table_norm_gly_enz.gene_name)));
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_k_node_ind    = find(strcmpi(num2str(cur_k),node_names));
    adj_mat(cur_name_node_ind,cur_k_node_ind)        = 1;
    adj_mat(cur_k_node_ind,cur_name_node_ind)        = 1;
end

% do gly->types edges
offset = 2;
for i=1:numel(DE_gly_nodes)
    cur_name          = DE_gly_nodes{i};
    cur_ind           = find(strcmpi(cur_name,gly_anno.composition));
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_types         = gly_anno.Properties.VariableNames(find(gly_anno{cur_ind,1+offset:end}==1)+offset);
    for j=1:numel(cur_types)
        cur_type_ind    = find(strcmpi(cur_types{j},node_names));
        adj_mat(cur_name_node_ind,cur_type_ind)        = 1;
        adj_mat(cur_type_ind,cur_name_node_ind)        = 1;
    end
end

% do enz->types edges
offset = 1;
for i=1:numel(enz_nodes)
    cur_name          = enz_nodes{i};
    cur_ind           = find(strcmpi(cur_name,gene_by_gly_type.Var1));
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_types         = gene_by_gly_type.Properties.VariableNames(find(gene_by_gly_type{cur_ind,1+offset:end}==1)+offset);
    for j=1:numel(cur_types)
        cur_type_ind    = find(strcmpi(cur_types{j},node_names));
        adj_mat(cur_name_node_ind,cur_type_ind)        = 1;
        adj_mat(cur_type_ind,cur_name_node_ind)        = 1;
    end
end


%% create graph
% Create a graph from the adjacency matrix
G = graph(adj_mat);
G.Nodes.Name = node_names;

% % Plot the graph
% figure
% plot(G, 'NodeLabel', G.Nodes.Name);

%% find cycles and annotate them by type and k clust
[cyc_tab,cyc_per_type_tab] = cyc_summary(G,DE_gly_nodes,k_nodes ,enz_nodes,gly_type_nodes);
cyc_per_type_tab

%% Generate random graphs as Null distribution
N_rand = 500;
rand_cycs_tab = table(gly_type_nodes);
for i=1:N_rand
    G_rand                               = generate_random_adj_mat(DE_gly_nodes,k_nodes ,enz_nodes,gly_type_nodes,node_names,adj_mat);
    [cyc_tab_rand,cyc_per_type_tab_rand] = cyc_summary(G_rand,DE_gly_nodes,k_nodes ,enz_nodes,gly_type_nodes);
    eval(['rand_cycs_tab.rand_run_' num2str(i) ' =  cyc_per_type_tab_rand{:,2};']);
end

%% Z score and plot
z_mat = normalize([rand_cycs_tab{:,2:end} cyc_per_type_tab.cyc_per_type],2);
cyc_per_type_tab.Z_score_cyc_num = z_mat(:,end);

figure;
bar(cyc_per_type_tab.Z_score_cyc_num);
xticklabels(gly_type_nodes)
ylabel('Z-score')

%% Save results
run_name = ['nopartialgal_type_list_K_' num2str(N_Ks)];
writetable(cyc_per_type_tab,['Z_score_tab_' run_name '.csv']);
writetable(mean_by_stage_table_norm_gly_enz,['combined_gly_enz_heatmap_' run_name '.csv']);



