%% load Nanostring data
var_list = {'roi_type_EVT','counts_norm_tab_DATA_cl','DATA_cl'};
save_path = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/';

for i = 1:length(var_list)
    var_name = var_list{i};
    filename = fullfile(save_path, [var_name '.csv']);
    naming = ['VariableNamingRule' ,',', 'preserve'];
    eval([var_name ' = readtable(filename, ' naming ');']);
end
roi_type_EVT = sortrows(roi_type_EVT,'sample_names','ascend');

%% correct roi_type_EVT to unite villi
ind_VCT_FV = find(ismember(roi_type_EVT.SegmentDisplayName,{'VCT','FV'}));
roi_type_EVT.SegmentDisplayName(ind_VCT_FV) = {'Villi'};

%% Calc median by EVT type
[G,id]          = findgroups(roi_type_EVT.SegmentDisplayName);
med_by_type     = splitapply(@median,counts_norm_tab_DATA_cl{:,2:end}',G);
temp            = array2table(med_by_type');
med_by_type_tab = [counts_norm_tab_DATA_cl(:,1) temp ];
med_by_type_tab.Properties.VariableNames = ['TargetName' cellstr(id')];

%% Make Deseq2 input files
roi_type_cur = roi_type_EVT;
counts_cl_cur = DATA_cl;

cat1s = {'Villi','Villi','Villi','EVT_A','EVT_A','EVT_I'}';
cat2s = {'EVT_A','EVT_I','EVT_E','EVT_I','EVT_E','EVT_E'}';

cat_col_name            = 'SegmentDisplayName';
Tr_keep                 = 4;
set_below_loq_const     = 0;
time_col_name           = 'GA';

comp_tab = table(cat1s,cat2s);
for c = 1:size(comp_tab,1)
    c
    cat1_name             = comp_tab.cat1s{c};
    cat2_name             = comp_tab.cat2s{c};
    
    indin_samples         = find((roi_type_cur.(cat_col_name)==cat1_name | roi_type_cur.(cat_col_name)==cat2_name ) & roi_type_cur.(time_col_name)>0);

    counts                = counts_cl_cur{:,2:end};
    bin_detected          = counts(:,indin_samples)>repmat(roi_type_cur.LOQ(indin_samples)',size(counts(:,indin_samples),1),1);

    N_segs_detected_in    = sum(bin_detected,2);

    indin_genes           = find(N_segs_detected_in>=Tr_keep);
    norm_tab_cur          = counts_cl_cur(indin_genes,[1 indin_samples'+1]);
    roi_type_test         = roi_type_cur(indin_samples,:);
    
        
    if set_below_loq_const==1
        bin_detected_cl = bin_detected(indin_genes, :);
        global_loq      = geomean(roi_type_cur.LOQ_norm);
        set_global_loq  = find(bin_detected_cl==0);
        norm_tab_cur_counts = norm_tab_cur{:,2:end};
        norm_tab_cur_counts(set_global_loq) = global_loq;
        norm_tab_cur{:,2:end} = norm_tab_cur_counts;
    end
    
    name = [cat1_name '_' cat2_name ]; %'_' time_col_name
    
    writetable(roi_type_test,['/Users/innaa/Documents/PEGBM_ns/DESeq2 inputs/meta_data_' name '.csv'],...
    'Delimiter',',')

    writetable(norm_tab_cur,['/Users/innaa/Documents/PEGBM_ns/DESeq2 inputs/counts_' name '.csv'],...
    'Delimiter',',')
end

%% Run Deseq2  in R - then continue 

%% set TRs
p_tr   = 0.05;
fc_tr  = 0.5;
q_tr   = 0.25;

%% Load DESeq2 results and select DEG
res_path = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Deseq2 output files/';
list     = dir([res_path '*.csv']); 
DEG_all  = {};
for i = 1:length(list)
    file_name = list(i).name;
    file_name
    deseq_table = readtable([res_path file_name]);
    DEG_cur = deseq_table.Var1(find(deseq_table.pvalue<=p_tr & deseq_table.padj<=q_tr & abs(deseq_table.log2FoldChange)>=fc_tr));
    DEG_all = [DEG_all;DEG_cur];
end

DEG_all = unique(DEG_all);

%% Get DEG expression data
indin      = get_gene_inds_of_list1_in_list2(DEG_all,med_by_type_tab.TargetName);
exp_dat    = med_by_type_tab(indin,:);

inter_traj = [exp_dat.Villi, exp_dat.EVT_A, exp_dat.EVT_I, exp_dat.EVT_E];
gene_name  = exp_dat.TargetName;

%% Drop lowly expressed DEG
exp_tr  = 50; 
max_exp = max(inter_traj,[],2);
indin   = find(max_exp>=exp_tr);

DEG_all = DEG_all(indin);
inter_traj = inter_traj(indin,:);
gene_name  = gene_name(indin);

%% calculate COM for trending genes
max_vec   = max(abs(inter_traj),[],2);
trend_mat = inter_traj./repmat(max_vec,1,size(inter_traj,2));
com       = sum(trend_mat.*repmat(1:size(trend_mat,2),size(trend_mat,1),1),2)./sum(trend_mat,2);

%% Make all gene traj tab
all_genes_traj      = med_by_type_tab(:,[5,2,4,3]);
all_genes_traj_norm = array2table( normalize(all_genes_traj{:,:},2));
all_traj_norm = [med_by_type_tab.TargetName all_genes_traj_norm];
all_traj_norm.Properties.VariableNames{1} = 'gene_name';

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

%% Do K-means
N_k_man  = 4;
dist     = 'sqeuclidean'; %
plot_mat = mean_by_stage_table_norm{:,2:5};

rng default; % For reproducibility
[idx,C,sumd] = kmeans(plot_mat,N_k_man,'Distance',dist,'Replicates',100,'OnlinePhase','on');

% re-order clusters by stage of peak
peak_stage_mean = nan(N_k_man,1);
for i=1:N_k_man
    inds = find(idx==i);
    peak_stage_mean(i) = peak_loc(plot_mat(inds,:));
end
cur_cluster = (1:N_k_man)';
switch_table = table(cur_cluster, peak_stage_mean);
switch_table = sortrows(switch_table,'peak_stage_mean','ascend');
idx_new      = zeros(size(idx));
for i=1:N_k_man
    inds = find(idx==switch_table.cur_cluster(i));
    idx_new(inds) = i;
end

% Finish trajectory table
mean_by_stage_table_norm.k_cluster = idx_new;
 
%% plot clustering heatmap
save_flag = 0;
mean_by_stage_table_norm    = sortrows(mean_by_stage_table_norm,{'k_cluster','com'},{'ascend','ascend'});
f_names                     = mean_by_stage_table_norm.gene_name;
for i=1:numel(f_names)
    f_names{i} = regexprep(f_names{i},'_',' ');
end
cluster_ch = find(mean_by_stage_table_norm.k_cluster(2:end)- mean_by_stage_table_norm.k_cluster(1:end-1));
lines_ch   = cluster_ch+0.5;
figure('Renderer', 'painters', 'Position', [100 100 700 1500])
imagesc(mean_by_stage_table_norm{:,2:5});
hold on
ylm = ylim;
for i=1:numel(lines_ch)
    plot(ylm,[lines_ch(i) lines_ch(i)],'k','linewidth',2);
end
% 
% yticks(1:numel(idx));
% yticklabels(f_names);

xticks(1:5);
xticklabels({'Villi','EVT-A','EVT-I','EVT-E'});
xtickangle(45)

colorbar;
set(gca,'LineWidth',2)
box on
set(gca,'FontSize',20)

if save_flag==1
    path  = '/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/Glycan EVT paper/Inna tables code/Analysis code with input file/';
    name = [path 'EVT_DEG_heatmap_DESeq2_updated'];
    h=gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(h,name,'-dpdf','-r0')
%     close all
    writetable(mean_by_stage_table_norm,[name '.csv'],'Delimiter',',','WriteRowNames',true)
end

