function [mean_by_stage_table_norm] = k_heatmap(N_k_man,mean_by_stage_table_norm,types,plot_flag,plot_row_names)
dist     = 'sqeuclidean'; %
plot_mat = mean_by_stage_table_norm{:,2:numel(types)+1};

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
try
    mean_by_stage_table_norm           = sortrows(mean_by_stage_table_norm,{'k_cluster','com'},{'ascend','ascend'});
catch
    mean_by_stage_table_norm           = sortrows(mean_by_stage_table_norm,{'k_cluster'},{'ascend'});
end

% plot
if plot_flag==1
    
    cluster_ch = find(mean_by_stage_table_norm.k_cluster(2:end)- mean_by_stage_table_norm.k_cluster(1:end-1));
    lines_ch   = cluster_ch+0.5;
    figure('Renderer', 'painters', 'Position', [100 100 700 1500])
    % imagesc(plot_mat);
    imagesc(mean_by_stage_table_norm{:,2:numel(types)+1});
    hold on
    ylm = ylim;
    for i=1:numel(lines_ch)
        plot(ylm,[lines_ch(i) lines_ch(i)],'k','linewidth',2);
    end
    %
    yticks(1:numel(idx));
    
    % plot row names
    if plot_row_names==1
        f_names                     = cellstr(mean_by_stage_table_norm{:,1});
        
        for i=1:numel(f_names)
            f_names{i} = regexprep(f_names{i},'_',' ');
        end
        yticklabels(f_names);
        %
    end
    
    xticks(1:numel(types));
    xticklabels(types);
    xtickangle(45)
    
    colorbar;
    set(gca,'LineWidth',2)
    box on
    set(gca,'FontSize',20)
end

end

