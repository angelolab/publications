function [cyc_tab,cyc_per_type_tab] = cyc_summary(G,DE_gly_nodes,k_nodes ,enz_nodes,gly_type_nodes)

%% find cycles and annotate them by type and k clust
[cycles,edgecycles] = allcycles(G,'MaxCycleLength',4);
bin_good_cyc        = zeros(numel(cycles),1);
k_cyc               = zeros(numel(cycles),1);
type_cyc            = {};
enz_cyc             = {};
gly_cyc             = {};
% filter cycles
for i=1:numel(cycles)
    cur_cyc = cycles{i};
    N_enz   = sum(ismember(cur_cyc,enz_nodes));
    N_gly   = sum(ismember(cur_cyc,DE_gly_nodes ));
    N_k     = sum(ismember(cur_cyc,k_nodes ));
    if N_enz==1 && N_gly==1 && N_k==1
        bin_good_cyc(i)=1;
    else
        type_cyc{i} = 'nan';
        enz_cyc{i}  = 'nan';
        gly_cyc{i}  = 'nan';
        continue
    end
    k_str   = cur_cyc{ismember(cur_cyc,k_nodes)};
    k_cyc(i)    = str2double(k_str);
    type_cyc{i} = cur_cyc{ismember(cur_cyc,gly_type_nodes)};
    enz_cyc{i}  = cur_cyc{ismember(cur_cyc,enz_nodes)};
    gly_cyc{i}  = cur_cyc{ismember(cur_cyc,DE_gly_nodes)};
end

cycles     = cycles(bin_good_cyc==1);
edgecycles = edgecycles(bin_good_cyc==1);
k_cyc      = k_cyc(bin_good_cyc==1);
type_cyc   = type_cyc(bin_good_cyc==1);
enz_cyc    = enz_cyc(bin_good_cyc==1);
gly_cyc    = gly_cyc(bin_good_cyc==1);

type_cyc   = type_cyc';
enz_cyc    = enz_cyc';
gly_cyc    = gly_cyc';

cyc_tab = table(cycles,edgecycles,enz_cyc,gly_cyc,k_cyc,type_cyc);
cyc_tab = sortrows(cyc_tab,{'k_cyc','type_cyc'},{'ascend','ascend'});

cyc_per_type = zeros(numel(gly_type_nodes),1);
for i=1:numel(gly_type_nodes)
    cyc_per_type(i) = numel(find(strcmpi(gly_type_nodes{i},cyc_tab.type_cyc)));
end

cyc_per_type_tab = table(gly_type_nodes,cyc_per_type);

end