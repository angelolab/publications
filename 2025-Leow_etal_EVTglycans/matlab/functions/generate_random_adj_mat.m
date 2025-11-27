function [G ] = generate_random_adj_mat(DE_gly_nodes,k_nodes ,enz_nodes,gly_type_nodes,node_names,adj_mat_orig)

%% create adj mat
adj_mat     = adj_mat_orig;%zeros(numel(node_names),numel(node_names)); 
ind_k_nodes = get_gene_inds_of_list1_in_list2(k_nodes,node_names);


% do gly->k edges 
for i=1:numel(DE_gly_nodes)
    cur_name         = DE_gly_nodes{i};
    cur_k            = randi([1,numel(k_nodes)],1,1);
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_k_node_ind   = find(strcmpi(num2str(cur_k),node_names));
    % zero out prev k
    adj_mat(cur_name_node_ind,ind_k_nodes) = zeros(1,numel(ind_k_nodes));
    adj_mat(ind_k_nodes,cur_name_node_ind) = zeros(numel(ind_k_nodes),1);
    % set new one
    adj_mat(cur_name_node_ind,cur_k_node_ind)        = 1;
    adj_mat(cur_k_node_ind,cur_name_node_ind)        = 1;
end

% do enz->k edges 
for i=1:numel(enz_nodes)
    cur_name          = enz_nodes{i};
    cur_k             = randi([1,numel(k_nodes)],1,1);
    cur_name_node_ind = find(strcmpi(cur_name,node_names));
    cur_k_node_ind    = find(strcmpi(num2str(cur_k),node_names));
    % zero out prev k
    adj_mat(cur_name_node_ind,ind_k_nodes) = zeros(1,numel(ind_k_nodes));
    adj_mat(ind_k_nodes,cur_name_node_ind) = zeros(numel(ind_k_nodes),1);
    % set new one
    adj_mat(cur_name_node_ind,cur_k_node_ind)        = 1;
    adj_mat(cur_k_node_ind,cur_name_node_ind)        = 1;
end

% % do gly->types edges
% for i=1:numel(DE_gly_nodes)
%     cur_name          = DE_gly_nodes{i};
%     cur_name_node_ind = find(strcmpi(cur_name,node_names));
%     N_types           = sum(adj_mat_orig(cur_name_node_ind,:))-1;
%     rand_type_inds    = randi([1,numel(gly_type_nodes)],1,N_types );
%     cur_types         = gly_type_nodes(rand_type_inds );
%     for j=1:numel(cur_types)
%         cur_type_ind    = find(strcmpi(cur_types{j},node_names));
%         adj_mat(cur_name_node_ind,cur_type_ind)        = 1;
%         adj_mat(cur_type_ind,cur_name_node_ind)        = 1;
%     end
% end
% 
% % do enz->types edges
% for i=1:numel(enz_nodes)
%     cur_name          = enz_nodes{i};
%     cur_name_node_ind = find(strcmpi(cur_name,node_names));
%     N_types           = sum(adj_mat_orig(cur_name_node_ind,:))-1;
%     rand_type_inds    = randi([1,numel(gly_type_nodes)],1,N_types );
%     cur_types         = gly_type_nodes(rand_type_inds );
%     for j=1:numel(cur_types)
%         cur_type_ind    = find(strcmpi(cur_types{j},node_names));
%         adj_mat(cur_name_node_ind,cur_type_ind)        = 1;
%         adj_mat(cur_type_ind,cur_name_node_ind)        = 1;
%     end
% end

G = graph(adj_mat);
G.Nodes.Name = node_names;

end