%% path
addpath('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/scripts')  

%% Read tables
lasso_tab = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/Lasso classifier tables/classifier_coeffs_summary.csv');
test_tab = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/feature_table.csv','PreserveVariableNames',true);
DCIS_Metadata = readtable('/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/DCIS_Metadata.csv','VariableNamingRule' , 'preserve');

%% Calculate feature skewness
ind_first_pred = 6;
out_method     = "median";
thresh         = 1.75;
tab_all        = test_tab;
Predictor      = tab_all.Properties.VariableNames(ind_first_pred:end)';
skew_tab = table(Predictor, nan(numel(Predictor),1), nan(numel(Predictor),1));
skew_tab.Properties.VariableNames{2} = 'skewness_p';
skew_tab.Properties.VariableNames{3} = 'skewness_np';
for f=1:size(skew_tab,1)
    cur_ft         = skew_tab.Predictor{f};
    cur_ft_vals    = table(tab_all.(cur_ft),categorical(tab_all.event_recur_type));
    cur_ft_vals.Properties.VariableNames{1} = 'values';
    cur_ft_vals.Properties.VariableNames{2} = 'labels';
    cur_ft_vals(isnan(cur_ft_vals.values),:) = [];
    overall_mean          = mean(cur_ft_vals.values);
    overall_median        = median(cur_ft_vals.values);
    cur_ft_vals.out       = isoutlier(cur_ft_vals.values,out_method,'ThresholdFactor',thresh );
    cur_ft_p_vals         = cur_ft_vals(find(cur_ft_vals.labels=='Invasive_Ipsilateral'),:);
    cur_ft_np_vals        = cur_ft_vals(find(cur_ft_vals.labels=='Non-progressor'),:);
    skew_tab.skewness_p(f)  = sum(cur_ft_p_vals.out & cur_ft_p_vals.values>overall_median)/numel(cur_ft_p_vals.out);
    skew_tab.skewness_np(f) = sum(cur_ft_np_vals.out & cur_ft_np_vals.values>overall_median)/numel(cur_ft_np_vals.out);
end

sk_tr = 0.05;
p_all = signrank(skew_tab.skewness_p,skew_tab.skewness_np,'tail','left')
skew_tab.skew_diff = skew_tab.skewness_np-skew_tab.skewness_p;
frac_all = sum(sign(skew_tab.skew_diff).*skew_tab.skewness_np>sk_tr)/numel(skew_tab.skew_diff)

skew_tab_lasso_ft = innerjoin(lasso_tab(:,{'Predictor','mean_coeff'}),skew_tab,'keys','Predictor');
p_lasso = signrank(skew_tab_lasso_ft.skewness_p,skew_tab_lasso_ft.skewness_np,'tail','left')
frac_lasso = sum(sign(skew_tab_lasso_ft.skew_diff).*skew_tab_lasso_ft.skewness_np>sk_tr)/numel(skew_tab_lasso_ft.skew_diff)

skew_tab_lasso_ft = sortrows(skew_tab_lasso_ft,'mean_coeff','descend');

%% Classify skewness per feature
skew_tab_lasso_ft.is_pos_skewed_p  = zeros(size(skew_tab_lasso_ft,1),1);
skew_tab_lasso_ft.is_pos_skewed_np = zeros(size(skew_tab_lasso_ft,1),1);

skew_tab_lasso_ft.is_pos_skewed_p(find(sign(skew_tab_lasso_ft.skewness_p).*skew_tab_lasso_ft.skewness_p>sk_tr)) = 1;
skew_tab_lasso_ft.is_pos_skewed_np(find(sign(skew_tab_lasso_ft.skewness_np).*skew_tab_lasso_ft.skewness_np>sk_tr)) = 1;

skew_tab_lasso_ft.skew_cat(find(skew_tab_lasso_ft.is_pos_skewed_np==1 & skew_tab_lasso_ft.is_pos_skewed_p==1)) = {'both'};
skew_tab_lasso_ft.skew_cat(find(skew_tab_lasso_ft.is_pos_skewed_np==0 & skew_tab_lasso_ft.is_pos_skewed_p==0)) = {'none'};
skew_tab_lasso_ft.skew_cat(find(skew_tab_lasso_ft.is_pos_skewed_np==0 & skew_tab_lasso_ft.is_pos_skewed_p==1)) = {'Progressors only'};
skew_tab_lasso_ft.skew_cat(find(skew_tab_lasso_ft.is_pos_skewed_np==1 & skew_tab_lasso_ft.is_pos_skewed_p==0)) = {'Non-progressors only'};

%% Select skewed features 
skewed_fts_lasso = skew_tab_lasso_ft.Predictor(find(strcmp(skew_tab_lasso_ft.skew_cat,'Non-progressors only')));

%% save skew table
writetable(skew_tab_lasso_ft,'/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/tables/high_outliers_in_lasso_fts.csv');

%% Definitions for LR fts 
p_large            = 0.05; 
keep_lasso_fts     = 1;
rem_cell_cell_dist = 0;

% Clean up feature table 
test_tab_l   = test_tab;
fts_lasso = skewed_fts_lasso;
if rem_cell_cell_dist==1
    featuresWithDot = fts_lasso(contains(fts_lasso, '.') & ~contains(fts_lasso, '0.9'));
    fts_lasso_sel = setdiff(fts_lasso,featuresWithDot);
else
    fts_lasso_sel = fts_lasso;
end

ind_lasso_fts = get_gene_inds_of_list1_in_list2(fts_lasso_sel, test_tab_l.Properties.VariableNames);
ind_lasso_fts = ind_lasso_fts(ind_lasso_fts>0);
test_tab_l = test_tab_l(:,[1:(ind_first_pred-1) ind_lasso_fts' ]);
    
%% Remove non differential cols
end_col_ind     = size(test_tab_l,2);
label_col       = 'Status';
rem = [];
for i=ind_first_pred:end_col_ind
    values = (test_tab_l{:,i});
    labels = test_tab_l.(label_col);
    try
    indfinite = find(isfinite(values));
    catch
%         cur_mc
        continue
    end
    values = values(indfinite);
    labels = labels(indfinite);
    p_val = kruskalwallis(values,labels,'off');
    if p_val>=p_large
        rem = [rem i];
    end
end
test_tab_l(:,rem) = [];

% rank and test features
ind_first_ft = ind_first_pred;

all_test_tab_l_fts = size(test_tab_l,2)-ind_first_pred+1;

%% Create normalized feature table
T_norm_man        = test_tab_l;
for i=ind_first_pred:size(T_norm_man,2)
    T_norm_man{:,i}   = normalize(T_norm_man{:,i});
    indnan = find(isnan(T_norm_man{:,i}));
    T_norm_man{indnan,i} = 0;
end

%% Quantify risk score performance by combination

predictors = T_norm_man.Properties.VariableNames(ind_first_pred:end)';
% Prepare the results table
results_table = table();

% Counter for storing performance results
performance_counter = 1;

% Step 4: Measure performance for subsets of predictors
num_p = numel(predictors);
for num_to_keep = [1:num_p]  % Include full sets

    % Get all combinations of markers being discarded
    combinations = nchoosek(1:num_p, num_to_keep);
    
    % Iterate through each combination
    for i = 1:size(combinations, 1)
        
        % Create a new list of predictors by excluding the discarded ones
        remaining_predictors = combinations(i,:);
        
        % Proceed only if there are remaining predictors
        if ~isempty(remaining_predictors)
            predictor_subset =predictors(remaining_predictors);
            
            % Step 5: Measure performance with the current subset of predictors
            current_performance = get_performance(predictor_subset,T_norm_man);
            results_table(performance_counter, :) = {numel(predictor_subset), current_performance, predictor_subset};
            performance_counter = performance_counter + 1;
        end
    end
end

% Set appropriate variable names for results table
results_table.Properties.VariableNames = {'Num_predictors', 'Performance','Predictors_selected'};

%% Plot predictor num dependence
res_tab_un = unique(results_table(:,[1 2]));
G = groupsummary(res_tab_un,["Num_predictors"],"max","Performance");
figure
% bar(G.max_Performance)
bar(G.max_Performance,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
ylabel('Predictive performance')
xlabel('Number of predictors')
set(gca,'FontSize',18 )
set(gca,'LineWidth',2)
set(gcf,'color','w');
box on

%% Select optimal combination for risk score
results_table = sortrows(results_table,'Performance','descend');
low_r_fts = results_table.Predictors_selected{1,1};
cols_keep = [T_norm_man.Properties.VariableNames(1:(ind_first_pred-1)) low_r_fts(1:end)'];
T_norm_sel = T_norm_man(:,{cols_keep{1:end}});

%% Calc risk score per patient and plot
T_norm_sel.score = sum(T_norm_sel{:,ind_first_pred:end},2);
T_norm_sel.is_low_risk = T_norm_sel.score>max(T_norm_sel.score(find(strcmp(T_norm_sel.event_recur_type,'Invasive_Ipsilateral'))));

figure;
violinplot(T_norm_sel.score,T_norm_sel.event_recur_type);
p_val       = kruskalwallis(T_norm_sel.score,T_norm_sel.event_recur_type,'off');
title(['Risk Score, p = ' num2str(p_val)])
xtickangle(45)
% xticklabels({'Invasive Progressor','Non Progressor'})
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(gca, 'TickLabelInterpreter', 'none')
% xlabel('Progression status')
ylabel('scaled score')
set(gca, 'FontSize', 18);
set(gcf,'color','w');

%% Add risk score to meta data table
try
    DCIS_Metadata = removevars(DCIS_Metadata, {'is_low_risk','score'});
catch
end
DCIS_Metadata     = outerjoin(DCIS_Metadata,T_norm_sel(:,{'fov','is_low_risk','score'}),'LeftKeys','MIBI_ID','RightKeys','fov','type','left');
DCIS_Metadata     = removevars(DCIS_Metadata, 'fov');


