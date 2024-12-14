import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import pandas as pd
import numpy as np
import quiche as qu
import anndata
import quiche as qu
import sys
import ast
import time
def load_anndata_with_retry(file_path, retries=5, wait_time=2):
    """
    Tries to load an AnnData object with retries.
    
    Args:
        file_path (str): Path to the AnnData file.
        retries (int): Number of retry attempts.
        wait_time (int): Time to wait between retries (in seconds).
    
    Returns:
        AnnData object or None if it fails.
    """
    for attempt in range(retries):
        try:
            adata = anndata.read_h5ad(file_path)  # Attempt to load the file
            print(f"Successfully loaded AnnData object on attempt {attempt + 1}.")
            return adata
        except OSError as e:  # Catch file access issues
            print(f"Failed to load AnnData file: {e}. Retrying ({attempt + 1}/{retries})...")
            time.sleep(wait_time)
    print("Failed to load AnnData object after multiple retries.")
    return None
##read in arguments for analysis
spatial_method = eval(sys.argv[1])
spatial_method_params = ast.literal_eval(str(sys.argv[2]))

eval_method = eval(sys.argv[3])
eval_method_params = ast.literal_eval(str(sys.argv[4]))

print(spatial_method)
print(eval_method)

base_dir = os.path.join('data', 'simulated', 'structured')
save_directory = os.path.join(base_dir, 'metrics')

try:
    qu.pp.make_directory(save_directory)
except:
    pass

method_id = spatial_method.__name__
metric = eval_method.__name__

if method_id == 'run_quiche':
    if spatial_method_params['khop'] != None:
        param_id = 'khop'
    elif spatial_method_params['delaunay'] == True:
        param_id = 'delaunay'
    elif spatial_method_params['n_neighbors'] != None:
        param_id = 'knn'
    else:
        param_id = 'radius'
    if spatial_method_params['merge'] == True:
        param_merge = 'merged'
    else:
        param_merge = 'original'
    param_id = param_id + '_' + spatial_method_params['label_scheme'] + '_' + spatial_method_params['sig_key'] + '_' + param_merge
elif (method_id == 'evaluate_cell_charter') | (method_id == 'evaluate_kmeans'):
    if spatial_method_params['n_clusters'] is not None:
        param_id = str(spatial_method_params['n_clusters'])
    else:
        param_id = 'auto'
else:
    param_id = 'default'

pct_change_list = [0.05, 0.08, 0.1, 0.2, 0.25]
radius_list = [100, 250, 500]
prevalence_list = [0.2, 0.4, 0.6, 0.8, 1.0]
n_regions = 1
evaluation_df = pd.DataFrame()
for trial in range(0, 5):
    for pct_change in pct_change_list:
        for radius in radius_list:
            for prevalence in prevalence_list:
                try:
                    if pct_change in [0.08, 0.2] and radius != 500:
                        continue
                    print(trial, pct_change, radius, prevalence)
                    prev = int(prevalence*100)
                    pct = int(pct_change*100)
                    fig_id = f'prev{prev}_r{radius}_p{pct}_trial{trial}'
                    adata = load_anndata_with_retry(os.path.join(base_dir, 'adata', f'adata_{fig_id}.h5ad'))
                    benchmarker = qu.tl.benchmark(adata = adata, spatial_method = spatial_method, spatial_method_params = spatial_method_params,
                                            eval_method = eval_method, eval_method_params = eval_method_params)
                    sig_niches, eval_df = benchmarker.benchmark()
                    eval_df['pct_change'] = pct_change
                    eval_df['radius'] = radius
                    eval_df['prevalence'] = prevalence
                    eval_df['trial'] = trial
                    evaluation_df = pd.concat([evaluation_df, eval_df], axis = 0)
                    del adata
                except Exception as e:
                    continue
evaluation_df.to_csv(os.path.join(save_directory, f'{method_id}_{metric}_{param_id}.csv'))
