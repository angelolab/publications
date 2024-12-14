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
run = str(sys.argv[5])
print(spatial_method)
print(eval_method)
base_dir = os.path.join('data', 'simulated', 'unstructured')
if run == 'balanced':
    param_dir = os.path.join('n5000', 't20', 'balanced')
    da_vec_A = ['A', 'C', 'E']
    da_vec_B = ['B', 'D']
elif run == 'balanced_uneven':
    param_dir = os.path.join('n5000', 't20', 'balanced', 'uneven', 'AB')
    da_vec_A = ['A', 'B']
    da_vec_B = ['C', 'D']
elif run == 'rare_even_ACE':
    param_dir = os.path.join('n5000', 't20', 'rare', 'even', 'ACE')
    da_vec_A = ['A', 'C', 'E']
    da_vec_B = ['B', 'D']
elif run == 'rare_even_BD':
    param_dir = os.path.join('n5000', 't20', 'rare', 'even', 'BD')
    da_vec_A = ['A', 'C', 'E']
    da_vec_B = ['B', 'D']
elif run == 'rare_uneven_AB':
    param_dir = os.path.join('n5000', 't20', 'rare', 'uneven', 'AB')
    da_vec_A = ['A', 'B']
    da_vec_B = ['C', 'D']
elif run == 'rare_uneven_AE':
    param_dir = os.path.join('n5000', 't20', 'rare', 'uneven', 'AE')
    da_vec_A = ['A', 'E']
    da_vec_B = ['C', 'D']
elif run == 'rare_uneven_AB_corrected':
    param_dir = os.path.join('n5000', 't20', 'rare', 'uneven', 'AB')
    da_vec_A = ['A', 'B']
    da_vec_B = ['C', 'D']
elif run == 'balanced_uneven_corrected':
    param_dir = os.path.join('n5000', 't20', 'balanced', 'uneven', 'AB')
    da_vec_A = ['A', 'B']
    da_vec_B = ['C', 'D']

save_directory = os.path.join(base_dir, 'metrics', param_dir)

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
    if 'corrected' in run:
        param_id = param_id + '_' + spatial_method_params['label_scheme'] + '_' + spatial_method_params['sig_key'] + '_' + param_merge + '_corrected'
    else:
        param_id = param_id + '_' + spatial_method_params['label_scheme'] + '_' + spatial_method_params['sig_key'] + '_' + param_merge
elif (method_id == 'evaluate_cell_charter') | (method_id == 'evaluate_kmeans'):
    if spatial_method_params['n_clusters'] is not None:
        param_id = str(spatial_method_params['n_clusters'])
    else:
        param_id = 'auto'
elif (method_id == 'run_multiscale_quiche'):
    param_id = 'multiscale_' + spatial_method_params['sig_key']
else:
    param_id = 'default'

n_regions = 1
evaluation_df = pd.DataFrame()
for trial in range(0, 5):
    ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
    grid_size_list = [3, 4, 5, 6, 7, 8, 9, 10, 14]
    for grid_size in grid_size_list:
        for ratio in ratio_list:
            try:
                A_id_join = ''.join(da_vec_A)
                B_id_join = ''.join(da_vec_B)
                ratio_id = str(ratio).replace('.', '_')
                fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
                adata = load_anndata_with_retry(os.path.join(base_dir, 'adata', param_dir, f'adata_simulated_unstructured_{fig_id}.h5ad'), retries=5, wait_time=2)
                benchmarker = qu.tl.benchmark(adata = adata, spatial_method = spatial_method, spatial_method_params = spatial_method_params,
                                          eval_method = eval_method, eval_method_params = eval_method_params)
                sig_niches, eval_df = benchmarker.benchmark()
                eval_df['grid_size'] = grid_size
                eval_df['ratio'] = ratio
                eval_df['trial'] = trial
                evaluation_df = pd.concat([evaluation_df, eval_df], axis = 0)
            except Exception as e:
                continue
evaluation_df.to_csv(os.path.join(save_directory, f'{method_id}_{metric}_{param_id}.csv'))
