import os
import numpy as np
import quiche as qu
import os
import gc

### 20 patients, 5000 cells, grid size, balanced cell type distribution
save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
save_directory_ = os.path.join(save_directory_, 'balanced')
da_vec_A = ['A', 'C', 'E']
da_vec_B = ['B', 'D']
n_regions = 1
n_patients_condA = 10
n_patients_condB = 10
ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
grid_size_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 14]
sample_size_A = {'A': 1000, 'B': 1000, 'C': 1000, 'D':1000, 'E': 1000}
sample_size_B = {'A': 1000, 'B': 1000, 'C': 1000, 'D':1000, 'E': 1000}

n_niches_A = np.array(list(sample_size_A.values())).sum()
n_niches_B = np.array(list(sample_size_B.values())).sum()

random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
for trial in range(0, 5):
    save_directory = os.path.join(save_directory_, f'trial{trial}')
    for grid_size in grid_size_list:
        for ratio in ratio_list:
            A_id_join = ''.join(da_vec_A)
            B_id_join = ''.join(da_vec_B)
            ratio_id = str(ratio).replace('.', '_')
            fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
            qu.pp.make_directory(save_directory)
            adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
                                                          n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
                                                            random_state_list_A = random_state_list_A, scale = 2048,
                                                            random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
                                                            fig_id = fig_id, save_directory=save_directory)
            adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
            del adata_simulated
            gc.collect()

# ### 20 patients, 5000 cells, grid size, rare cell type distribution, even, ACE
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'rare', 'even', 'ACE')
# da_vec_A = ['A', 'C', 'E']
# da_vec_B = ['B', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size_list = [3, 4, 5, 6, 7, 8, 9, 10, 14]
# sample_size_A = {'A': 250, 'B': 1000, 'C': 1000, 'D':1000, 'E': 1000}
# sample_size_B = {'A': 250, 'B': 1000, 'C': 1000, 'D':1000, 'E': 1000}
# n_niches_A = np.array(list(sample_size_A.values())).sum()
# n_niches_B = np.array(list(sample_size_B.values())).sum()

# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for grid_size in grid_size_list:
#         for ratio in ratio_list:
#             print(grid_size, ratio)
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)    
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### 20 patients, 5000 cells, grid size, rare cell type distribution, even, BD
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'rare', 'even', 'BD')
# da_vec_A = ['A', 'C', 'E']
# da_vec_B = ['B', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size_list = [4, 5, 6, 7, 8, 9, 10, 14]
# sample_size_A = {'A': 1000, 'B': 250, 'C': 1000, 'D':1000, 'E': 1000}
# sample_size_B = {'A': 1000, 'B': 250, 'C': 1000, 'D':1000, 'E': 1000}
# n_niches_A = np.array(list(sample_size_A.values())).sum()
# n_niches_B = np.array(list(sample_size_B.values())).sum()
# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for grid_size in grid_size_list:
#         for ratio in ratio_list:
#             print(grid_size, ratio)
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B, 
#                                                             fig_id = fig_id, save_directory=save_directory)     
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### 20 patients, 5000 cells, grid size, rare cell type distribution, differences in total abundance, AB
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'rare', 'uneven', 'AB')
# da_vec_A = ['A', 'B']
# da_vec_B = ['C', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size_list = [4, 5, 6, 7, 8, 9, 10, 14]
# sample_size_A = {'A': 250, 'B': 1000, 'C': 250, 'D':1000, 'E': 1000}
# sample_size_B = {'A': 250, 'B': 1000, 'C': 250, 'D':1000, 'E': 2000}
# n_niches_A = np.array(list(sample_size_A.values())).sum()
# n_niches_B = np.array(list(sample_size_B.values())).sum()

# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for grid_size in grid_size_list:
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)        
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### 20 patients, 5000 cells, grid size, rare cell type distribution with over enriched niche, differences in total abundance, AE
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'rare', 'uneven', 'AE')
# da_vec_A = ['A', 'E']
# da_vec_B = ['C', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size_list = [4, 5, 6, 7, 8, 9, 10, 14]
# sample_size_A = {'A': 250, 'B': 1000, 'C': 250, 'D':1000, 'E': 2000}
# sample_size_B = {'A': 250, 'B': 1000, 'C': 250, 'D':1000, 'E': 1000}
# n_niches_A = np.array(list(sample_size_A.values())).sum()
# n_niches_B = np.array(list(sample_size_B.values())).sum()

# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for grid_size in grid_size_list:
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)   
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### 20 patients, cell size (1000, 2500, 5000, 7500, 10000), even distribution ACE, BD, 5 trials
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'cell_size')
# da_vec_A = ['A', 'C', 'E']
# da_vec_B = ['B', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size = 5 #moderately sized niche across all conditions
# n_niches_list = [1000, 2500, 5000, 7500, 10000]
# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for n_niches in n_niches_list:
#         n_niches_A = n_niches
#         n_niches_B = n_niches
#         sample_size_A = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}
#         sample_size_B = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_niches{n_niches}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)   
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### sample size, even, 5000 cells, 20 patients (10/10), 50 patients (25/25), 100 patients (50/50)
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'sample_size', 'even')
# da_vec_A = ['A', 'C', 'E']
# da_vec_B = ['B', 'D']
# n_regions = 1
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size = 5 #moderately sized niche across all conditions
# n_niches_list = [1000, 2500, 5000, 7500, 10000]
# random_state_list_A_total = [[58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828],
#                         [18255, 375048, 841289, 358212, 865561, 910340, 932356, 991428, 452710, 939069, 824313, 860195, 257761,  46185, 247190, 202698, 330544, 121697, 648207, 283553, 169522, 486357, 858414, 304794, 205472],
#                         [408012, 458676, 717973, 668781, 403886, 606634, 260306, 939276, 623881, 255186, 367705, 528345, 226844, 728294, 639584, 685932, 910892, 486827, 637358, 814176, 698718, 713803, 166737, 881948, 468107, 821048, 309228, 222704, 161452, 923765, 952348, 928353, 378920, 163929, 984270, 492804, 755216, 436558, 318103, 607161, 615113, 399706, 433870, 499257, 598094, 387245, 674469, 536081, 414155,  96618]]
# random_state_list_B_total = [[51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420],
#                         [4350, 986258, 852314, 102436, 879390, 348535, 455353, 761904, 476135, 961161, 355725,  87752, 232359, 617969, 216322, 901868, 751788, 166704, 285442, 399510, 715552, 755666, 478920, 955944, 548483],
#                         [997559,  54808, 293390,  94168, 788780,  70251, 490678, 443910, 526926, 873664, 654329, 744592, 264408, 514690, 896512, 440709, 473844, 485323, 932347, 732826, 439922, 289889, 141581, 394896, 971390, 266777, 498694, 934209, 262397, 597341, 975235, 476137, 322417, 757786, 146101,  10980, 750593, 705199, 134491, 238890, 260925, 457256, 247559, 156366, 690471, 974356, 531052, 691124, 755871, 587338]]
# n_niches = 5000
# n_niches_A = n_niches
# n_niches_B = n_niches
# sample_size_A = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}
# sample_size_B = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}

# sample_size_list = [20, 50, 100]

# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for j in range(0, len(sample_size_list)):
#         print(trial, j)
#         sample_size = sample_size_list[j]
#         n_patients_condA = int(sample_size/2)
#         n_patients_condB = int(sample_size/2)
#         random_state_list_A = random_state_list_A_total[j]
#         random_state_list_B = random_state_list_B_total[j]
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_s{sample_size}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### sample size, uneven, 5000 cells, 100 patients (10/90, 20/80, 30/70, 40/60, 50/50)
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'sample_size', 'uneven')
# da_vec_A = ['A', 'C', 'E']
# da_vec_B = ['B', 'D']
# n_regions = 1
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size = 5 #moderately sized niche across all conditions
# n_niches = 5000
# n_niches_A = n_niches
# n_niches_B = n_niches
# sample_size_A = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}
# sample_size_B = {'A': int(n_niches/5), 'B': int(n_niches/5), 'C': int(n_niches/5), 'D':int(n_niches/5), 'E': int(n_niches/5)}

# n_patients_condA_list = [10, 20, 30, 40, 50]
# n_patients_condB_list = [90, 80, 70, 60, 50]

# random_state_list_A_total = [[240493, 395437, 949113, 242241, 505889, 185155, 213423, 799357, 177846, 141381],
#                              [413861, 349394, 516693, 678970,  16658, 814168, 965118, 219289, 99359, 647259, 799261, 350303, 816069, 309915, 986835,  79676, 647856, 988874, 768169, 463018],
#                              [698038, 351199, 550538, 700690, 177952, 674634, 147269, 794813, 509790, 298059, 936270, 325970,  85210, 545660, 620427, 874482, 221946, 388787,  90398, 431601, 255886, 170482, 791677, 754641,299848,   1132, 445568, 120415, 876597, 205897],
#                              [383711, 192552, 884589, 444274, 500596, 680477, 543995, 449702, 98161, 467592, 204710, 968154, 221870, 669431, 371414, 562171, 163868, 791125, 721641, 689501, 420910, 930694, 133836, 651289, 637379,  66283, 702001, 607649, 746632, 819674, 473316, 442550, 722067, 272730, 773553, 511921, 361239,  33873,  93175, 213860],
#                              [848220, 927783, 436358, 393034, 290584, 346851, 584048, 292324, 772694,  85461, 588152, 828031, 625201, 569494, 884678, 922410, 32609,  87289, 183729,  31852, 220947, 825632, 383676, 487401, 812266, 721738, 323297, 666727, 684737, 311361, 719237, 517523, 780623, 448674, 356669, 220739, 379194, 493796, 111230, 867165, 102437, 402966, 815466, 117721,  14149, 538800, 601240, 265469, 688721, 591421]]


# random_state_list_B_total = [[61802, 414200,  20197, 293304, 512599, 741984, 229224, 278317, 250480, 241110,  25618, 401229, 155034,  81538, 105018, 760766, 232141, 651535, 287083, 274530, 812626,  76827, 339341, 297657, 160544, 859692, 193674, 639873,  57680,  90786, 831413, 915333, 71259, 563580, 625740, 191691, 323496, 732286, 239377, 483896, 347613, 588464, 657482, 121438, 802330, 206688, 662808, 761130, 530189, 784919, 726053,  85726, 745948, 672210, 382170, 289213, 526205, 417845, 120708, 853704, 142103, 377233, 482941, 266932,613082, 248783, 653193, 989563, 493920, 242516, 898643, 559490, 840194, 564683,  49732, 495212, 991987, 977246, 283154, 916101, 925336, 532583, 909261, 640759, 471242, 988223, 675821, 505717,73817, 856696],
#                              [485387, 992361, 511881, 189214, 518096, 213581, 653281, 286891, 521372, 319780, 129910, 299823, 653856, 640596,  28695, 945730, 780724, 413971, 240258, 731740, 498076, 940057, 411469, 480085, 18911, 909985, 272588, 792964, 330715, 123420, 798392, 236039, 995011, 280878, 305705, 784310, 864351, 668517,  75106, 627050, 316897, 186907, 509909, 211946, 111777,  19943, 401213, 241122, 827007, 231140,  39956, 699512, 343228, 146887, 265478,   8537, 131627, 304460, 448877, 679082, 154824, 847025, 714171, 278313, 117797, 714513, 694812, 851913, 602395, 953011, 176222, 238522, 699949,  30855, 377731, 953191, 520921, 498301, 961603, 460237],
#                              [20255, 783098, 510646, 258881, 618949, 719710, 686936,  52240, 253871, 593534, 723463, 104974, 137611, 707107, 576473, 240522, 79029, 277864, 792963, 671086, 424303, 347881, 129531, 990201, 461269, 769505, 254001, 587175, 881861, 375850, 963339, 745837, 953318,  77041, 879977, 875850, 400844, 992128, 429900, 803638, 922389, 321646, 358693, 757608, 485085, 185538, 129185, 441015, 758491,   3511, 948575, 606790, 310362, 628820, 789512, 943630, 680024, 377328, 791435, 493682, 555570,  59164,  65847,  19524, 297697,  75108, 167191, 223313, 482020, 432454],
#                              [296452, 843693, 447635, 900452, 612330,  19363, 723769, 186457, 990400, 484646, 574270, 431561, 438256, 921178, 713983, 224570, 647193, 734795, 915604, 664428, 944849, 960489, 844272, 859931, 256300, 459001, 878128, 429649, 259785,  45524, 158554, 331419, 851949,  31504,  70714, 681962, 480359,  93007, 323777, 512935, 420134, 925511, 324753, 905563, 900765, 454024, 205436, 525031, 182752, 479156, 999501, 850433,  65928, 913099, 557852, 214287,19474, 693876, 470806, 935838],
#                              [7073, 5411, 9234, 9677, 4828, 6462, 5897, 1378, 7219, 8395, 2577, 7873, 9876, 6683, 9043,  353,  235, 1264, 4222, 5172, 4488, 7001, 4307, 4132, 2599, 4077, 1216, 9100,  915, 7614, 4203, 3192, 8493, 615, 3225, 7712, 1171, 6333, 2335, 1804, 8258, 4239, 2979, 6446, 8977, 3711, 9795, 4635, 2403, 2307]]

# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for j in range(0, len(n_patients_condA_list)):
#         n_patients_condA = n_patients_condA_list[j]
#         n_patients_condB = n_patients_condB_list[j]
#         random_state_list_A = random_state_list_A_total[j]
#         random_state_list_B = random_state_list_B_total[j]
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_s{n_patients_condA}_s{n_patients_condB}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                           n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                             random_state_list_A = random_state_list_A, scale = 2048,
#                                                             random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                             fig_id = fig_id, save_directory=save_directory)  
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()

# ### 20 patients, 5000 cells, grid size, rare cell type distribution, differences in total abundance, AB
# save_directory_ = r'data/simulated/unstructured/adata/n5000/t20'
# save_directory_ = os.path.join(save_directory_, 'balanced', 'uneven', 'AB')
# da_vec_A = ['A', 'B']
# da_vec_B = ['C', 'D']
# n_regions = 1
# n_patients_condA = 10
# n_patients_condB = 10
# ratio_list = [0.2, 0.4, 0.6, 0.8, 1.0]
# grid_size_list = [4, 5, 6, 7, 8, 9, 10, 14]
# sample_size_A = {'A': 1000, 'B': 1000, 'C': 1000, 'D':1000, 'E': 1000}
# sample_size_B = {'A': 1000, 'B': 1000, 'C': 1000, 'D':1000, 'E': 2000}

# n_niches_A = np.array(list(sample_size_A.values())).sum()
# n_niches_B = np.array(list(sample_size_B.values())).sum()

# random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828]
# random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420]
# for trial in range(0, 5):
#     save_directory = os.path.join(save_directory_, f'trial{trial}')
#     for grid_size in grid_size_list:
#         for ratio in ratio_list:
#             A_id_join = ''.join(da_vec_A)
#             B_id_join = ''.join(da_vec_B)
#             ratio_id = str(ratio).replace('.', '_')
#             fig_id = A_id_join+'_'+B_id_join+f'_grid{grid_size}_ratio{ratio_id}_trial{trial}'
#             qu.pp.make_directory(save_directory)
#             adata_simulated = qu.tl.simulate_unstructured(n_patients_condA = n_patients_condA, n_patients_condB = n_patients_condB, num_grids_x = grid_size, num_grids_y = grid_size, ratio = ratio, n_niches_A = n_niches_A, n_niches_B = n_niches_B,
#                                                         n_regionsA = n_regions, n_regionsB = n_regions, da_vec_A = da_vec_A, da_vec_B = da_vec_B,
#                                                         random_state_list_A = random_state_list_A, scale = 2048,
#                                                         random_state_list_B = random_state_list_B, sample_size_A = sample_size_A, sample_size_B = sample_size_B,
#                                                         fig_id = fig_id, save_directory=save_directory)        
#             cells = adata_simulated[adata_simulated.obs['Patient_ID'] == 'B0'][adata_simulated[adata_simulated.obs['Patient_ID'] == 'B0'].obs['cell_cluster'] == 'E'].obs_names[:1000]
#             adata_simulated = adata_simulated[~np.isin(adata_simulated.obs_names, cells)]    
#             adata_simulated.write_h5ad(os.path.join(save_directory_, f'adata_simulated_unstructured_{fig_id}.h5ad'))
#             gc.collect()