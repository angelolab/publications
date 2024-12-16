import sys
import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator
import quiche as qu
import time

class benchmark(BaseEstimator):
    def __init__(
        self,
        adata = None,
        spatial_method = None,
        spatial_method_params = None,
        eval_method = None,
        eval_method_params = None,
        **kwargs
    ):
        """Class for benchmarking spatial enrichment methods

        Parameters
        adata: anndata.AnnData (default = None)
            annotated data object containing preprocessed single-cell data (dimensions = cells x features)
        spatial method: function (default = None)
            function housed in quiche.tl.evaluate script that specifies which type of spatial enrichment to perform. Can be one of the following, or you can provide your own.
                qu.tl.run_quiche
                qu.tl.evaluate_pairwise
                qu.tl.evaluate_graphcompass
                qu.tl.evaluate_kmeans
                qu.tl.evaluate_cell_charter
        spatial_method_params: dictionary (default = None)
            dictionary referring to the spatial enrichment method hyperparameters. For more information on method-specific hyperparameters, see the quiche.tl.evaluate script for the method of interest.
                quiche example: {'radius': 200,'labels_key':'cell_cluster', 'spatial_key':'spatial', 'fov_key':'Patient_ID', 'patient_key':'Patient_ID','khop':3, 'n_neighbors': 10, 'delaunay': False, 'min_cells':5, 'k_sim':100, 'design':'~condition', 'model_contrasts':'conditionA-conditionB', 'sketch_size':None, 'nlargest': 5, 'annotation_key':'quiche_niche', 'n_jobs':-1, 'label_scheme':'fov_norm', 'sig_key':'SpatialFDR', 'merge':False}
        eval_method: function (default = None)
            function housed in the quiche.tl.metrics script that specifies the evaluation method to perform. Can be one of the following (or you may provide your own):
                niche recall: qu.tl.group_recall
                niche purity: qu.tl.evaluate_purity
        eval_method_params: dictionary (default = None)
            dictionary referring to the evaluation method hyperparameters
                group_recall example for quiche: {'method_type':'spatial_cluster', 'feature_key':'spatial_nhood', 'ground_key': 'DA_group', 'labels_key': 'quiche_niche', 'ground_truth_niches':['A_C_E', 'B_D']}
                group_purity example for quiche: {'annot_key':'quiche_niche', 'labels_key':'cell_cluster', 'fov_key':'Patient_ID', 'condition_key':'DA_group','feature_key':'spatial_nhood'}
        ----------
        Attributes
        benchmarker.perform_enrichment()
            Performs spatial enrichment according to method of interest
            Returns:
                adata: anndata.AnnData
                    annotated data object containing enrichment results for spatial clustering methods
                sig_niches: pd.DataFrame
                    dataframe containing significant niches and their associated p-values
        benchmarker.evaluate_enrichment()
            Performs spatial enrichment according to method of interest, then performs evaluation according to metric of interest
            Returns:
                eval_df: pd.DataFrame
                    dataframe containing performance scores
        ----------
        """
        self.adata = adata
        self.spatial_method = spatial_method
        self.spatial_method_params = spatial_method_params
        self.eval_method = eval_method
        self.eval_method_params = eval_method_params
        self.kwargs = kwargs
        
        if self.kwargs is None:
            self.kwargs = {}

        if self.spatial_method_params is None:
            self.spatial_method_params = {}

        if self.eval_method_params is None:
            self.eval_method_params = {}
            
    def perform_enrichment(self):
        sys.stdout.write('performing enrichment: {}'.format(self.spatial_method.__name__)+'\n')
        if self.spatial_method.__name__ == 'evaluate_pairwise':
            self.sig_niches = self.spatial_method(adata = self.adata, **self.spatial_method_params)
            return self.sig_niches
        else:
            self.adata, self.sig_niches = self.spatial_method(adata = self.adata, **self.spatial_method_params)
            return self.adata, self.sig_niches

    def evaluate_enrichment(self):
        sys.stdout.write('performing evaluation: {}'.format(self.eval_method.__name__)+'\n')

        self.eval_df = self.eval_method(mdata = self.adata, scores_df = self.sig_niches, **self.eval_method_params)
        return self.eval_df

    def benchmark(self):
        if self.spatial_method.__name__ == 'evaluate_pairwise':
            self.sig_niches = self.perform_enrichment()
        else:
            self.adata, self.sig_niches = self.perform_enrichment()
        self.eval_df = self.evaluate_enrichment()
        self._aggregate() 
        return self.sig_niches, self.eval_df

    def compute_runtime(self):
        tic = time.perf_counter()
        self.perform_enrichment()
        toc = time.perf_counter()
        time_df = pd.DataFrame([toc-tic], columns = ['runtime'])
        time_df = time_df.melt()
        return time_df

    def _aggregate(self):
        p = str(self.spatial_method.__name__)
        self.eval_df['method'] = p