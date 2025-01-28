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