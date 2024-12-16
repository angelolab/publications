import os
from sklearn.preprocessing import StandardScaler
import numpy as np

def make_directory(directory: str = None):
    """Creates a directory at the specified path if one doesn't exist.

    Parameters
    ----------
    directory : str
        A string specifying the directory path.

    Returns
    -------
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def standardize(x =  None):
    """Standardizes data by removing the mean and scaling to unit variance.

    Parameters
    x: pd.DataFrame (default = None)
        data matrix (dimensions = cells x features)
    ----------

    Returns
    X: pd.DataFrame
        standardized data matrix (dimensions = cells x features)
    ----------
    """
    scaler = StandardScaler(with_mean = True, with_std = True)
    X = scaler.fit_transform(x)
    return X

def compute_percentile(df, p = 70):
    scores = df.values
    scores = scores[~np.isnan(scores)]
    perc = np.percentile(scores, p)
    return perc

def filter_fovs(adata, patient_key, threshold):
    """Filters samples according to the number of cells/niches specifed.

    Parameters
    adata: (default = None)
        anndata object
    patient_key: str 
        string indicating filtering key
    threshold: int
        integer referring to the minimum number of niches per sample
    ----------

    Returns
    adata:
        filtered anndata object
    ----------
    """
    n_niches = adata.obs[patient_key].value_counts(sort=False)
    adata = adata[~np.isin(adata.obs[patient_key], n_niches[n_niches < threshold].index)]    
    return adata