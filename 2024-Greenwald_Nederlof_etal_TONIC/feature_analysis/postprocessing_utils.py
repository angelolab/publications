import os
import warnings
from typing import Callable

import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests
from python_files.utils import compare_populations

TIMEPOINT_NAMES = ['primary', 'baseline', 'pre_nivo', 'on_nivo']

# define the timepoint pairs to use, these index into harmonized_metadata
TIMEPOINT_PAIRS = [
    ("primary", "baseline"),
    ("baseline", "pre_nivo"),
    ("baseline", "on_nivo"),
    ("pre_nivo", "on_nivo")
]

# define the timepoint columns as they appear in timepoint_df
TIMEPOINT_COLUMNS = [
    "primary__baseline", "baseline__pre_nivo", "baseline__on_nivo", "pre_nivo__on_nivo"
]


def euclidean_timepoint(tp_one_data: pd.Series, tp_two_data: pd.Series) -> float:
    """Compute the Euclidean distance between two timepoint data

    Args:
        tp_one_data (pd.Series):
            The data for the first timepoint
        tp_two_data (pd.Series):
            The data for the second timepoint

    Returns:
        float:
            The Euclidean distance between the two timepoint datapoints
    """
    # combine the two series into one df
    tp_combined: pd.DataFrame = pd.concat([tp_one_data, tp_two_data], axis=1)

    # drop nans across both columns
    tp_combined = tp_combined.dropna(axis=0)

    # unit normalize each column
    tp_combined = tp_combined.apply(lambda x: (x / np.linalg.norm(x)), axis=0)

    # return Euclidean distance
    return np.linalg.norm(tp_combined.values[:, :1] - tp_combined.values[:, 1:])


def generate_patient_paired_timepoints(
        harmonized_metadata: pd.DataFrame, timepoint_df: pd.DataFrame,
        distance_metric: Callable[[pd.Series, pd.Series], float], tissue_id_col: str = "Tissue_ID",
        patient_id_col: str = "Patient_ID", timepoint_col: str = "Timepoint",
        feature_to_compare: str = "normalized_mean",
        feature_to_create: str = "euclidean_distance_normalized_mean_filtered_features",
        timepoint_pairs = TIMEPOINT_PAIRS, timepoint_names=TIMEPOINT_NAMES
) -> pd.DataFrame:
    """Compute distance metric between timepoints aggregated across all evolution table features.

    Both raw_value (direct output of the distance metric) and normalized_value (the output z-scored) are computed.

    This table gets appended to the evolution table and combined into the final timepoint features table as an
    additional feature.

    Args:
        harmonized_metadata (pd.DataFrame):
            Maps each FOV and Tissue ID to the corresponding patient and timepoint
        timepoint_df (pd.DataFrame):
            Maps the features measured for each Tissue ID
        distance_metric (Callable[[pd.Series, pd.Series], float]):
            A custom distance metric used to compute the distance between timepoint data
        tissue_id_col (str):
            The column to index into the tissue ID
        patient_id_col (str):
            The column to index into the patient ID
        timepoint_col (str):
            The column to index into the timepoint
        feature_to_compare (str):
            The feature to compare across different timepoints in timepoint_col
        feature_to_create (str):
            The name of the feature to append to the existing evolution table

    Returns:
        pd.DataFrame:
            The table defining feature_to_create across all other features, appended to the evolution table as
            an aggregated feature.
    """
    # define a DataFrame that stores the distance metric feature between all other features
    col_names = [
        "feature_name_unique", "Patient_ID", "comparison", "raw_value"
    ]
    timepoint_comparisons = pd.DataFrame(columns=col_names)

    # group the metadata by patient ID
    patient_groups = harmonized_metadata[
        [tissue_id_col, patient_id_col, timepoint_col]
    ].groupby(patient_id_col)

    # iterate through each patient and their timepoint data
    for patient_id, patient_data in patient_groups:
        # get the unique tissue samples for each timepoint
        patient_data_dedup = patient_data[
            patient_data[timepoint_col].isin(
                timepoint_names
            )
        ].drop_duplicates()

        # define which tissue ID maps to which timepoint, this will help with sorting
        tissue_id_timepoint_map = dict(
            zip(patient_data_dedup[tissue_id_col], patient_data_dedup[timepoint_col])
        )

        # get the corresponding timepoint data
        timepoint_subset = timepoint_df.loc[
                           timepoint_df[tissue_id_col].isin(patient_data_dedup[tissue_id_col].values), :
                           ]

        # in the case there aren't any corresponding tissue IDs, continue
        # NOTE: this can happen because the tissue IDs between harmonized_metadata and timepoint_df
        # don't always match up
        if len(timepoint_subset) == 0:
            warnings.warn(f"Skipping patient {patient_id}, no corresponding timepoint values")
            continue

        # group into specific columns by tissue, then rename columns to corresponding timepoint
        wide_timepoint = pd.pivot(
            timepoint_subset, index="feature_name_unique", columns=tissue_id_col,
            values=feature_to_compare
        ).rename(tissue_id_timepoint_map, axis=1)

        # if a specific timepoint pair exists, then compute the mean difference across all features
        for tp in timepoint_pairs:
            if tp[0] in wide_timepoint.columns.values and tp[1] in wide_timepoint.columns.values:
                col_difference = distance_metric(
                    wide_timepoint.loc[:, tp[0]], wide_timepoint.loc[:, tp[1]]
                )
                feature_0 = tp[0]
                feature_1 = tp[1]
                timepoint_comparisons = pd.concat(
                    [
                        timepoint_comparisons,
                        pd.DataFrame(
                            [[
                                feature_to_create,
                                patient_id,
                                f"{feature_0}__{feature_1}",
                                col_difference
                            ]],
                            columns=col_names
                        )
                    ]
                )

    # z-score the raw_values and store in normalized_value
    timepoint_comparisons["normalized_value"] = (
                                                        timepoint_comparisons["raw_value"] - timepoint_comparisons[
                                                    "raw_value"].mean()
                                                ) / timepoint_comparisons["raw_value"].std()

    return timepoint_comparisons


def combine_features(analysis_dir, harmonized_metadata, timepoint_features, timepoint_features_agg,
                     patient_metadata, metadata_cols=['Patient_ID'], timepoint_columns=TIMEPOINT_COLUMNS,
                     timepoint_names=TIMEPOINT_NAMES, timepoint_pairs=TIMEPOINT_PAIRS, file_suffix='',
                     drop_immune_agg=True):
    evolution_dfs = []
    # generate evolution df based on difference in timepoints
    for evolution_col in timepoint_columns:
        timepoint_1, timepoint_2 = evolution_col.split('__')
        evolution_df = timepoint_features_agg[timepoint_features_agg[evolution_col]].copy()
        evolution_df = evolution_df.loc[evolution_df.Timepoint.isin([timepoint_1, timepoint_2])]

        # get the paired features
        evolution_df_wide = evolution_df.pivot(index=['feature_name_unique', 'Patient_ID'], columns='Timepoint',
                                               values=['normalized_mean', 'raw_mean'])
        evolution_df_wide.columns = ['_'.join(col).strip() for col in evolution_df_wide.columns.values]

        # calculate difference between normalised and raw values across timepoints
        evolution_df_wide = evolution_df_wide.reset_index()
        evolution_df_wide = evolution_df_wide.dropna(axis=0)
        evolution_df_wide['comparison'] = evolution_col
        evolution_df_wide['normalized_value'] = evolution_df_wide['normalized_mean_' + timepoint_2] - evolution_df_wide[
            'normalized_mean_' + timepoint_1]
        evolution_df_wide['raw_value'] = evolution_df_wide['raw_mean_' + timepoint_2] - evolution_df_wide[
            'raw_mean_' + timepoint_1]
        evolution_df_wide = evolution_df_wide[
            ['feature_name_unique', 'Patient_ID', 'comparison', 'normalized_value', 'raw_value']]

        evolution_dfs.append(evolution_df_wide)

    evolution_df = pd.concat(evolution_dfs)

    # add the Euclidean distance between all normalized and raw features across all patients
    aggregate_euclidean = generate_patient_paired_timepoints(
        harmonized_metadata, timepoint_features,
        distance_metric=euclidean_timepoint, timepoint_pairs=timepoint_pairs, timepoint_names=timepoint_names
    )
    aggregate_euclidean = evolution_df[["Patient_ID", "comparison"]].drop_duplicates().merge(
        aggregate_euclidean, on=["Patient_ID", "comparison"]
    )[aggregate_euclidean.columns]
    evolution_df = pd.concat([evolution_df, aggregate_euclidean])
    evolution_df.to_csv(os.path.join(analysis_dir, f"timepoint_evolution_features{file_suffix}.csv"), index=False)

    # create combined df
    timepoint_features = pd.read_csv(os.path.join(analysis_dir, f'timepoint_features_filtered{file_suffix}.csv'))
    timepoint_features = timepoint_features.merge(
        harmonized_metadata[['Patient_ID', 'Tissue_ID', 'Timepoint'] + timepoint_columns].drop_duplicates(),
        on='Tissue_ID')
    timepoint_features = timepoint_features.merge(patient_metadata[metadata_cols].drop_duplicates(), on='Patient_ID',
                                                  how='left')

    # Hacky, remove once metadata is updated
    #timepoint_features = timepoint_features.loc[timepoint_features.Patient_ID < 200, :]
    timepoint_features = timepoint_features.loc[timepoint_features.Timepoint.isin(timepoint_names), :]
    timepoint_features = timepoint_features[
        ['Tissue_ID', 'feature_name', 'feature_name_unique', 'raw_mean', 'raw_std', 'normalized_mean', 'normalized_std',
         'Timepoint'] + metadata_cols]

    # look at evolution
    evolution_df = pd.read_csv(os.path.join(analysis_dir, f'timepoint_evolution_features{file_suffix}.csv'))
    evolution_df = evolution_df.merge(patient_metadata[metadata_cols].drop_duplicates(), on='Patient_ID', how='left')
    evolution_df = evolution_df.rename(
        columns={'raw_value': 'raw_mean', 'normalized_value': 'normalized_mean', 'comparison': 'Timepoint'})
    evolution_df = evolution_df[['feature_name_unique', 'raw_mean', 'normalized_mean', 'Timepoint'] + metadata_cols]

    # combine together into single df
    combined_df = timepoint_features.copy()
    combined_df = combined_df[['feature_name_unique', 'raw_mean', 'normalized_mean', 'Timepoint'] + metadata_cols]
    combined_df = pd.concat([combined_df, evolution_df[['feature_name_unique', 'raw_mean', 'normalized_mean', 'Timepoint'] + metadata_cols]])
    combined_df['combined_name'] = combined_df.feature_name_unique + '__' + combined_df.Timepoint

    # drop any immune_agg features
    if drop_immune_agg==True:
        combined_df = combined_df[~combined_df.feature_name_unique.str.contains('immune_agg')]

    combined_df.to_csv(os.path.join(analysis_dir, f'timepoint_combined_features{file_suffix}.csv'), index=False)


def generate_feature_rankings(analysis_dir, combined_df, feature_metadata, pop_col='Clinical_benefit',
                              populations=['No', 'Yes'], timepoint_names=TIMEPOINT_NAMES, file_suffix=''):
    method = 'ttest'
    pop1, pop2 = populations

    # placeholder for all values
    total_dfs = []

    for comparison in combined_df.Timepoint.unique():
        population_df = compare_populations(feature_df=combined_df, pop_col=pop_col,
                                            timepoints=[comparison], pop_1=pop1, pop_2=pop2, method=method)

        if np.sum(~population_df.log_pval.isna()) == 0:
            continue
        long_df = population_df[['feature_name_unique', 'log_pval', 'mean_diff', 'med_diff']]
        long_df['comparison'] = comparison
        long_df = long_df.dropna()
        long_df['pval'] = 10 ** (-long_df.log_pval)
        long_df['fdr_pval'] = multipletests(long_df.pval, method='fdr_bh')[1]
        total_dfs.append(long_df)

    # summarize hits from all comparisons
    ranked_features_df = pd.concat(total_dfs)
    ranked_features_df['log10_qval'] = -np.log10(ranked_features_df.fdr_pval)

    # get ranking of each row by FDR p-value and correlation
    ranked_features_df['pval_rank'] = ranked_features_df.fdr_pval.rank(ascending=True)
    ranked_features_df['cor_rank'] = ranked_features_df.med_diff.abs().rank(ascending=False)
    ranked_features_df['combined_rank'] = (
                                                      ranked_features_df.pval_rank.values + ranked_features_df.cor_rank.values) / 2

    # generate importance score
    max_rank = len(~ranked_features_df.med_diff.isna())
    normalized_rank = ranked_features_df.combined_rank / max_rank
    ranked_features_df['importance_score'] = 1 - normalized_rank

    ranked_features_df = ranked_features_df.sort_values('importance_score', ascending=False)

    # generate signed version of score
    ranked_features_df['signed_importance_score'] = ranked_features_df.importance_score * np.sign(
        ranked_features_df.med_diff)

    # add feature type
    ranked_features_df = ranked_features_df.merge(feature_metadata, on='feature_name_unique', how='left')

    feature_type_dict = {'functional_marker': 'phenotype', 'linear_distance': 'interactions',
                         'density': 'density', 'cell_diversity': 'diversity', 'density_ratio': 'density',
                         'mixing_score': 'interactions', 'region_diversity': 'diversity',
                         'compartment_area_ratio': 'compartment', 'density_proportion': 'density',
                         'morphology': 'phenotype', 'pixie_ecm': 'ecm', 'fiber': 'ecm', 'ecm_cluster': 'ecm',
                         'compartment_area': 'compartment', 'ecm_fraction': 'ecm'}
    ranked_features_df['feature_type_broad'] = ranked_features_df.feature_type.map(feature_type_dict)

    # get ranking of each feature
    ranked_features_df['feature_rank_global_evolution'] = ranked_features_df.importance_score.rank(ascending=False)

    # get ranking of non-evolution features
    ranked_features_no_evo = ranked_features_df.loc[
                             ranked_features_df.comparison.isin(timepoint_names), :]
    ranked_features_no_evo['feature_rank_global'] = ranked_features_no_evo.importance_score.rank(ascending=False)
    ranked_features_df = ranked_features_df.merge(
        ranked_features_no_evo.loc[:, ['feature_name_unique', 'comparison', 'feature_rank_global']],
        on=['feature_name_unique', 'comparison'], how='left')

    # get ranking for each comparison
    ranked_features_df['feature_rank_comparison'] = np.nan
    for comparison in ranked_features_df.comparison.unique():
        # get subset of features from given comparison
        ranked_features_comp = ranked_features_df.loc[ranked_features_df.comparison == comparison, :]
        ranked_features_comp['temp_comparison'] = ranked_features_comp.importance_score.rank(ascending=False)

        # merge with placeholder column
        ranked_features_df = ranked_features_df.merge(
            ranked_features_comp.loc[:, ['feature_name_unique', 'comparison', 'temp_comparison']],
            on=['feature_name_unique', 'comparison'], how='left')

        # replace with values from placeholder, then delete
        ranked_features_df['feature_rank_comparison'] = ranked_features_df['temp_comparison'].fillna(
            ranked_features_df['feature_rank_comparison'])
        ranked_features_df.drop(columns='temp_comparison', inplace=True)

    # saved formatted df
    ranked_features_df.to_csv(os.path.join(analysis_dir, f'feature_ranking{file_suffix}.csv'), index=False)


def prediction_preprocessing(df_feature, prediction_dir):
    all_features = set(df_feature['feature_name_unique'])

    for each_time_point in np.unique(df_feature['Timepoint']):
        print(each_time_point)
        # loop over df_feature
        pts_dict = {}
        feature_list = set()
        for i in range(df_feature.shape[0]):
            row = df_feature.iloc[i]
            pts_id = row['Patient_ID']
            timepoint = row['Timepoint']
            if timepoint == each_time_point: #'baseline':
                if pts_id not in pts_dict:
                    pts_dict[pts_id] = {}
                feature_name = row['feature_name_unique']
                if feature_name not in pts_dict[pts_id]:
                    pts_dict[pts_id][feature_name] = row['normalized_mean']
                    # Note that we are currently using row['normalized_mean'] instead of row['raw_mean' ]

        # loop over all fov and features
        for pts_id in pts_dict:
            for feature_name in all_features:
                if feature_name not in pts_dict[pts_id]:
                    pts_dict[pts_id][feature_name] = 0

        df_matrix = pd.DataFrame.from_dict(pts_dict, orient='index')
        file_name = os.path.join(prediction_dir, 'df_matrix_' + each_time_point + '.csv')
        df_matrix.to_csv(file_name)
