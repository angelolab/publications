import pandas as pd
import numpy as np

# read data
#df_feature = pd.read_csv('data/genomics_df.csv')
df_feature = pd.read_csv('data/genomics/processed_genomics_features_final.csv')
print(df_feature.head())
print(np.unique(df_feature['Timepoint']))

#all_features = set(df_feature['feature_name'])
# RNA features only
all_features = set(df_feature[df_feature['data_type'] == 'DNA']['feature_name'])
print(df_feature.shape)
print(len(all_features))

for each_time_point in ['baseline']:#np.unique(df_feature['Timepoint']):
    print(each_time_point)
    # loop over df_feature
    pts_dict = {}
    feature_list = set()
    for i in range(df_feature.shape[0]):
        row = df_feature.iloc[i]
        pts_id = row['Patient_ID']
        timepoint = row['Timepoint']
        type = row['data_type']
        if type == 'DNA':
            if timepoint == each_time_point: #'baseline':
                if pts_id not in pts_dict:
                    pts_dict[pts_id] = {}
                feature_name = row['feature_name']
                if feature_name not in pts_dict[pts_id]:
                    pts_dict[pts_id][feature_name] = row['normalized_value' ]

    # loop over all fov and features
    for pts_id in pts_dict:
        for feature_name in all_features:
            if feature_name not in pts_dict[pts_id]:
                pts_dict[pts_id][feature_name] = 0

    df_matrix = pd.DataFrame.from_dict(pts_dict, orient='index')
    print(df_matrix.head())
    print(df_matrix.shape)
    file_name = './data/genomics/dna_df_matrix_' + each_time_point + '.csv'
    print(file_name)
    df_matrix.to_csv(file_name)
