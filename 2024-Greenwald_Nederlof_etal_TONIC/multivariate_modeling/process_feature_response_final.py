import pandas as pd
import numpy as np

# read data
df_feature = pd.read_csv('./data/mibi/combined_df.csv')
print(df_feature.head())
#df_meta = pd.read_csv('data/harmonized_metadata.csv')
print(np.unique(df_feature['Timepoint']))

all_features = set(df_feature['feature_name_unique'])
#print(all_features)
print(df_feature.shape)
print(len(all_features))
# print(df_meta.shape)

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
    print(df_matrix.head())
    print(df_matrix.shape)
    file_name = './data/mibi/processed_data/df_matrix_' + each_time_point + '.csv'
    print(file_name)
    df_matrix.to_csv(file_name)
    #df_matrix.to_csv('processed_data/df_matrix_baseline_1030.csv')
