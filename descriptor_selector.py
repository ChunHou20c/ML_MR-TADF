from sklearn.datasets import load_diabetes
from sklearn.ensemble import RandomForestRegressor
from boruta import BorutaPy
import pandas as pd
import numpy as np

res = pd.read_csv('./data/generated_data/updated_a_new-8-Jan.csv')
features = pd.read_csv('./data/generated_data/descriptors.csv')


# merge the features and result first
merged_df = features.merge(res, left_on='Name', right_on='Filename', how='outer')

# cleaining the dataset before training

# then create the X and Y set
y = merged_df[['DeltaEST']]
print(y)

y = y.apply(lambda x: x.replace(np.inf, x[x != np.inf].max())
                         .replace(-np.inf, x[x != -np.inf].min()))
y = y.squeeze()
feature_list = features.columns.values.tolist()
feature_list.remove('Name') # to remove the name from the list

X = merged_df[feature_list].astype(np.float32)

X = X.apply(lambda x: x.replace(np.inf, x[x != np.inf].max())
                         .replace(-np.inf, x[x != -np.inf].min()))
# # let's initialize a RF model 
model = RandomForestRegressor(n_estimators=100, max_depth=5, random_state=42, n_jobs=-1)

# # let's initialize Boruta
feat_selector = BorutaPy(
    verbose=2,
    estimator=model,
    n_estimators='auto',
    max_iter=50,  # number of iterations to perform
    perc = 90
)

# # train Boruta
# # N.B.: X and y must be numpy arrays
feat_selector.fit(np.array(X), np.array(y))

# print support and ranking for each feature
print("\n------Support and Ranking for each feature------")
print(feat_selector.support_)

filter_list = zip(feat_selector.support_, X.columns)
keep_features = filter(lambda x: x[0], filter_list)
keep_feature_list = [i[1] for i in keep_features]

selected_x = merged_df[keep_feature_list].to_csv('./target/training_data/input.csv', index=False)
selected_y = merged_df[['DeltaEST']].to_csv('./target/training_data/output.csv', index=False)

# for i in range(len(feat_selector.support_)):
#     if feat_selector.support_[i]:
#         print("Passes the test: ", X.columns[i],
#               " - Ranking: ", feat_selector.ranking_[i])
#     else:
#         continue
        # print("Doesn't pass the test: ",
        #       X.columns[i], " - Ranking: ", feat_selector.ranking_[i])
