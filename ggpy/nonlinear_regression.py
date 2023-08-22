
import os
import sys
import pandas as pd
import numpy as np

from sklearn.model_selection import (RandomizedSearchCV,
                                     StratifiedKFold,
                                     StratifiedShuffleSplit)
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance

class FeatureImportance:

    def __init__(self,
                 features_file,
                 all_ggi_results,
                 n_jobs = 4,
                 n_estimators = 500,
                 boots = 1,
                 test_prop = 0.35,
                 n_repeats = 30, # iterationsf for FI
                #  cv_random_search = 2,
                 cv_size = 5, # hyperparameter space
                 suffix = 'demo',
                 out_folder = "../demo",
                 gini_based = False # FI from gini impurity metric
                ) -> None:
        
        # self.features_file = features_file
        # self.all_ggi_results = all_ggi_results
        
        self.new_df = pd.read_csv(features_file, sep='\t')
        self.ggi_pd = pd.read_csv(all_ggi_results, sep = '\t')

        self.n_jobs = n_jobs
        self.n_estimators = n_estimators
        self.boots = boots
        self.test_prop = test_prop
        self.n_splits = round(1/self.test_prop) # est


        self.n_repeats = n_repeats
        # self.cv_random_search = cv_random_search
        self.cv_size = cv_size


        self.param_grid = {
            'max_depth'       : np.linspace(   5, 100, self.cv_size, dtype=int),
            'min_samples_split': np.linspace(  2,  50, self.cv_size, dtype=int),
            'max_leaf_nodes': np.linspace(   100, 250, self.cv_size, dtype=int) 
        }

        # self.suffix = suffix
        # self.out_folder = out_folder

        self.out_file = os.path.join(out_folder, 'rf_FI_%s.csv' % suffix) # est
        self.params_file = os.path.join(out_folder, 'rf_params_%s.txt' % suffix)

        self.gini_based = gini_based

    def get_params(self):
        return self.param_grid
    
    def set_params(self, max_depth, min_samples_split, max_leaf_nodes):
        self.param_grid = {
            'max_depth'       : np.linspace(   max_depth[0], max_depth[1], self.cv_size, dtype=int),
            'min_samples_split': np.linspace(  min_samples_split[0],  min_samples_split[1], self.cv_size, dtype=int),
            'max_leaf_nodes': np.linspace(   max_leaf_nodes[0], max_leaf_nodes[1], self.cv_size, dtype=int) 
        }

    def do_resampling_dis(self, X_train_new, y_train):
        """
        balacing labels in function of the p-values
        """
        # train_num2.shape
        # tree id 1
        train_labels = np.argmax(y_train, axis=1) == 0

        pos_features = X_train_new[ train_labels]
        neg_features = X_train_new[~train_labels]
        

        pos_labels = y_train[ train_labels]
        neg_labels = y_train[~train_labels]


        if len(pos_features) < len(neg_features):

            ids = np.arange(len(pos_features))

            # taking as much as neg features are
            # available
            choices = np.random.choice(ids, len(neg_features)) 

            pos_features  = pos_features[choices]
            pos_labels    = pos_labels[choices]

        if len(pos_features) > len(neg_features):

            ids = np.arange(len(neg_features))

            # taking as much as pos features are
            # available
            choices = np.random.choice(ids, len(pos_features)) 

            neg_features = neg_features[choices]
            neg_labels   = neg_labels[choices]

        # res_pos_features.shape
        resampled_features = np.concatenate([pos_features, neg_features], axis=0)
        resampled_labels   = np.concatenate([pos_labels, neg_labels], axis=0)

        order = np.arange(len(resampled_labels))
        
        np.random.shuffle(order)

        resampled_features = resampled_features[order]
        resampled_labels   = resampled_labels[order]

        return (resampled_features, 
                resampled_labels  ,)

    def _scaler(self, ref, dat, include_clip = True):

        standarizer = StandardScaler().fit(ref)
        sted = standarizer.transform(dat)

        if include_clip:
            return np.clip(sted, -5, 5)
        
        else:
            return sted

    def transform_data(self, ref, dat):

        return self._scaler(ref, dat, True)

    def make_au_labels(self, ggi_pd, set_pd):
        # ggi_pd,set_pd = ggi_pd, new_df
        ggi_pd['tree_id'] = ggi_pd['tree_id'].astype(int  )
        ggi_pd['au_test'] = ggi_pd['au_test'].astype(float)

        out_labels = []
        has_ggi    = []
        
        for seq in set_pd['aln_base']:
            # seq
            tmp_df = ggi_pd[ ggi_pd['alignment'] == seq ]

            if len(tmp_df) < 2:
                has_ggi += [False]
                continue

            has_ggi += [True]

            out_labels.append(

                tmp_df[['tree_id', 'au_test']]
                    .sort_values('tree_id')
                    ['au_test']
                    .tolist()
            )

        labs_2d = np.array(out_labels)
        return labs_2d,set_pd.iloc[has_ggi,:]

    def data_prep(self, new_df_num, all_labels):

        all_labels_dis = np.argmax( all_labels, axis=1 ) == 0

        split = StratifiedShuffleSplit(n_splits = 1,
                                        test_size = self.test_prop, 
                                        random_state = 42)
        
        for train_index, _ in split.split(new_df_num, all_labels_dis):
            # train_index, test_index
            X_train = new_df_num.iloc[train_index,:]
            # X_test  = new_df_num.iloc[test_index,:]

            y_train = all_labels[train_index]
            # y_test  = all_labels[test_index]

        X_train_new = self.transform_data(X_train, X_train)
        # X_test_new  = transform_data(X_train, X_test)

        # do_resampling_dis(X_train_new, y_train)
        return self.do_resampling_dis(X_train_new, y_train)

    def hyperparameter_tunning(self, new_df_num, all_labels):
        # X_train_new, y_train = data_prep(new_df_num, all_labels)

        (resampled_features, 
         resampled_labels   ) = self.data_prep(new_df_num, all_labels)

        base_estimator = RandomForestRegressor( n_estimators=self.n_estimators,
                                                n_jobs=self.n_jobs, 
                                                verbose=1 )
        
        rsearch = RandomizedSearchCV(base_estimator, 
                                     self.param_grid, 
                                     cv = self.n_splits).fit(resampled_features,
                                                                     resampled_labels)

        best_rf_params = rsearch.best_params_
        # hyperparameter tunning

        # params_file = os.path.join(out_folder, 'rf_params_%s.txt' % suffix)
        with open(self.params_file, 'w') as f:
            f.write(  str(best_rf_params) + "\n" )

        return best_rf_params

    def regression_beta(self):
        # n_splits = round(1/self.test_prop) # est

        all_labels, new_df = self.make_au_labels( self.ggi_pd, self.new_df )
        new_df_num = new_df.drop(["aln_base"], axis = 1)

        # hyperparameter tunning
        best_rf_params = self.hyperparameter_tunning(new_df_num, all_labels)

        total_iter = self.boots*self.n_splits
        # n_epochs = 15
        cvscores = np.zeros((total_iter,new_df_num.shape[1]))
        sys.stdout.write(f'\n\n')

        cv  = 0
        # rank ones
        all_labels_dis = np.argmax( all_labels, axis=1 ) == 0

        for b in range(self.boots):
            # b = 1
            kfold = StratifiedKFold(
                n_splits = self.n_splits, 
                shuffle = True,
                random_state = None
            )
            for train, test in kfold.split( new_df, all_labels_dis ):
                # train,test
                # len(test)/len(train)
                sys.stdout.write(f'boot: {cv + 1}/{total_iter}\n')

                X_train = new_df_num.iloc[train,:]
                y_train = all_labels[train]

                X_test = new_df_num.iloc[test,:]
                y_test = all_labels[test]

                X_train_new = self.transform_data(X_train, X_train)
                X_test_new  = self.transform_data(X_train, X_test)
                resampled_features, resampled_labels = self.do_resampling_dis(X_train_new, y_train)

                base_estimator = self.get_base_estimator(
                                        X = resampled_features,
                                        y = resampled_labels,
                                        hyperparams = best_rf_params,
                                    )
                # permutation-based
                cvscores[cv,:] = self.FI(base_estimator, X_test_new, y_test)
                cv += 1

        self.write_FI(cvscores)
    
    def write_FI(self, cvscores):

        feature_names = self.new_df.columns.to_list()[1::]
    
        stacked = np.dstack((
                        cvscores.mean(axis=0), 
                        cvscores.std(axis=0)**2)
                )[0]
        fi = pd.DataFrame(data = stacked, columns=['mean', 'var'], index=feature_names)
    
        fi['features'] = feature_names

        (fi[['features', 'mean', 'var']]
        .sort_values('mean', ascending=False)
        .to_csv(self.out_file, sep = '\t', index=False)
        )

    def FI(self, base_estimator, X_test_new = None, y_test = None):
        if not self.gini_based:
            # permutation-based
            r = permutation_importance(base_estimator, 
                                            X_test_new, 
                                            y_test,
                                            n_repeats=self.n_repeats,
                                            random_state=0, 
                                            n_jobs=self.n_jobs)
            return r['importances_mean']
        
        else:
            # gini-based
            return base_estimator.feature_importances_

    def get_base_estimator(self, X, y, hyperparams):

        base_estimator = RandomForestRegressor(
                            n_estimators=self.n_estimators, 
                            n_jobs=self.n_jobs, 
                            verbose=0, 
                            **hyperparams 
                        )
        base_estimator.fit(X, y)
        return base_estimator


# # data --------------------------------------
# suffix = 'GAAA'
# out_folder = "../demo"

# # Features: qcutils results
# features_file = '../demo/features_stats_demo.tsv'
# # Labels: GGI results
# all_ggi_results = '../demo/out_ggi_demo.txt'
# # data --------------------------------------

# self = FeatureImportance(
#     features_file = features_file,
#     all_ggi_results = all_ggi_results,
#     n_jobs = 4,
#     n_estimators = 10,
#     boots = 1,
#     test_prop = 0.35,
#     n_repeats = 30, # iter for Permutations
#     cv_random_search = 2, 
#     cv_size = 5,
#     suffix = suffix,
#     out_folder = out_folder,
#     gini_based = True # FI from gini impurity metric
# )
# # self.get_params()
# self.regression_beta()

# self.new_df


