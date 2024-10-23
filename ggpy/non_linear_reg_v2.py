import numpy as np
from phylokrr.inspection import  pdp_1d, pdp_2d, permutation_importance
from phylokrr.utils import split_data

import tensorflow as tf
from tensorflow import keras
import keras_tuner



def order_labels(all_ggi_results):
    """
    Order the labels in the same order as the features
    It skips the first line of the all_ggi_results file

    The ordering is based on the ids, which is the second column of the all_ggi_results file
    What is being ordered is the last column of the all_ggi_results file

    :param all_ggi_results: path to all_ggi_results file
    :return: ordered_results
    """

    # read tsv file all_ggi_results
    ggi_results = {}
    all_ids = []

    with open(all_ggi_results) as f:
        lines = f.readlines()
        for n,line in enumerate(lines):

            if n == 0:
                continue

            line = line.strip().split('\t')
            seq = line[0]
            id = int(line[1])
            test = float(line[-1])

            if seq not in ggi_results:
                ggi_results[seq] = {id: test}

            else:
                ggi_results[seq][id] = test

            if id not in all_ids:
                all_ids.append(id)

    ordered_ids = sorted(all_ids)

    ordered_results = {}
    for k in ggi_results.keys():
        # print(k,v)
        ordered_results[k] = [ggi_results[k][id] for id in ordered_ids]

    return ordered_results

def get_Xy(features_file, all_ggi_results):
    """
    Get the features and labels from the files

    :param features_file: path to the features file
    :param all_ggi_results: path to the all_ggi_results file
    :return: X, y, header
    """

    pre_labels = order_labels(all_ggi_results)

    y = []
    X = []
    header = []
    with open(features_file) as f:
        lines = f.readlines()
        for n,line in enumerate(lines):
            line = line.strip().split('\t')
            if n == 0:
                header += line[1:]
                continue
            # print(line)
            seq = line[0]
            y.append(pre_labels[seq])
            X.append([float(x) for x in line[1:]])

    X = np.array(X)
    y = np.array(y)

    return X, y, header



def build_model(hp, p, params):
    r_layers        = params['layers']
    r_units         = params['units']
    r_drop_out      = params['drop_out']
    r_learning_rate = params['learning_rate']
    r_decay_rate    = params['decay_rate']
    

    model = keras.Sequential()
    model.add( keras.layers.InputLayer( input_shape = p ) )

    # Tune the number of layers.
    for i in range(hp.Int("num_layers", 
                          r_layers[0],
                          r_layers[1])):
        
        drop_out = hp.Float(f"drop_{i}",
                             min_value=r_drop_out[0],
                             max_value=r_drop_out[1], sampling="log")
        
        model.add( keras.layers.BatchNormalization() )
        model.add( keras.layers.Dropout(drop_out) )
        model.add(
            keras.layers.Dense(
                # Tune number of units separately.
                units=hp.Int(f"units_{i}", 
                                min_value=r_units[0],
                                max_value=r_units[1], step=3),
                kernel_initializer = 'lecun_normal',
                use_bias = False
            )
        )
        model.add( keras.layers.Activation('selu') )

    model.add( keras.layers.Dense( units = 2, activation = 'selu' ) )

    learning_rate = hp.Float("lr", 
                             min_value=r_learning_rate[0],
                             max_value=r_learning_rate[1],
                            sampling="log")
    decay_rate = hp.Float("decay", 
                          min_value=r_decay_rate[0], 
                          max_value=r_decay_rate[1], sampling="log")
    
    optimizer = keras.optimizers.legacy.SGD(
        learning_rate = learning_rate,
        momentum = 0.90,
        nesterov = True,
        decay = decay_rate,
        )

    model.compile(
        optimizer=optimizer,
        loss = 'mse',
        metrics=[ tf.keras.metrics.CosineSimilarity(axis=1) ]
    )

    return model


def tune_dnn( X_train, y_train, X_test, y_test, project_name = 'test', params = {}):

    max_trials = params['max_trials']
    n_epochs   = params['n_epochs']

    mBM = lambda hp: build_model(hp, p = X_train.shape[1], params = params)

    tuner = keras_tuner.BayesianOptimization(
        hypermodel=mBM,
        objective="val_loss",
        max_trials=max_trials,
        overwrite=True,
        project_name=project_name, # random folder name
    )

    early_stopping_cb = keras.callbacks.EarlyStopping('val_loss', 
                                patience =100,
                                restore_best_weights=True,
                                mode = 'min'
                        )

    tuner.search(
        x = X_train,
        y = y_train,
        epochs=n_epochs,
        validation_data=(
            X_test, 
            y_test, 
        ),
        callbacks=[
            early_stopping_cb,
            # onecycle
        ],
    )




# data --------------------------------------
suffix = 'GAAA'
out_folder = "../demo"

# Features: qcutils results
features_file = '../demo/features_stats_demo.tsv'
# Labels: GGI results
all_ggi_results = '../demo/out_ggi_demo.txt'
# data --------------------------------------


X, y, header = get_Xy(features_file, all_ggi_results)




