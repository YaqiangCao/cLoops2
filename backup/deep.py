#!/usr/bin/env python3
#--coding:utf-8 --
"""
cLoops2 deep-learning related core modules. 
"""

#sklearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import mean_absolute_error as MAE

#keras
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras import metrics
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import load_model, Sequential, Model
from tensorflow.keras.callbacks import ModelCheckpoint, ReduceLROnPlateau, EarlyStopping
from tensorflow.keras.layers import Flatten, Dense, Dropout, BatchNormalization, Activation

#warning settings
import warnings
warnings.filterwarnings("ignore")


#gloabl settings
class hyperparameters:
    dim = 2
    batch_size = 4096
    learning_rate = 1e-3
    weight_decay = 0.0005
    train_steps = 1000  # will be changed accoridng to sample numbers
    vali_steps = 100
    test_steps = 100
    epochs = 100
    reduce_lr_patience = 10
    early_stop_patience = 20
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'


PARA = hyperparameters()


def getRegrModel(input_shape,
                 checkpoint=None,
                 lr=1e-3,
                 reduce_lr_patience=5,
                 early_stop_patience=10):
    """
    Get the deep-learning regression model.
    """
    if os.path.isfile(checkpoint):
        print("loading existing model")
        model = load_model(checkpoint)
    else:
        model = Sequential()
        model.add(Dense(128, input_dim=input_shape))
        model.add(BatchNormalization())
        model.add(Activation(activation="relu"))
        model.add(Dense(64, activation="relu"))
        model.add(Dropout(0.3))
        model.add(Dense(1, activation="linear"))
        model.compile(
            loss='mse',
            optimizer=Adam(lr=lr),
            metrics=['mse', 'mae'],
        )
    reduce_lr = ReduceLROnPlateau(monitor="val_loss",
                                  patience=reduce_lr_patience)
    early_stop = EarlyStopping(monitor='val_loss',
                               patience=early_stop_patience)
    callbacks = [reduce_lr, early_stop]
    if checkpoint is not None:
        cp = ModelCheckpoint(checkpoint,
                             monitor='val_loss',
                             verbose=1,
                             save_weights_only=False,
                             save_best_only=True,
                             mode='min')
    callbacks.append(cp)
    return model, callbacks
