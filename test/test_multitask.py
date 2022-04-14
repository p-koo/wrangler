import os, h5py
import numpy as np
from tensorflow import keras
import tensorflow as tf



def custom_model(input_shape, num_labels, params=[128, 196, 256, 512], activation='relu'):
    l2_reg = tf.keras.regularizers.L2(1e-6)

    inputs = keras.layers.Input(shape=input_shape)

    # layer 1
    nn = keras.layers.Conv1D(filters=128, kernel_size=19, use_bias=False, padding='same',
                             kernel_regularizer=l2_reg)(inputs)     
    nn = keras.layers.BatchNormalization()(nn)
    nn = keras.layers.Activation(activation)(nn)
    nn = keras.layers.MaxPool1D(pool_size=20)(nn)
    nn = keras.layers.Dropout(0.1)(nn)

    # layer 2
    nn = keras.layers.Conv1D(filters=196, kernel_size=7, use_bias=False, padding='same', activation=None,
                             kernel_regularizer=l2_reg)(nn)      
    nn = keras.layers.BatchNormalization()(nn)
    nn = keras.layers.Activation('relu')(nn)
    nn = keras.layers.MaxPool1D(pool_size=5)(nn)
    nn = keras.layers.Dropout(0.2)(nn)

    # layer 3
    nn = keras.layers.Conv1D(filters=256, kernel_size=5, use_bias=False, padding='valid', activation=None,
                             kernel_regularizer=l2_reg)(nn)      
    nn = keras.layers.BatchNormalization()(nn)
    nn = keras.layers.Activation('relu')(nn)
    nn = keras.layers.MaxPool1D(pool_size=2)(nn)
    nn = keras.layers.Dropout(0.3)(nn)

    nn = keras.layers.Flatten()(nn)

    # layer 4 - Fully-connected 
    nn = keras.layers.Dense(512, activation=None, use_bias=False, 
                            kernel_regularizer=l2_reg)(nn)  
    nn = keras.layers.BatchNormalization()(nn)
    nn = keras.layers.Activation('relu')(nn) 
    nn = keras.layers.Dropout(0.5)(nn)

    # Output layer
    logits = keras.layers.Dense(num_labels, activation='linear', use_bias=True)(nn)
    outputs = keras.layers.Activation('sigmoid')(logits)

    # create keras model
    return keras.Model(inputs=inputs, outputs=outputs)


cell_line = 'HCT116'
filepath = '/home/koolab/peter/data/tf/binary/'+cell_line+'.h5'
dataset = h5py.File(filepath, 'r')

with h5py.File(filepath, 'r') as dataset:
    x_train = np.array(dataset['x_train']).astype(np.float32)
    y_train = np.array(dataset['y_train']).astype(np.float32)
    x_valid = np.array(dataset['x_valid']).astype(np.float32)
    y_valid = np.array(dataset['y_valid']).astype(np.int32)
    x_test  = np.array(dataset['x_test']).astype(np.float32)
    y_test  = np.array(dataset['y_test']).astype(np.int32)

N, L, A = x_train.shape
num_labels = y_valid.shape[1]
print(x_train.shape)


model = custom_model(input_shape=(L,A), num_labels=num_labels,
                     params=[128, 196, 256, 512], activation='exponential')

# set up optimizer and metrics
auroc = keras.metrics.AUC(curve='ROC', name='auroc')
aupr = keras.metrics.AUC(curve='PR', name='aupr')
optimizer = keras.optimizers.Adam(learning_rate=0.0005)
loss = keras.losses.BinaryCrossentropy(from_logits=False, label_smoothing=0.0)
model.compile(optimizer=optimizer,
                loss=loss,
                metrics=[auroc, aupr])


# early stopping callback
es_callback = keras.callbacks.EarlyStopping(monitor='val_aupr', #'val_aupr',#
                                            patience=10, 
                                            verbose=1, 
                                            mode='max', 
                                            restore_best_weights=True)
# reduce learning rate callback
reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='val_aupr', 
                                                factor=0.2,
                                                patience=3, 
                                                min_lr=1e-7,
                                                mode='max',
                                                verbose=1) 

# train model
history = model.fit(x_train, y_train, 
                    epochs=100,
                    batch_size=64, 
                    shuffle=True,
                    validation_data=(x_valid, y_valid), 
                    callbacks=[es_callback, reduce_lr])


model.evaluate(x_test, y_test)

