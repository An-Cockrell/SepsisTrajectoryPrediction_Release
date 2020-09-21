# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 05/15/2020
# Script that will generate MLP Cytokine prediction models with a permanent dropout layer

import numpy as np
from keras import backend as K
from keras.layers.core import Lambda
import os
import math
import sys
import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.layers import Flatten
# import tensorflow.python.util.deprecation as deprecation
# deprecation._PRINT_DEPRECATION_WARNINGS = False #suppress depreciation warnings


print("\n\ninputPredictionParallel\n\n")


taskID = int(sys.argv[1])
print(taskID)
numTasks = int(sys.argv[2])

def PermaDropout(rate):
    return Lambda(lambda x: K.dropout(x, level=rate))


# 5 observations for the network to predict on
numForPredictions = 5

# percent of training data to use
percent = 100

# creating a directory to store the models in
pathToModelsDir = "./inputPredictionParallelMLPModels/"
pathToTrainingDataDir = "./inputPredictionParallelMLPTrainingMetrics/"

try:
    os.mkdir(pathToModelsDir)
except:
    pass
try:
    os.mkdir(pathToTrainingDataDir)
except:
    pass

# parameters for making the MLP neural network
def makeMLP():
    model = Sequential()
    model.add(Dense(units=1000,input_shape=(numForPredictions,11)))
    model.add(PermaDropout(0.01))
    model.add(Dense(units = 500))
    model.add(Flatten())
    model.add(Dense(units = 500))
    model.add(Dense(units = 500))
    # model.add(Dense(units = 112))
    model.add(Dense(units = 1))
    # opt = keras.optimizers.Adam(epsilon=0.01)
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_absolute_error'])
    model.summary()

    return model

# -------------------------------------#
# main section

data = np.load("trainingData.npy")
# ignore oxygen deficit
data = data[1:,:,:]
data = np.transpose(data,(2,1,0))
# expect shape = (X,6,11)
np.random.shuffle(data)

# make the dataset smaller
percent = percent/100
endPoint = math.floor(data.shape[0]*percent)
data = data[0:endPoint,:,:]

print(data.shape)

for i in range(1):


    MLPModel = makeMLP()
    print("created MLP" + str(taskID))


    checkpoint_filepath = str(pathToTrainingDataDir + "model" + str(taskID) + "_"+ str(i) + "-{epoch:02d}.h5")
    model_checkpoint_callback = keras.callbacks.ModelCheckpoint(
        filepath=checkpoint_filepath,
        save_weights_only=False)

    features = data[:,1:,:] # shape = (X, 5, 11)
    labels = data[:,0:1,taskID] # shape = (X, 1)

    hist = MLPModel.fit(features,labels, validation_split = .1, epochs = 7, batch_size = 10000, verbose = 2, callbacks=[model_checkpoint_callback])

    # np.save(pathToTrainingDataDir + "valLoss_" + str(taskID) + "_" + str(i), hist.history['val_loss'])
    # np.save(pathToTrainingDataDir + "valMAE_" + str(taskID) + "_" + str(i), hist.history['val_mean_absolute_error'])
    # np.save(pathToTrainingDataDir + "loss_" + str(taskID) + "_" + str(i), hist.history['loss'])
    # np.save(pathToTrainingDataDir + "MAE_" + str(taskID) + "_" + str(i), hist.history['mean_absolute_error'])

    MLPModel.save(pathToModelsDir + "model" + str(taskID) + ".h5")
