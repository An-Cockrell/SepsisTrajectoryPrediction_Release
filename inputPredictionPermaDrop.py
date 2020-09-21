# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 05/15/2020
# Script that will generate LSTM Cytokine prediction models with a permanent dropout layer

import numpy as np
from keras import backend as K
from keras.layers.core import Lambda
import os
import math

import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.layers import Flatten
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False #suppress depreciation warnings



# make regressor do rolling prediction as well
# make lstm for just oxydef



def PermaDropout(rate):
    return Lambda(lambda x: K.dropout(x, level=rate))

# 5 observations for the network to predict on
numForPredictions = 5

# percent of training data to use
percent = 10

# creating a directory to store the models in
pathToModelsDir = "./inputPredictionPermaDropout/"

try:
    os.mkdir(pathToModelsDir)
except:
    pass

# parameters for making the LSTM neural network
def makeLSTM():
    model = Sequential()
    model.add(LSTM(units=100,return_sequences=True,input_shape=(numForPredictions,11)))
    model.add(LSTM(units=100,return_sequences=False,input_shape=(numForPredictions,11)))
    model.add(PermaDropout(0.1))
    model.add(Dense(units = 300))
    model.add(Dense(units = 200))
    # model.add(Dense(units = 112))
    model.add(Dense(units = 1))
    model.compile(loss='mean_squared_error', optimizer='Adam', metrics=['mean_absolute_error'])
    model.summary()

    return model

# -------------------------------------#
# main section

data = np.load("TrainingData.npy")
# ignore oxygen deficit
data = data[1:,:,:]
data = np.transpose(data,(2,1,0))
# expect shape = (X,6,11)
np.random.shuffle(data)

# make the dataset smaller
percent = percent/100
endPoint = math.floor(data.shape[0]*percent)
data = data[0:endPoint,:,:]


# for each network
for networkNum in range(11):

    try:
        lstmModel = keras.models.load_model(pathToModelsDir + "model" + str(networkNum) + ".h5")
        print("loaded lstm" + str(networkNum))
    except:
        lstmModel = makeLSTM()
        print("created lstm" + str(networkNum))

    features = data[:,1:,:] # shape = (X, 5, 11)
    labels = data[:,0:1,networkNum] # shape = (X, 1)

    lstmModel.fit(features,labels, validation_split = .1, epochs = 10, batch_size = 1028, verbose = 1)
    lstmModel.save(pathToModelsDir + "/model" + str(networkNum) + ".h5")
