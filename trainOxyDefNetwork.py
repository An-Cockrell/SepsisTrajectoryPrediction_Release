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

makingMLP = False



def PermaDropout(rate):
    return Lambda(lambda x: K.dropout(x, level=rate))

# 5 observations for the network to predict on
numForPredictions = 5

percent = 100
# creating a directory to store the dictionaries in

pathToModelsDir = "./"


# parameters for making the LSTM neural network
def makeLSTM():
    model = Sequential()
    model.add(LSTM(units=1000,return_sequences=True,input_shape=(numForPredictions,1)))
    model.add(LSTM(units=1000,return_sequences=False,input_shape=(numForPredictions,1)))
    model.add(PermaDropout(0.1))

    model.add(Dense(units = 2000))
    model.add(Dense(units = 3000))
    model.add(Dense(units = 1))

    model.compile(loss='mean_squared_error', optimizer='Adam', metrics=['mean_absolute_error'])
    model.summary()

    return model

def makeNN():
    model = Sequential()
    model.add(Dense(units=1000,input_dim=(numForPredictions)))
    model.add(PermaDropout(.01))
    model.add(Dense(units = 150))
    model.add(Dense(units = 150))

    model.add(Dense(units = 1))

    model.compile(optimizer = 'adam', loss = 'mean_squared_error', metrics = ['mean_absolute_error'])
    model.summary()

    return model


# -------------------------------------#
# main section

data = np.load("trainingData.npy")
data = np.transpose(data,(2,1,0))
# print(data[0:3,:,:])
np.random.shuffle(data)
# print(data[0:3,:,:])
print("shuffled")

percent = percent/100
endPoint = math.floor(data.shape[0]*percent)
data = data[0:endPoint,:,:]


if makingMLP:
    try:
        lstmModel = keras.models.load_model(pathToModelsDir + "OxyDefOnlyMLP.h5")
        print("loaded MLP Model")
    except:
        lstmModel = makeNN()
        print("created NN")
else:
    try:
        lstmModel = keras.models.load_model(pathToModelsDir + "OxyDefOnlyLSTMp1.h5")
        print("loaded LSTM Model")
    except:
        lstmModel = makeLSTM()
        print("created LSTM")

# features = np.load("./trimmedFeatures/" + str(file))
# labels = np.load("./trimmedSingleLabels/" + str(file))
features = data[:,1:,0]
labels = data[:,0:1,0]
# print(labels.shape)
# print(features.shape)
temp = math.floor(features.shape[0]/10)

lstmModel.fit(features,labels, validation_split=.1, epochs = 20, batch_size = 4000, verbose = 1)
if makingMLP:
    lstmModel.save(pathToModelsDir + "OxyDefOnlyMLP.h5")
else:
    lstmModel.save(pathToModelsDir + "OxyDefOnlyLSTMp1.h5")
