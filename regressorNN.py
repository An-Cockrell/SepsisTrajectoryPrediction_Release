import numpy as np
from keras import backend as K
import random
import os
import keras
import math
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout

percent = 100
# creating a directory to store the dictionaries in
pathToModelsDir = "./"

data = np.load("trainingData.npy")
data = np.transpose(data,(2,1,0))
print(data.shape)
# print(data[0:3,:,:])
np.random.shuffle(data)
# print(data[0:3,:,:])
print("shuffled")

percent = percent/100
endPoint = math.floor(data.shape[0]*percent)
data = data[0:endPoint,:,:]



def makeNN():
    model = Sequential()
    model.add(Dense(units=1000,input_dim=(11)))
    model.add(Dense(units = 1000))

    # model.add(Dropout(0.3))
    model.add(Dense(units = 1500))
    model.add(Dense(units = 1500))

    model.add(Dense(units = 1))

    model.compile(optimizer = 'adam', loss = 'mean_squared_error', metrics = ['mean_absolute_error'])

    return model


try:
    lstmModel = keras.models.load_model(pathToModelsDir + "regressorMLP.h5")
    print("loaded MLP")
except:
    lstmModel = makeNN()
    print("created MLP")


# features = np.load("./trimmedFeatures/" + str(file))
# labels = np.load("./trimmedSingleLabels/" + str(file))
features = data[:,1,1:]
labels = data[:,1,0]


lstmModel.fit(features,labels, validation_split = .1, epochs = 20, batch_size = 5000, verbose = 1)
lstmModel.save(pathToModelsDir + "regressorMLP.h5")
