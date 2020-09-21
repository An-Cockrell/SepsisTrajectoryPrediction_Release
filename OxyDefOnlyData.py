# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 06/09/2020
# Script that will generate random Oxydef Trajectories using 5 sample input of oxydef by means of a MLP NN

import numpy as np
import math
import random
import tensorflow as tf
import keras
# import matplotlib.pyplot as plt
import os
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False #suppress depreciation warnings


modelName = "./OxyDefOnlyMLP.h5"
import h5py
f = h5py.File(modelName,'r+')
data_p = f.attrs['training_config']
data_p = data_p.decode().replace("learning_rate","lr").encode()
f.attrs['training_config'] = data_p
f.close()



print("\n\nOxyDefOnly.py\n\n")

startNum = 1000

numToPredict = 1000
numInCone = 100
timeSinceLastResample = 10
resamples = 1
TrajDir = "./OxyDefOnlyFullTraj/"


# --------METHODS--------
# this function will take sequential data and make the (X, 6, 12) data used for prediction
# newest data point is at index [:,0,:], oldest at index [:,5,:]
def makeFeatures(currentData):
    features = np.ones((currentData.shape[0],6, currentData.shape[1]))
    features *= -1
    index = 0
    for j in range(7, currentData.shape[0]):
        feature = np.zeros((1,6, currentData.shape[1]))
        for k in reversed(range(6)):
            temp = currentData[j - k, :]
            feature[0,k, :] = temp

        features[index,:,:] = feature
        index += 1

    return features



try:
    os.mkdir(TrajDir)
except:
    pass



# --------MAIN--------
data = np.load("finalTestData.npy")
print("Data shape: " + str(data.shape))
# expect (20,X)
data = data[[0],:]
traj = data[0,:3000]
data = data.T
data = makeFeatures(data)

label = data[:,0,:]
features = data[:,1:,:]
print(features.shape)

model = keras.models.load_model(modelName)

startPoint = startNum
# int is the index of the point in time
# plt.figure(figsize=(10,6))
for i in range(resamples):
    dataToSave = np.zeros((startNum+numToPredict+1,numInCone))
    for z in range(numInCone):
        int = startNum
        # creating inputs to regressor model to predict oxydef
        # pre allocating space for speed, shape should be (X,5)
        inputs = np.zeros((numToPredict + 1,features.shape[1],1))

        # Start predicting from real observations
        for k in range(1):
            # j is the 5 observation time steps
            for j in range(5):
                # use the correct cytokine model to predict cytokine value for next time step
                predict = model.predict(features[int:int+1,:,0])
                # add it to inputs at the correct index
                inputs[k,j,0] = predict
                # move forward one time step
                int += 1
        # print(inputs[:10,:,0])
        # start predicting from previous predictions
        for k in range(numToPredict):
            print(z,k,end="\r")
            for j in range(1):
                predict = model.predict(inputs[k:k+1,:,0])
                # insert the prediction into the first observation slot in the next time step
                inputs[k+1,j,0] = predict
            # insert the 4 most recent predictions into the next 4 observation slots in the next time step
            inputs[k+ 1,-4:,:] = inputs[k,:4,:]
        # print(inputs[:10,:,0])

        # combine the true cytokines and the predicted cytokines
        inputs = np.vstack((features[:startNum,:,:],inputs))
        dataToSave[:,z] = inputs[:,0,0]
    #     plt.plot(traj)
    # plt.plot(dataToSave,color="red")
    print(dataToSave.shape)

    startNum += timeSinceLastResample
    np.save(TrajDir + "Trajectory" + str(i) + ".npy", dataToSave)

# plt.title("MLP Oxydef Only")
# plt.xlabel("Timesteps")
# plt.ylabel("Oxygen Deficit")
# plt.ylim(-900,10000)
# plt.xlim(-300,4000)
# plt.show()
