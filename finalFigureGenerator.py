# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 05/15/2020
# Script that will generate rolling prediction cones

import numpy as np
import math
import random
import keras
import os
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False #suppress depreciation warnings

print("\n\nfinalFigGenerator.py\n\n")

SAVE_FIG = True
PRINT_FIG = False

SAVE_DATA = True

# number of trajectories for the cloud
numTrajectories = 100

# number of different clouds to make
timesToResample = 300

# time steps inbetween clouds
timeSinceLastResample = 10
startNum = 100
numToPredict = 100

col = 'red'

ModelDir = "./inputPredictionParallelMLPModels/"

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


# --------Creating Directories--------
if SAVE_DATA:
    try:
        TrajDir ="./TrajectoryDataFinal/"
        os.mkdir(TrajDir)
        print("Created Trajectory Dir")
    except:
        pass

# --------MAIN--------
data = np.load("finalTestData.npy")
print("Data shape: " + str(data.shape))
# expect (20,X)
data = data[[0,2,3,4,5,12,13,14,15,16,17,18],:]
traj = data[0,:]
# traj.shape = X
data = data[1:,:]
data = data.T
# data.shape = (X,11)



data = makeFeatures(data)
# data.shape = (X,6,11)


features = data[:,1:,:]

print("Features shape: " + str(features.shape))
# expect (X, 5, 11)

# ---LOAD MODELS---
# model to regress oxyDef after predicting cytokines
regressor = keras.models.load_model("Models/regressorMLP.h5")

# loading all the cytokine models into a vector
models = []
for i in range(11):
    tempModel = keras.models.load_model(ModelDir + "model" + str(i) + ".h5")
    models.append(tempModel)

print("models loaded")
print("startNum:  " + str(startNum))
resampleStart = startNum%10000
originalStart = startNum%10000


# for the number of prediction clouds specified
for x in range(timesToResample):
    # prediction start index for this cloud
    graphStart = startNum%10000
    # preallocating memory for speed
    featuresToSave = np.zeros((numToPredict + graphStart, features.shape[2] + 1, numTrajectories)) #(NumPredictions, cytokines+Oxy, numTrajectories)


    # for each individual trajectory in the cloud
    for z in range(numTrajectories):

        # int is the index of the current point in time
        int = startNum -1

        # creating inputs to regressor model to predict oxydef
        # pre allocating space for speed, shape should be (X,5,11)
        inputs = np.zeros((numToPredict + 1,features.shape[1], features.shape[2]))

        # Start predicting from real observations
        for k in range(1):
            # j is the 5 observation time steps
            for j in range(5):
                # i is each cytokine
                for i in range(11):
                    # use the correct cytokine model to predict cytokine value for next time step
                    predict = models[i].predict(features[int:int+1,:,:])
                    # add it to inputs at the correct index
                    inputs[k,-j,i] = predict
                # move forward one time step
                int += 1

        # start predicting from previous predictions
        for k in range(numToPredict):
            for j in range(1):
                for i in range(11):
                    predict = models[i].predict(inputs[k:k+1,:,:])
                    # insert the prediction into the first observation slot in the next time step
                    inputs[k+1,j,i] = predict
                # insert the 4 most recent predictions into the next 4 observation slots in the next time step
                inputs[k+1,-4:,:] = inputs[k,:4,:]

        # combine the true cytokines and the predicted cytokines
        inputs = np.vstack((features[startNum - graphStart:startNum,:,:],inputs))

        # regress oxydef from the inputs
        predictions = regressor.predict(inputs[:,0,:])

        # for each prediction
        for a in range(inputs.shape[0] - 1):
            # for each cytokine
            for b in range(11):
                featuresToSave[a,b,z] = inputs[1 + a,0,b]
            featuresToSave[a,-1,z] = predictions[a]
            # (number of samples, [cytokines,oxydef], trajectory)


    print("Predicted Trajectories Shape: " + str(featuresToSave.shape))

    if SAVE_DATA:
        np.save(TrajDir + "/Trajectory" + str(x) + ".npy", featuresToSave)
    # increment the start point by the amount indicated
    startNum += timeSinceLastResample
# print graphs
plt.show()
