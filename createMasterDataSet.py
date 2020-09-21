import numpy as np
import os
from sklearn.preprocessing import normalize

combinedSamples = 11000000
trainingSamples = 10000000

files = os.listdir("/sepsisData/")
np.random.shuffle(files)

numForPredictions = 5

combinedData = np.zeros((12,combinedSamples))
trainingData = np.zeros((12,numForPredictions+1,trainingSamples))
combinedIndex = 0
trainingIndex = 0

lowerLimit = 300
upperLimit = 3000


numRandomPerTraj = 200


fullFlag = False
trainingFull = False
count = 0
fullFileCount = 0
for file in files:
    count += 1
    print(f'{"Files Read: " + str(count) + "/" + str(len(files))}\r', end = "")
    if file.startswith("Data") and not fullFlag:
        try:
            # print(str(file))
            currentData = np.load("/sepsisData/" + str(file))
            currentData = currentData[[0,2,3,4,5,12,13,14,15,16,17,18],lowerLimit:,:] #12,10000,150
            # print(currentData)
            summedData = np.sum(currentData[1:,:,:], 0)
            # print(summedData.shape)

            for i in range(summedData.shape[1]):
                stopIndex = numForPredictions+1
                while summedData[stopIndex, i] > 1 and stopIndex < upperLimit:
                    stopIndex += 1
                if stopIndex == upperLimit:
                    fullFileCount += 1
                if stopIndex > numForPredictions+1:
                    # print(currentData[:,stopIndex-2:stopIndex+1,i], stopIndex)
                    for f in range(numRandomPerTraj):
                        k = np.random.randint(numForPredictions, stopIndex)
                        if not trainingFull:
                            for j in range(0,numForPredictions+1):
                                trainingData[:,j,trainingIndex] = currentData[:,k-j,i]
                            trainingIndex += 1

                        for j in range(0, combinedData.shape[0]):
                            combinedData[j,combinedIndex] = currentData[j,k,i]
                        combinedIndex += 1

        except OSError as e:
            os.rename("sepsisData/" + str(file), "/CorruptFiles/" + str(file))
            print("moved " + str(file) + " to /CorruptFiles")

        except IndexError as e:
            try:
                print(index)
                print(trainingData[:,:, index-2:index+2])
            except:
                pass
        if trainingIndex >= trainingSamples and not trainingFull:
            trainingFull = True
            print("Training Full            ")
            print("Files included in Training: " + str(count) + "/"+ str(len(files)) + "\n")
        if combinedIndex >= combinedSamples:
            fullFlag = True
            print("Combined Full            \n")
            break

print("Total Files Read: " + str(count) + "/" + str(len(files)) + "\n")
print("Num Trajectories Clipped: " + str(fullFileCount))
# print(combinedData[:,-10:])
# print(trainingData[:,:,-10:])

print("TrainingData shape: " + str(trainingData.shape))

combinedData = combinedData[:,:combinedIndex]
trainingData = trainingData[:,:,:trainingIndex]

np.save("./combinedData.npy", combinedData)
np.save("./trainingData.npy", trainingData)
