import numpy as np
import keras
import matplotlib.pyplot as plt


data = np.load("finalTestData.npy")
print("Data shape: " + str(data.shape))
# expect (20,X)
data = data[[0,2,3,4,5,12,13,14,15,16,17,18],:]
traj = data[0,:]
# traj.shape = (X)
data = data[1:,:]
data = data.T
# data.shape = (X,11)


model = keras.models.load_model("./regressorMLP.h5")
predictedTraj = model.predict(data)

plt.figure(figsize=(10,6))

plt.plot(traj,c="blue")
plt.plot(predictedTraj, c="red")
plt.show()
