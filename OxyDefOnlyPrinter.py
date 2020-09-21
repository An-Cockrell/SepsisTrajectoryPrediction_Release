# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 05/15/2020
# Script that will print graphs of previously created rolling prediction cones


import matplotlib.pyplot as plt
import numpy as np
import os

colorNum = 0

TrajDataLocation = './OxyDefOnlyTrajectoryData/'

trueData = np.load("finalOxyDefTraj.npy")


files = os.listdir(TrajDataLocation)
files = sorted(files, key = lambda name : int(name[10:-4]))
# print(files)


count = 0
colorNum = -1

plt.figure(figsize=(10,6))
for file in files:
    count += 1
    # if count%5 != 0:
    #     continue

    data = np.load(TrajDataLocation + str(file))
    if data.shape[0] != 1501:
        continue
    print(data.shape)
    print(file)
    colorNum += 1
    if colorNum%5 == 0:
        col = "red"
    elif colorNum%5 == 1:
        col = "magenta"
    elif colorNum%5 == 2:
        col = "purple"
    elif colorNum%5 == 3:
        col = "cyan"
    elif colorNum%5 == 4:
        col = "green"


    plt.title('Oxy Def Trajectory')
    plt.xlabel('Timesteps (7 min)')
    plt.ylabel('Oxydef Level')
    plt.plot(data, color = col)
    plt.plot(trueData[:2500], color='blue')
    plt.xlim((300,2500))
    plt.ylim((0,8200))
    plt.plot()
plt.show()
