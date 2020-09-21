# Written by Dale Larie, send questions to Dale.B.Larie@gmail.com
# Last updated 05/15/2020
# Script that will print graphs of previously created rolling prediction cones


import matplotlib.pyplot as plt
import numpy as np
import os

colorNum = 0

TrajDataLocation = './TrajectoryDataFinal/'

trueData = np.load("finalTestData.npy")
trueData = trueData[[0,2,3,4,5,12,13,14,15,16,17,18], :]



trueCytokines = trueData[1:, :]
trueOxy = trueData[0,:]

# putting oxydef at the end
trueData = np.vstack((trueCytokines,trueOxy))
files = os.listdir(TrajDataLocation)
# files = sorted(files, key = lambda name : int(name[10:-4]))
# print(files)

#
first = plt.figure(figsize=(10,6))
second = plt.figure(figsize=(10,6))
third = plt.figure(figsize=(10,6))
fourth = plt.figure(figsize=(10,6))
fifth = plt.figure(figsize=(10,6))
sixth = plt.figure(figsize=(10,6))
seventh = plt.figure(figsize=(10,6))
eigth = plt.figure(figsize=(10,6))
nineth = plt.figure(figsize=(10,6))
tenth = plt.figure(figsize=(10,6))
eleventh = plt.figure(figsize=(10,6))
twelth = plt.figure(figsize=(10,6))
# plt.hold()
count = 0
for file in files:
    # if count%5 != 0:
    #     continue
    count += 1
    print(file)
    data = np.load(TrajDataLocation + str(file))
    print(data.shape)

    # for i in range(data.shape[1]):
    #     data[:,i,:] = (data[:,i,:] * (normParams[i,1] - normParams[i,0])) + normParams[i,0]
    # print(data[:,-1,0])
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
    #

    ax1 = first.add_subplot(1, 1, 1)
    ax1.set_title('Cytokine Trajectory')
    ax1.set_xlabel('Timesteps (7 min)')
    ax1.set_ylabel('Cytokine Level')

    ax2 = second.add_subplot(1, 1, 1)
    ax2.set_title('Cytokine Trajectory')
    ax2.set_xlabel('Timesteps (7 min)')
    ax2.set_ylabel('Cytokine Level')

    ax3 = third.add_subplot(1, 1, 1)
    ax3.set_title('Cytokine Trajectory')
    ax3.set_xlabel('Timesteps (7 min)')
    ax3.set_ylabel('Cytokine Level')

    ax4 = fourth.add_subplot(1, 1, 1)
    ax4.set_title('Cytokine Trajectory')
    ax4.set_xlabel('Timesteps (7 min)')
    ax4.set_ylabel('Cytokine Level')

    ax5 = fifth.add_subplot(1, 1, 1)
    ax5.set_title('Cytokine Trajectory')
    ax5.set_xlabel('Timesteps (7 min)')
    ax5.set_ylabel('Cytokine Level')

    ax6 = sixth.add_subplot(1, 1, 1)
    ax6.set_title('Cytokine Trajectory')
    ax6.set_xlabel('Timesteps (7 min)')
    ax6.set_ylabel('Cytokine Level')

    ax7 = seventh.add_subplot(1, 1, 1)
    ax7.set_title('Cytokine Trajectory')
    ax7.set_xlabel('Timesteps (7 min)')
    ax7.set_ylabel('Cytokine Level')

    ax8 = eigth.add_subplot(1, 1, 1)
    ax8.set_title('Cytokine Trajectory')
    ax8.set_xlabel('Timesteps (7 min)')
    ax8.set_ylabel('Cytokine Level')

    ax9 = nineth.add_subplot(1, 1, 1)
    ax9.set_title('Cytokine Trajectory')
    ax9.set_xlabel('Timesteps (7 min)')
    ax9.set_ylabel('Cytokine Level')

    ax10 = tenth.add_subplot(1, 1, 1)
    ax10.set_title('Cytokine Trajectory')
    ax10.set_xlabel('Timesteps (7 min)')
    ax10.set_ylabel('Cytokine Level')

    ax11 = eleventh.add_subplot(1, 1, 1)
    ax11.set_title('Cytokine Trajectory')
    ax11.set_xlabel('Timesteps (7 min)')
    ax11.set_ylabel('Cytokine Level')



    ax1.plot(data[:,0,:], color=col)
    ax1.plot(trueData[0,7:data.shape[0]], color='blue')
    ax2.plot(data[:,1,:], color=col)
    ax2.plot(trueData[1,7:data.shape[0]], color='blue')
    ax3.plot(data[:,2,:], color=col)
    ax3.plot(trueData[2,7:data.shape[0]], color='blue')
    ax4.plot(data[:,3,:], color=col)
    ax4.plot(trueData[3,7:data.shape[0]], color='blue')
    ax5.plot(data[:,4,:], color=col)
    ax5.plot(trueData[4,7:data.shape[0]], color='blue')
    ax6.plot(data[:,5,:], color=col)
    ax6.plot(trueData[5,7:data.shape[0]], color='blue')
    ax7.plot(data[:,6,:], color=col)
    ax7.plot(trueData[6,7:data.shape[0]], color='blue')
    ax8.plot(data[:,7,:], color=col)
    ax8.plot(trueData[7,7:data.shape[0]], color='blue')
    ax9.plot(data[:,8,:], color=col)
    ax9.plot(trueData[8,7:data.shape[0]], color='blue')
    ax10.plot(data[:,9,:], color=col)
    ax10.plot(trueData[9,7:data.shape[0]], color='blue')
    ax11.plot(data[:,10,:], color=col)
    ax11.plot(trueData[10,7:data.shape[0]], color='blue')


    ax12 = twelth.add_subplot(1, 1, 1)
    ax12.set_title('oxyDef Trajectory')
    ax12.set_xlabel('Timesteps (7 min)')
    ax12.set_ylabel('Oxydef Level')
    ax12.plot(data[:,-1,:], color = col)
    ax12.plot(trueData[11,6:data.shape[0]], color='blue')
    ax12.set_ylim([0,10000])

    first.canvas.draw()
    second.canvas.draw()
    third.canvas.draw()
    fourth.canvas.draw()
    fifth.canvas.draw()
    sixth.canvas.draw()
    seventh.canvas.draw()
    eigth.canvas.draw()
    nineth.canvas.draw()
    tenth.canvas.draw()
    eleventh.canvas.draw()
    twelth.canvas.draw()


plt.show()
