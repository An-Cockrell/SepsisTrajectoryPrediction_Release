import warnings
import numpy as np
import matplotlib.pyplot as plt
import random



SUPPRESS_WARNINGS = True
if SUPPRESS_WARNINGS:
    warnings.filterwarnings("ignore")




data = np.load("./combinedData.npy")
print("loaded")
print(data.shape)
# data = np.load("./Data0_25.npy")
# data = data[[0,2,3,4,5,12,13,14,15,16,17,18],:]
# data = np.reshape(data, (12, data.shape[1] * data.shape[2]))
data = data[:,:8000000]
# data = data[np.random.choice(data.shape[1], 200000, replace=False), :]
# data = random.sample(data,data.shape[1]/10)
print("Subset made")
print(data.shape)
data = np.round(data, 0)
labels = data[0,:]

feats = data[1:,:]

summedCytokines = np.sum(feats, 0)


labels = np.array(labels)
summedCytokines = np.array(summedCytokines)
summedDict = {}
Stats = np.zeros((7, labels.shape[0]))


for i in range(labels.shape[0]):
    try:
        summedDict[str(summedCytokines[i])].append(labels[i])
    except KeyError:
        summedDict[str(summedCytokines[i])] = [labels[i]]

# print(sorted(summedDict.keys()))

index = 0
for key in sorted(summedDict.keys()):
    if float(key) > 1:
        data = summedDict[key]

        lower = 0
        upper = 0

        Q1 = np.quantile(data, .25)
        Q2 = np.quantile(data, .5)
        Q3 = np.quantile(data, .75)

        mean = np.mean(data)
        for j in range(len(data)):
            if data[j] > mean:
                lower = data[0:j-1]
                upper = data[j-1:]
                break
        lowerMean = np.mean(lower)
        upperMean = np.mean(upper)

        Stats[0,index] = key
        Stats[1,index] = Q1
        Stats[2,index] = Q2
        Stats[3,index] = Q3
        Stats[4,index] = mean
        Stats[5,index] = upperMean
        Stats[6,index] = lowerMean

        index += 1

index = 0

print("plotting")
# print(Stats[0][0:30])
# print(Stats[4][0:30])
fig = plt.figure(figsize = (5,3))
ax1 = fig.add_subplot(111)


ax1.scatter(summedCytokines, labels, c='k', marker="o", s = 1, alpha = .005)
ax1.scatter(Stats[0], Stats[1], c='c', marker="o", s = 1, label="quartile") #1quartile
ax1.scatter(Stats[0], Stats[2], c='g', marker="o", s = 1, label="median") #median
ax1.scatter(Stats[0], Stats[3], c='c', marker="o", s = 1) #3quartile
ax1.scatter(Stats[0], Stats[4], c='r', marker="o", s = 1, label = "mean") #mean
# ax1.scatter(Stats[0], Stats[5], c='c', marker="o", s = 1, label = "mean of upper half") #mean upper half
# ax1.scatter(Stats[0], Stats[6], c='c', marker="o", s = 1, label = "mean of lower half") # mean lower half
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2)
ax1.tick_params(direction='in', length=6, width=2)

ax1.legend(loc="lower right",prop={'size': 10, 'weight':'bold'}, markerscale=4)


# plt.title('Cytokine Total vs Oxygen Deficit')
ax1.set_xlabel('Cytokine Total',fontweight='bold')
ax1.set_ylabel('Oxygen Deficit (arb. units)',fontweight='bold')
plt.savefig("./cytokineOxyDefError.tif", dpi=400,bbox_inches='tight')
# plt.show()




# plt.show()
