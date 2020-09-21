Quick explanation of all the scripts and their general function

Use the following command to compile the c++ into an executable simulation
g++ -fPIC -Wall -shared -o IIRABM_Trajectory.so ./iirabm_c++/*.cpp -std=c++11 -O30

To generate a full data set, run trajSweepNonMPI.py. Currently set to run with 100 processors, can be easily changed to run on a single processor. Final data set will be in the hundreds of gigabytes. Each individual file will be the fulloutput of 150 iirabm runs. 

createMasterDataSet.py will combine all the .csv files into a training data set and a set of combined data. The training data will have 6 sequential 12D points, combined data would only be a single 12D point.
createFinalTestData.py will create a data file of a single run as the final comparison data.

regressorNN.py will create either an LSTM or MLP NN to regress oxygen deficit from the 11 monitored cytokines
graphRegressor.py will graph the regressor NN to compare to true data.

inputPredictionPermaDrop.py will create 11 LSTM networks that will predict the next value for each individual cytokine based on the state of the previous 11 cytokine values while including stocastic variations from a permanent dropout layer.
inputPredictionParallel.py does the same thing as inputPredictionPermaDrop.py, but is set to run on 11 processors and uses an MLP instead.

finalFigureGenerator.py will create the predicted final trajectory cones based on the finalTestData, using the regressor and inputPrediction networks.
finalFigurePrinter.py displays the predicted trajectories.

trainOxyDefNetwork.py creates a NN that will predict oxygen deficit values based only on previous oxygen deficit values.
OxyDefOnlyData.py uses the oxydef NN to create oxy def trajectories.
OxyDefOnlyPRint.py displays these trajectories
