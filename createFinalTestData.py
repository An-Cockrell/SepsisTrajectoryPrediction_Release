# NOTE THAT THIS CODE DOES NOT WORK WITH CURRENT VERIONS OF NUMPY, USE 1.14 instead
# This is due to some bug with the numpy/cytpes interface and they arent going to fix it soon
# https://github.com/numpy/numpy/pull/11277

import ctypes
import numpy as np
import sys
from numpy.ctypeslib import ndpointer
import os
import io

import subprocess

_IIRABM = ctypes.CDLL('./IIRABM_Trajectory.so')
_IIRABM.mainSimulation.argtypes = (ctypes.c_float, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                   ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_float))
_IIRABM.mainSimulation.restype = ndpointer(
    dtype=ctypes.c_float, shape=(20, 10000))
    #(numruns*20,10000)

# get class function takes in all the parameters and returns the simulation values with 10000 time steps
def getClass(NIR, IS, oxyHeal, NRI, injNum, seed, internalParameterization):
    answer = 0
    c_float_p = ctypes.POINTER(ctypes.c_float)
    NIR = int(NIR)
    NRI = int(NIR)
    IS = int(IS)
    injNum = int(injNum)
    seed = int(seed)
    NIR = ctypes.c_int(NIR)
    IS = ctypes.c_int(IS)
    NRI = ctypes.c_int(NRI)
    injNum = ctypes.c_int(injNum)
    seed = ctypes.c_int(seed)
    oxyHeal = ctypes.c_float(oxyHeal)
    answer = _IIRABM.mainSimulation(
        oxyHeal, IS, NRI, NIR, injNum, seed, 9,
        internalParameterization.ctypes.data_as(c_float_p))
    return answer

# setting up all the different param combinations
injNum = np.array([5, 10, 15, 20, 25, 30, 35, 40])
NIR = np.array([1, 2, 3, 4])
IS = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
NRI = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
oxyHeal = np.linspace(0.05, 1, 20)
internalParameterization = np.array(
    [1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.float32)




data = [0, 0, 0, 0, 0]
index = 0
iteration = 0

NIR = 2
IS = 4
oxyHeal = .05
NRI = 2
injNum = 27
seed = 1

input = [2, 4, .05, 2, 27, 1]
data = getClass(input[0],input[1],input[2],input[3],input[4], input[5], internalParameterization)


np.save("./finalTestData", data)
np.save("./finalTestParams", input)
