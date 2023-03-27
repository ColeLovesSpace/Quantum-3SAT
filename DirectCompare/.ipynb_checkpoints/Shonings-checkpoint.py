# Needed for functions
import time
import numpy as np
from numpy import pi
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import random

#import packages for classical algorithm
import ctypes
from numpy import ctypeslib as npct

# Stochastic algorithm

#initialize variables for classical solver
array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.int_, ndim=1, flags='CONTIGUOUS')

c_lib_numpy = npct.load_library("lib/SATProbCalc.so", ".")

c_lib_numpy.SolveSATbpp.restype = ctypes.c_int #ctypes.POINTER(ctypes.c_double * (10**2))
c_lib_numpy.SolveSATbpp.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int]

c_lib_numpy.Shoning.restype = ctypes.c_int #ctypes.POINTER(ctypes.c_double * (10**2))
c_lib_numpy.Shoning.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int]

def SolveSATbpp(n, sat, clauses, iterations):
    steps = c_lib_numpy.SolveSATbpp(n, sat.flatten(), clauses, iterations)
    return steps

def Shoning(n, sat, clauses, iterations):
    steps = c_lib_numpy.Shoning(n, sat.flatten(), clauses, iterations)
    return steps

# def classicalStatistics(n, SAT, c, i, numRuns):
#     stepList = [] 
#     for x in range(numRuns):
#         stepList.append(SolveSATbpp(n, SAT, c, i))
#         print(stepList[-1])
# #     counter = Counter(stepList)
# #     res = {item: count for item, count in counter.items()}
# #     print(res)
#     return stepList

# def classicalStatisticsS(n, SAT, c, i, numRuns):
#     stepList = []
#     for x in range(numRuns):
#         stepList.append(Shoning(n, SAT, c, i))
#         print(stepList[-1])
# #     counter = Counter(stepList)
# #     res = {item: count for item, count in counter.items()}
# #     print(res)
#     return stepList

def classicalSolve(n,SAT,c,i,numIt):

    SAT = np.array(SAT).astype('double')
#     for x in range(len(SAT)):
#         SAT[x] = -SAT[x]

    start = time.time()
    stepList = []
    for x in range(numIt):
        stepList.append(SolveSATbpp(n, SAT, c, i))
#         print(stepList[-1])
#     counter = Counter(stepList)
#     res = {item: count for item, count in counter.items()}
#     print(res)
    print("C++ time to solve: " , start - time.time())

    return stepList

def classicalSolveS(n,SAT,c,i,numIt):

    SAT = np.array(SAT).astype('double')
#     for x in range(len(SAT)):
#         SAT[x] = -SAT[x]

    start = time.time()
    
    stepList = []
    for x in range(numIt):
        stepList.append(Shoning(n, SAT, c, i))
#         print(stepList[-1])
#     res = classicalStatisticsS(n, SAT, c, i, numIt)
    
    
    print("C++ time to solve: " , start - time.time())

    return stepList

# def classicalSolve(n,SAT,c,i,numIt):

#     SAT = np.array(SAT).astype('double')
# #     for x in range(len(SAT)):
# #         SAT[x] = -SAT[x]
#     # print(SAT)

#     start = time.time()
#     res = classicalStatistics(n, SAT, c, i, numIt)
#     print("C++ time to solve: " , start - time.time())

#     return res

## Shared functions/dependencies 

# Build 3SAT instance with only one solution
def buildSatSingleSolution(numVar):
    sat = [[-1,-2,-3],[-1,2,3],[1,-2,3],[1,2,-3],[-1,-2,3],[-1,2,-3],[1,-2,-3]]
    for n in range(4,numVar+1):
        sat.append([2,3,n])
        sat.append([-2,-3,-n])
        sat.append([-2,-3,n])
        
    return sat

def SATset(n,nc,ns):
    SATlist = []
    variables = list(range(n+1))
    sign = [-1,1]
    variables.remove(0)
    
    while len(SATlist)<ns:
        SAT=[]
        while len(SAT)<nc:
            clause = list(set(random.choices(variables,k=3)))
            for v in range(len(clause)):
                clause[v]=clause[v]*random.choice(sign)
            clause = sorted(clause)
            if len(clause) == 3 and clause not in SAT:
                SAT.append(clause)
        SAT = sorted(SAT)
        if SAT not in SATlist:
            SATlist.append(SAT)
    return SATlist