import time
import numpy as np
import ctypes
from numpy import ctypeslib as npct
import matplotlib.pyplot as plt
# from pathos.multiprocessing import ProcessingPool as Pool
# from bppSAT import *

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.int_, ndim=1, flags='CONTIGUOUS')

# numpy.ctypeslib.as_array(obj, shape=None)
x = 250#75
i = 1844674407370955161#int(np.ceil((x**2))) 
c = 1065 #375

c_lib_numpy = npct.load_library("lib/SATProbCalc.so", ".")
# c_lib_numpy.SolveSAT.restype = ctypes.c_double#ctypes.POINTER(ctypes.c_double * (10**2))
# c_lib_numpy.SolveSAT.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int, ctypes.c_int]
# c_lib_numpy.SolveSATactual.restype = npct.ndpointer(dtype=ctypes.c_double, shape=(i,))
# c_lib_numpy.SolveSATactual.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int, ctypes.c_int]
# c_lib_numpy.SolveSATiterations.restype = ctypes.c_int
# c_lib_numpy.SolveSATiterations.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int, ctypes.c_int]
# c_lib_numpy.SolveSATgates.restype = npct.ndpointer(dtype=ctypes.c_double, shape=(c,))
# c_lib_numpy.SolveSATgates.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_int, ctypes.c_int]
c_lib_numpy.SolveSATbpp.restype = None#ctypes.c_double#ctypes.POINTER(ctypes.c_double * (10**2))
c_lib_numpy.SolveSATbpp.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_int, ctypes.c_ulonglong]

def buildSatSingleSolution(numVar):
    sat = [[1,2,3],[-1,2,3],[1,-2,3],[1,2,-3],[-1,-2,3],[-1,2,-3],[1,-2,-3]]
    for n in range(4,numVar+1):
        sat.append([2,3,n])
        sat.append([2,3,-n])
        sat.append([-2,-3,n])

    return sat
# def SolveSAT(n, sat, clauses, iterations, variables):
#     probs = c_lib_numpy.SolveSAT(n, sat.flatten(), clauses, iterations, variables)
#     return probs
# def SolveSATactual(n, sat, clauses, iterations, variables):
#     probs = c_lib_numpy.SolveSATactual(n, sat.flatten(), clauses, iterations, variables)
#     return probs
# def SolveSATiterations(n, sat, clauses, iterations, variables):
#     probs = c_lib_numpy.SolveSATiterations(n, sat.flatten(), clauses, iterations, variables)
#     return probs
# def SolveSATgates(n, sat, clauses, iterations, variables):
#     probs = c_lib_numpy.SolveSATgates(n, sat.flatten(), clauses, iterations, variables)
#     return probs
def SolveSATbpp(n, sat, clauses, iterations):
    c_lib_numpy.SolveSATbpp(n, sat.flatten(), clauses, iterations)
    return 5
# def PythonSolve():

# def getInstance():
#     instance = []
#     with open('/Users/colecoughlin/Desktop/QSAT/instance2.txt', 'r', encoding='utf-8') as g:
#         data = g.readlines()
#         for line in data:
#             line = line.strip().split(' ')
#             for var in line:
#                 instance.append(-int(var.strip()))
#     return instance



def getSAT(s):
    filename = "./SATS/UF250.1065.100/uf250-0" + str(s + 11) + ".cnf"#"./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" ./SATS/UF75.325.100/uf75-024.cnf
    #filename = "./SATS/UF75.325.100/uf75-0" + str(s + 75) + ".cnf"#"./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" ./SATS/UF75.325.100/uf75-024.cnf

    print(filename)

    with open(filename) as f:
        lines = f.readlines()

    dets = lines[7].strip().split(" ")
    n = int(dets[2])
    c = int(dets[4])
    lines = lines[8:-3]

    print(n, " ", c)

    SAT = []
    for line in lines:
        SAT.append(line.strip().replace(" 0","").split(" "))

    SAT = np.array(SAT).astype('double')
    for x in range(len(SAT)):
        SAT[x] = -SAT[x]
    # print(SAT)

    start = time.time()
    probs = SolveSATbpp(n, SAT, c, i + 1)
    print("C++ time to solve: " , start - time.time())

    # SAT = np.array(SAT).astype('integer')

    # start = time.time()
    # probs = PythonSolve(n, SAT, c, i+1, 3*c)#int(n*n/2)
    # print("Python time to solve: " + start - time.time())

    

    # w = n**2
    # probs = SolveSATgates(n, SAT, c, w, 3 * c)  # int(n*n/2)

    return probs

def getSATprime(s):
    # filename = "./SATS/UF250.1065.100/uf250-0" + str(s + 70) + ".cnf"#"./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" ./SATS/UF75.325.100/uf75-024.cnf
    filename = "./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" #./SATS/UF75.325.100/uf75-024.cnf

    print(filename)

    with open(filename) as f:
        lines = f.readlines()

    dets = lines[0].strip().split(" ")
    n = int(dets[2])
    c = int(dets[3])
    lines = lines[1:]

    print(n, " ", c, " numLines: ", len(lines))

    SAT = []
    for line in lines:
        SAT.append(line.strip().replace(" 0","").split(" "))

    SAT = np.array(SAT).astype('double')
    for x in range(len(SAT)):
        SAT[x] = -SAT[x]
    # print(SAT)

    start = time.time()
    probs = SolveSATbpp(n, SAT, c, i + 1)
    print("C++ time to solve: " , start - time.time())

    return probs

def average():


    probs = getSAT(1)


    return probs


probs = average()
