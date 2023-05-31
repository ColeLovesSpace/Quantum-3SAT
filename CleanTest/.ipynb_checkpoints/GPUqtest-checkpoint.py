# Needed for functions
import time
import numpy as np
from numpy import pi
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import random
from QSAT import *
from Shonings import *

#import packages for classical algorithm
import ctypes
from numpy import ctypeslib as npct

#import packages for quantum algorithm
from qiskit import *
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import qiskit.quantum_info as qi

from qiskit import Aer, transpile
# import qiskit.providers.aer
from qiskit_aer import AerError

from scipy.optimize import minimize


def checkSAT(SAT, state, c, n):
        
    for t in range(c):
        # //get the variables the clause is interested in
        seive = (1 << (abs(int(SAT[t][0])) - 1)) + (1 << (int(abs(SAT[t][1])) - 1)) + (1 << (int(abs(SAT[t][2])) - 1))
        # //get the values in the clause to xor
        gate = 0
            
        if (SAT[t][0] > 0):
            gate += 1 << (abs(int((SAT[t][0]))) - 1)
        if (SAT[t][1] > 0):
            gate += 1 << (abs(int((SAT[t][1]))) - 1)
        if (SAT[t][2] > 0):
            gate += 1 << (abs(int((SAT[t][2]))) - 1)

        a = seive & state
        b = gate ^ a

        if (b == 0):
            return False
    return True

def exhausting(SAT,c,n):
    solutions = []
#     for t in range(c):
#         print((int(SAT[t][0]))," ",(int(SAT[t][1]))," ",(int(SAT[t][2])))
    
    for s in range(0,2**n):
        state = int(s)
        if(checkSAT(SAT, state, c, n)):
            solutions.append(state)
    return solutions
            
    
def setup(n,nc,nSAT):
    par = {}
    par['nQ'] = n
    par['nC'] = nc
    par['nSAT'] = nSAT
    par['sat'] = buildSatSingleSolution(par['nQ'])
    par['nC'] = len(par['sat'])
    par['satList'] = SATset(n,nc,nSAT)
    par['numRepeats'] = 1000
    par['numIt'] = par['nC']
    par['maxIt'] = par['nC']
    
#     par['theta'] = -np.pi/2
#     par['lambda'] = np.pi/2
    
    par['gamma'] = 2.008e+00 
    par['beta'] = -2.395e-04     
   
    
    #Plotting Conditionals
    par['latex'] = True
    par['statevector'] = True
    par['saveEnd'] = False
    par['figName'] = "plots/QAOASAT/{:0>3}.png".format(par['nQ'])
    par['xAsis'] = "$\frac{c}{v}$"
    par['yAsis'] = "$P(\text{solved})_{i<"+"+{}".format(par['maxIt'])+"}$"
    
    par['measure'] = False
    par['had'] = True
    
    par['gpu'] = False
    
    state = ''
    for q in range(n):
        state = state + '0'
    par['state'] = state
    
    return par

par = setup(7,16,1)

solutions = exhausting(par['sat'],par['nC'],par['nQ'])
print(solutions)
print(par['sat'])

def toMin2(phases):
    an = par['nQ']
    probList3 = []
        
    par['gamma'] = phases[0]
    par['beta'] = phases[1]
#     par['numIt']  = phases[2]

    qc2 = buildCircuitCost(par,False)

    if(par['statevector']):
        # Execute and get counts
        DM2 = qi.DensityMatrix.from_instruction(qc2)
        # probs = DM.probabilities_dict()
        DM22 = qi.partial_trace(DM2,[par['nQ']])
        probs2 = DM22.probabilities_dict()
#         entropyList2.append(qi.entropy(DM22))

#     print("Iteration " + str(it) + " took " + str(start-time.time()))

    names2 = list(probs2.keys())
    values2 = list(probs2.values())
#         print(names2)
#         print(values2)

    totalProb2 = 0
#         print(solutions)
    for s in solutions:
#             print(s)
        totalProb2 += probs2[f'{s:>0{an}b}']
#     probList2.append(totalProb2)

    return 1 - totalProb2

maxIt = par['maxIt']#c#n**2

probListGamma = []
probList2 = []
xs = np.linspace(1,maxIt,num = maxIt,dtype='int')
# par['sat'] = par['satList'][0]
# gammaList = np.linspace(0,2*np.pi,num = 6,dtype='float')
ndiv = 7
for g in range(ndiv):
    probListGamma.append([]) 
    
phaseList = []
gammaList = []
for g in range(ndiv):
    par['numIt'] = int(g*maxIt/ndiv)
    gammaList.append(int(g*maxIt/ndiv))
    sendIt = [1.0,1.0]
    res = minimize(toMin2,sendIt,method='COBYLA',options={'disp': True})#, 'maxiter':50
    phaseList.append(res.x)
    print(res.x)
an = par['nQ']
# print(gammaList)
for it in xs:
    par['numIt']  = it
    
#     shots = 10000 

    start = time.time()
    qc2 = buildCircuit(par,False)

    if(par['statevector']):
        # Execute and get counts
        DM2 = qi.DensityMatrix.from_instruction(qc2)
        # probs = DM.probabilities_dict()
        DM22 = qi.partial_trace(DM2,[par['nQ']])
        probs2 = DM22.probabilities_dict()

    names2 = list(probs2.keys())
    values2 = list(probs2.values())

    totalProb2 = 0
    
    for s in solutions:
        totalProb2 += probs2[f'{s:>0{an}b}']
    probList2.append(totalProb2)
    
    for g in range(len(probListGamma)):
        
        par['gamma'] = phaseList[g][0]
        par['beta'] = phaseList[g][1]
        
        
        qc2 = buildCircuitCost(par,False)
#         totalProb = 0
#         for s in solutions:
#             totalProb += list(probs.values())[s]
#         probList.append(totalProb)

        if(par['statevector']):
            # Execute and get counts
            DM2 = qi.DensityMatrix.from_instruction(qc2)
            # probs = DM.probabilities_dict()
            DM22 = qi.partial_trace(DM2,[par['nQ']])
            probs2 = DM22.probabilities_dict()

        

        names2 = list(probs2.keys())
        values2 = list(probs2.values())

        totalProb2 = 0
        for s in solutions:
            totalProb2 += probs2[f'{s:>0{an}b}']
        probListGamma[g].append(totalProb2)
      
    print("Iteration " + str(it) + " took " + str(start-time.time()))
    
for li in range(len(probListGamma)):
    name = "QAOA g="+ str(gammaList[li])
    plt.plot(xs, probListGamma[li] ,label=name)
    
print(len(probList2))
plt.plot(xs, probList2 ,label="QSAT")
    
plt.legend()

n = par['nQ']

name = "plots/clean/6qaOat3st{:0>3}.png".format(n)
plt.savefig(name)