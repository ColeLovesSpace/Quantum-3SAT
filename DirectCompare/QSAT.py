# Needed for functions
import time
import numpy as np
from numpy import pi
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import random

#import packages for quantum algorithm
from qiskit import *
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import qiskit.quantum_info as qi

from qiskit import Aer, transpile
# import qiskit.providers.aer
from qiskit_aer import AerError

# Quantum algorithm

def OR(a,b,c,d,e,circuit):
    circuit.x(a)
    circuit.x(b)
    circuit.x(c)
    circuit.ccx(a,b,d)
    circuit.ccx(c,d,e)
    circuit.x(a)
    circuit.x(b)
    circuit.x(c)
    circuit.x(e)
    circuit.reset(d)

def AND(a,b,c,d,e,circuit):
    circuit.ccx(a,b,d)
    circuit.ccx(c,d,e)
    circuit.reset(d)
    
def had3(e,a,b,c,circuit,par):
    circuit.ch(e,abs(a)-1)
    circuit.ch(e,abs(b)-1)
    circuit.ch(e,abs(c)-1)
    
def not3(e,a,b,c,circuit,par):
    circuit.cnot(e,abs(a)-1)
    circuit.cnot(e,abs(b)-1)
    circuit.cnot(e,abs(c)-1)
    
def hadAll(e,a,b,c,circuit,par):
    for b in range(par['nQ']):
        circuit.ch(e,b)
#     circuit.ch(e,abs(a)-1)
#     circuit.ch(e,abs(b)-1)
#     circuit.ch(e,abs(c)-1)
    
        
def remStates(a,b,c,d,circuit):
    if a < 0:
        circuit.x(abs(a)-1)
    if b < 0:
        circuit.x(abs(b)-1)
    if c < 0:
        circuit.x(abs(c)-1)
    
    circuit.rcccx(abs(a)-1,abs(b)-1,abs(c)-1,d)
        
    if a < 0:
        circuit.x(abs(a)-1)
    if b < 0:
        circuit.x(abs(b)-1)
    if c < 0:
        circuit.x(abs(c)-1)
         
def endRem(a,b,c,d,e,circuit):
    if a < 0:
        circuit.x(abs(a)-1)
    if b < 0:
        circuit.x(abs(b)-1)
    if c < 0:
        circuit.x(abs(c)-1)
        
def qft3(c,q1,q2,q3):
    c.h(abs(q3)-1)
    c.cp(pi/4, abs(q1)-1, abs(q3)-1)
    c.cp(pi/2, abs(q2)-1, abs(q3)-1)
    c.h(abs(q2)-1)
    c.cp(pi/2, abs(q1)-1, abs(q2)-1)
    c.h(abs(q1)-1)
    c.swap(abs(q1)-1,abs(q3)-1)
    
def qft3i(c,q1,q2,q3):
    c.swap(abs(q1)-1,abs(q3)-1)
    c.h(abs(q1)-1)
    c.cp(-pi/2, abs(q1)-1, abs(q2)-1)
    c.h(abs(q2)-1)
    c.cp(-pi/2, abs(q2)-1, abs(q3)-1)
    c.cp(-pi/4, abs(q1)-1, abs(q3)-1)
    c.h(abs(q3)-1)

def qft_rotations(circuit, n):
    """Performs qft on the first n qubits in circuit (without swaps)"""
    if n == 0:
        return circuit
    n -= 1
    circuit.h(n)
    for qubit in range(n):
        circuit.cp(pi/2**(n-qubit), qubit, n)
    # At the end of our function, we call the same function again on
    # the next qubits (we reduced n by one earlier in the function)
    qft_rotations(circuit, n)
    
def swap_registers(circuit, n):
    for qubit in range(n//2):
        circuit.swap(qubit, n-qubit-1)
    return circuit

def qft(circuit, n):
    """QFT on the first n qubits in circuit"""
    qft_rotations(circuit, n)
    swap_registers(circuit, n)
    return circuit

def inverse_qft(circuit, n, s):
    """Does the inverse QFT on the first n qubits in circuit"""
    # First we create a QFT circuit of the correct size:
    qft_circ = qft(QuantumCircuit(n), n)
    # Then we take the inverse of this circuit
    invqft_circ = qft_circ.inverse()
    # And add it to the first n qubits in our existing circuit
    circuit.append(invqft_circ, circuit.qubits[s:n+s])
    return circuit.decompose() # .decompose() allows us to see the individual gates
    
def buildCircuit(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations

    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
                
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c)
            had3(h,clause[0],clause[1],clause[2],c,par)
            c.reset(h)
            
            c.barrier(qr)  
            
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitNOT(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations

    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
                
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c)
            not3(h,clause[0],clause[1],clause[2],c,par)
            c.reset(h)
            
            c.barrier(qr)  
            
    if measure:
        c.measure(qr[:-1],cr)
            
    return c


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