# Needed for functions
import time
import numpy as np
from numpy import pi
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import random
from math import gcd as bltin_gcd

#import packages for quantum algorithm
from qiskit import *
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import qiskit.quantum_info as qi
from qiskit.quantum_info.operators import Operator

from qiskit import Aer, transpile
# import qiskit.providers.aer
from qiskit_aer import AerError

def coprime(a, b):
    return bltin_gcd(a, b) == 1

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
    
def or3(d,a,b,c,circuit,offset):
    if a > 0:
        circuit.x(abs(a)-1+offset)
    if b > 0:
        circuit.x(abs(b)-1+offset)
    if c > 0:
        circuit.x(abs(c)-1+offset)
    
    circuit.rcccx(abs(a)-1+offset,abs(b)-1+offset,abs(c)-1+offset,d)
        
    if a > 0:
        circuit.x(abs(a)-1+offset)
    if b > 0:
        circuit.x(abs(b)-1+offset)
    if c > 0:
        circuit.x(abs(c)-1+offset)
    circuit.x(d)
    
def measureFlip(e,a,b,c,circuit,par):
    circuit.measure(e,0)
    
    circuit.h(abs(a)-1).c_if(0,1)
    circuit.h(abs(b)-1).c_if(0,1)
    circuit.h(abs(c)-1).c_if(0,1)
    
def hadAll(e,circuit,par):
    for b in range(par['nQ']):
        circuit.ch(e,b)
#     circuit.ch(e,abs(a)-1)
#     circuit.ch(e,abs(b)-1)
#     circuit.ch(e,abs(c)-1)

def notAll(e,circuit,par):
    for b in range(par['nQ']):
        circuit.cx(e,b)
#     circuit.ch(e,abs(a)-1)
#     circuit.ch(e,abs(b)-1)
#     circuit.ch(e,abs(c)-1)

def addCost(a,b,c,circuit,par,theta):
    loc = 0
    if(a>0):
        loc+=1
    if(b>0):
        loc+=2
    if(c>0):
        loc+=4
#     ee = np.eye(2 ** 3,dtype=complex)*np.exp(-1j*theta)
#     ee[loc][loc] = 1
    ee = np.eye(2 ** 3,dtype=complex)
    ee[loc][loc] = np.exp(-1j*theta)
    op = Operator(ee)
    circuit.append(op,[abs(a)-1,abs(b)-1,abs(c)-1])
    
def addPhase(circuit,par,lamb):
    for c in range(par['nQ']):
        circuit.rx(-lamb,c)
        
def remStates(a,b,c,d,circuit,offset):
    if a < 0:
        circuit.x(abs(a)-1+offset)
    if b < 0:
        circuit.x(abs(b)-1+offset)
    if c < 0:
        circuit.x(abs(c)-1+offset)
    
    circuit.rcccx(abs(a)-1+offset,abs(b)-1+offset,abs(c)-1+offset,d)
        
#     if a < 0:
#         circuit.x(abs(a)-1+offset)
#     if b < 0:
#         circuit.x(abs(b)-1+offset)
#     if c < 0:
#         circuit.x(abs(c)-1+offset)
         
def endRem(a,b,c,circuit):
    if a < 0:
        circuit.x(abs(a)-1)
    if b < 0:
        circuit.x(abs(b)-1)
    if c < 0:
        circuit.x(abs(c)-1)
        
def phase3(circuit,a,b,c,h):
    num = (a/abs(a))*2**0 + (b/abs(b))*2**1 + (c/abs(c))*2**2
    circuit.crz(num*pi,h,abs(a)-1)
    circuit.crz(num*pi/2,h,abs(b)-1)
    circuit.crz(num*pi/4,h,abs(c)-1)
#     Works pretty good
#     num = (max(0,a)/a)*2**2 + (max(0,b)/b)*2**1 + (max(0,c)/c)*2**0
#     circuit.crz(-num*pi/4,h,abs(a)-1)
#     circuit.crz(-num*pi/2,h,abs(b)-1)
#     circuit.crz(-num*pi,h,abs(c)-1)

def phase32(circuit,a,b,c,h):
    num = int(max(0,a/abs(a))*2**0 + max(0,b/abs(b))*2**1 + max(0,c/abs(c))*2**2)
    co = random.randint(1,7)
    while(not coprime(num,co)):
        co = random.randint(1,7)
    circuit.crz(co*pi,h,abs(a)-1)
    circuit.crz(co*pi/2,h,abs(b)-1)
    circuit.crz(co*pi/4,h,abs(c)-1)
#     Works pretty good
#     num = (max(0,a)/a)*2**2 + (max(0,b)/b)*2**1 + (max(0,c)/c)*2**0
#     circuit.crz(-num*pi/4,h,abs(a)-1)
#     circuit.crz(-num*pi/2,h,abs(b)-1)
#     circuit.crz(-num*pi,h,abs(c)-1)
    
def phaseAll(h,circuit,par,a,b,c):
#     for p in range(par['nQ']):
#         circuit.crz(pi,h,p)  
#     num = (a/((a**2)**(1/2)))*2**a
    num = 2**abs(a)
           
    for p in range(par['nQ']):
        circuit.crz(num*pi/(2**p) + (max(0,-a)/(-a))*pi,h,p)
                                                               
#     num = (b/((b**2)**(1/2)))*2**b
    num = 2**abs(b) 
           
    for p in range(par['nQ']):
        circuit.crz(num*pi/(2**p) + (max(0,-b)/(-b))*pi,h,p)
#         circuit.crz(num*pi/(2**(par['nQ']-p-1)),h,p)
                                                               
#     num = (c/((c**2)**(1/2)))*2**c
    num = 2**abs(c) 
    for p in range(par['nQ']):
        circuit.crz(num*pi/(2**p) + (max(0,-c)/(-c))*pi,h,p)                                                        
                                                            
def phase31(circuit,a,b,c,h):
    circuit.crz(pi,h,abs(a)-1)
    circuit.crz(pi/2,h,abs(b)-1)
    circuit.crz(pi/4,h,abs(c)-1)
        
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
    circuit.append(invqft_circ, circuit.qubits[:n])
#     circuit.qubits[s:n+s])
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
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
            had3(h,clause[0],clause[1],clause[2],c,par)
            endRem(clause[0],clause[1],clause[2],c)
            c.reset(h)
            
            c.barrier(qr)  
            
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitCost(par,measure):
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
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
#             remStates(clause[0],clause[1],clause[2],h,c,0)
#             had3(h,clause[0],clause[1],clause[2],c,par)
#             c.reset(h)
            
#             c.barrier(qr) 
            
        addPhase(c,par,par['beta'])
            
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
            c.barrier(qr) 
            
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitQFT(par,measure):
    helper = 1 #number of ancillary qubits
    h = 2*par['nQ'] #qubit that can be used for calculations
    
    qr = QuantumRegister(2*par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
    
    #copy state onto other qubits
    for q in range(par['nQ']):
        c.cnot(q,par['nQ']+q)
            
    qft(c,par['nQ'])
                
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,par['nQ'])
#             hadAll(h,c,par)
            had3(h,clause[0],clause[1],clause[2],c,par)
            c.reset(h)
            
            c.barrier(qr) 
            
#         addPhase(c,par,par['beta'])
            
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
#             addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
#             c.barrier(qr) 
    inverse_qft(c,par['nQ'],0)
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitQFTor(par,measure):
    helper = 3 #number of ancillary qubits
    h = 2*par['nQ'] #qubit that can be used for calculations
    h1 = 2*par['nQ']   #qubit that can be used for calculations
    h2 = 2*par['nQ']+2  #qubit where 
    flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(2*par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
    c.x(h1+flip)
    
    #copy state onto other qubits
    for q in range(par['nQ']):
        c.cnot(q,par['nQ']+q)
            
    qft(c,par['nQ'])
                
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h2,c,par['nQ'])
#             hadAll(h,c,par)
            
#             c.reset(h)
            
            
            
            c.x(h2)
            c.ccx(h1+flip,h2,h1)
            h1 = h1 + flip
            flip = -flip
            c.reset(h1)
            c.reset(h2)
            
            c.barrier(qr) 
            
        hadAll(h1+flip,c,par)
#         addPhase(c,par,par['beta'])
            
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
#             addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
#             c.barrier(qr) 
    inverse_qft(c,par['nQ'],0)
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitQFTorNOT(par,measure):
    helper = 3 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations
    h1 = par['nQ']   #qubit that can be used for calculations
    h2 = par['nQ']+2  #qubit where 
    flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
    c.x(h1+flip)
    
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance

            remStates(clause[0],clause[1],clause[2],h2,c,0)
#             hadAll(h,c,par)
            
#             c.reset(h)
            
            
            
            c.x(h2)
            c.ccx(h1+flip,h2,h1)
            h1 = h1 + flip
            flip = -flip
            c.reset(h1)
            c.reset(h2)
            
            c.barrier(qr) 
            
#         hadAll(h1+flip,c,par)
#         addPhase(c,par,par['beta'])
            
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
#             addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
#             c.barrier(qr) 
    c.x(h1+flip)
    qft(c,par['nQ'])
    notAll(h1+flip,c,par)
    inverse_qft(c,par['nQ'],0)
    if measure:
        c.measure(qr[:-helper],cr)
            
    return c

def buildCircuitQFTNOT(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations
#     h1 = par['nQ']   #qubit that can be used for calculations
#     h2 = par['nQ']+2  #qubit where 
#     flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
#     c.x(h1+flip)
    
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
#         cnum = 0
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
#             hadAll(h,c,par)
            
#             c.reset(h)
            qft(c,par['nQ'])
            notAll(h,c,par)
#             hadAll(h,c,par)
#             notAll(h,c,par)
            inverse_qft(c,par['nQ'],0)
            endRem(clause[0],clause[1],clause[2],c)
            
            c.reset(h)
#             cnum+=1
            c.barrier(qr) 
            
            
    
    if measure:
        c.measure(qr[:-helper],cr)
            
    return c
                                                               
def buildCircuitQFTphase(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations
#     h1 = par['nQ']   #qubit that can be used for calculations
#     h2 = par['nQ']+2  #qubit where 
#     flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
#     c.x(h1+flip)
    
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
#         cnum = 0
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
#             hadAll(h,c,par)
            
#             c.reset(h)
            qft(c,par['nQ'])
            phaseAll(h,c,par,clause[0],clause[1],clause[2])
#             hadAll(h,c,par)
#             notAll(h,c,par)
            inverse_qft(c,par['nQ'],0)
            endRem(clause[0],clause[1],clause[2],c)
            c.reset(h)
#             cnum+=1
            c.barrier(qr) 
            
            
    
    if measure:
        c.measure(qr[:-helper],cr)
            
    return c

def buildCircuitQFT3phase(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations
#     h1 = par['nQ']   #qubit that can be used for calculations
#     h2 = par['nQ']+2  #qubit where 
#     flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
#     c.x(h1+flip)
    
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
#             hadAll(h,c,par)
            
#             c.reset(h)
            qft3(c,clause[0],clause[1],clause[2])
            phase32(c,clause[0],clause[1],clause[2],h)
            qft3i(c,clause[0],clause[1],clause[2])
            endRem(clause[0],clause[1],clause[2],c)
            c.reset(h)
            
            c.barrier(qr) 
            
#         hadAll(h1+flip,c,par)
#         addPhase(c,par,par['beta'])
            
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
#             addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
#             c.barrier(qr) 
    
    if measure:
        c.measure(qr[:-helper],cr)
            
    return c

def buildCircuitQFT3NOT(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations
#     h1 = par['nQ']   #qubit that can be used for calculations
#     h2 = par['nQ']+2  #qubit where 
#     flip = 1   #variable to change position of h1
    
    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
            
    # flip helper variable so first gate is added to rolling and
#     c.x(h1+flip)
    
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
#             hadAll(h,c,par)
            
#             c.reset(h)
            qft3(c,clause[0],clause[1],clause[2])
            not3(h,clause[0],clause[1],clause[2],c,par)
#             had3(h,clause[0],clause[1],clause[2],c,par)
#             not3(h,clause[0],clause[1],clause[2],c,par)
#             had3(h,clause[0],clause[1],clause[2],c,par)
#             not3(h,clause[0],clause[1],clause[2],c,par)
            qft3i(c,clause[0],clause[1],clause[2])
            endRem(clause[0],clause[1],clause[2],c)
            c.reset(h)
            
            c.barrier(qr) 
            
#         hadAll(h1+flip,c,par)
#         addPhase(c,par,par['beta'])
            
#         for clause in par['sat']: # loop through all clauses in the SAT instance
            
#             addCost(clause[0],clause[1],clause[2],c,par,par['gamma'])            
#             c.barrier(qr) 
    
    if measure:
        c.measure(qr[:-helper],cr)
            
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
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
            not3(h,clause[0],clause[1],clause[2],c,par)
            c.reset(h)
            
            c.barrier(qr)  
            
    if measure:
        c.measure(qr[:-1],cr)
            
    return c

def buildCircuitMeasure(par,measure):
    helper = 1 #number of ancillary qubits
    h = par['nQ'] #qubit that can be used for calculations

    qr = QuantumRegister(par['nQ']+helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
#     cr2 = ClassicalRegister(1, 'c2')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    if par['had']:
        for x in range(par['nQ']):
            c.h(x)
                
    for n in range(par['numIt']): # Repeat circuit numIt times
#         random.shuffle(par['sat'])
        for clause in par['sat']: # loop through all clauses in the SAT instance
            
            
            remStates(clause[0],clause[1],clause[2],h,c,0)
            measureFlip(h,clause[0],clause[1],clause[2],c,par)
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