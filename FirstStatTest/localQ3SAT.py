# Importing standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, Aer, IBMQ
from qiskit import *
# from qiskit.tools.jupyter import *
from qiskit.visualization import *
# from ibm_quantum_widgets import *
import qiskit.providers.aer
# import QasmSimulator
import matplotlib.pyplot as plt
import qiskit.quantum_info as qi

# Loading your IBM Quantum account(s)
# provider = IBMQ.load_account()

# Needed for functions
# import numpy as np
from numpy import pi
import time
# from copy import deepcopy
# import pandas as pd

# define functions for circuit
def OR(a, b, c, d, e, circuit):
    circuit.x(a)
    circuit.x(b)
    circuit.x(c)
    circuit.ccx(a, b, d)
    circuit.ccx(c, d, e)
    circuit.x(a)
    circuit.x(b)
    circuit.x(c)
    circuit.x(e)
    circuit.reset(d)


def AND(a, b, c, d, e, circuit):
    circuit.ccx(a, b, d)
    circuit.ccx(c, d, e)
    circuit.reset(d)


def had3(e, a, b, c, circuit, par):
    circuit.ch(e, abs(a) - 1)
    circuit.ch(e, abs(b) - 1)
    circuit.ch(e, abs(c) - 1)


def hadAll(e, a, b, c, circuit, par):
    for b in range(par['nQ']):
        circuit.ch(e, b)


#     circuit.ch(e,abs(a)-1)
#     circuit.ch(e,abs(b)-1)
#     circuit.ch(e,abs(c)-1)


def remStates(a, b, c, d, circuit):
    if a < 0:
        circuit.x(abs(a) - 1)
    if b < 0:
        circuit.x(abs(b) - 1)
    if c < 0:
        circuit.x(abs(c) - 1)

    circuit.rcccx(abs(a) - 1, abs(b) - 1, abs(c) - 1, d)

    if a < 0:
        circuit.x(abs(a) - 1)
    if b < 0:
        circuit.x(abs(b) - 1)
    if c < 0:
        circuit.x(abs(c) - 1)


def endRem(a, b, c, d, e, circuit):
    if a < 0:
        circuit.x(abs(a) - 1)
    if b < 0:
        circuit.x(abs(b) - 1)
    if c < 0:
        circuit.x(abs(c) - 1)


def qft3(c, q1, q2, q3):
    c.h(abs(q3) - 1)
    c.cp(pi / 4, abs(q1) - 1, abs(q3) - 1)
    c.cp(pi / 2, abs(q2) - 1, abs(q3) - 1)
    c.h(abs(q2) - 1)
    c.cp(pi / 2, abs(q1) - 1, abs(q2) - 1)
    c.h(abs(q1) - 1)
    c.swap(abs(q1) - 1, abs(q3) - 1)


def qft3i(c, q1, q2, q3):
    c.swap(abs(q1) - 1, abs(q3) - 1)
    c.h(abs(q1) - 1)
    c.cp(-pi / 2, abs(q1) - 1, abs(q2) - 1)
    c.h(abs(q2) - 1)
    c.cp(-pi / 2, abs(q2) - 1, abs(q3) - 1)
    c.cp(-pi / 4, abs(q1) - 1, abs(q3) - 1)
    c.h(abs(q3) - 1)


def qft_rotations(circuit, n):
    """Performs qft on the first n qubits in circuit (without swaps)"""
    if n == 0:
        return circuit
    n -= 1
    circuit.h(n)
    for qubit in range(n):
        circuit.cp(pi / 2 ** (n - qubit), qubit, n)
    # At the end of our function, we call the same function again on
    # the next qubits (we reduced n by one earlier in the function)
    qft_rotations(circuit, n)


def swap_registers(circuit, n):
    for qubit in range(n // 2):
        circuit.swap(qubit, n - qubit - 1)
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
    circuit.append(invqft_circ, circuit.qubits[s:n + s])
    return circuit.decompose()  # .decompose() allows us to see the individual gates


def buildSatSingleSolution(numVar):
    sat = [[1, 2, 3], [-1, 2, 3], [1, -2, 3], [1, 2, -3], [-1, -2, 3], [-1, 2, -3], [1, -2, -3]]
    for n in range(4, numVar + 1):
        sat.append([2, 3, n])
        sat.append([2, 3, -n])
        sat.append([-2, -3, n])

    return sat

# Build circuit
# SAT problem that will be implemented in circuit
# sat = [[-1,-2,-3],[-1,2,3],[1,-2,3],[1,2,-3],[-1,-2,3],[-1,2,-3],[1,-2,-3],[-2,-3,4],[-2,-3,-4],[-2,-3,5],[-2,-3,-5],[-2,-3,6],[-2,-3,-6],[-2,-3,7],[-2,-3,-7],[2,3,4],[2,3,5],[2,3,6],[2,3,7]]

def buildCircuit(par):
    helper = 1  # number of ancillary qubits
    h = par['nQ']  # qubit that can be used for calculations

    qr = QuantumRegister(par['nQ'] + helper, 'q')
    cr = ClassicalRegister(par['nQ'], 'c')
    c = QuantumCircuit(qr, cr)

    # Hadimard all qubits used in SAT
    for x in range(par['nQ']):
        c.h(x)

    for n in range(par['numIt']):  # Repeat circuit numIt times
        for clause in par['sat']:  # loop through all clauses in the SAT instance

            c.reset(h)
            remStates(clause[0], clause[1], clause[2], h, c)
            had3(h, clause[0], clause[1], clause[2], c, par)

            c.barrier(qr)

    return c


def setup(nQ):
    par = {}
    par['nQ'] = nQ

    # Plotting Conditionals
    par['latex'] = True
    par['statevector'] = True
    par['measure'] = False
    par['saveEnd'] = False
    par['figName'] = "Plots/TestSave2/{:0>6}.png"
    par['sat'] = buildSatSingleSolution(par['nQ'])  # [[1,2,3],[4,5,6]]
    par['numIt'] = nQ  # len(par['sat'])
    print(par['numIt'])
    state = ''
    for q in range(nQ):
        state = state + '0'
    par['state'] = state

    return par


par = setup(9)

# par['sat'] = [[-1,2,-3],[1,-3,-4],[2,3,4]]#
par['numIt'] = 1

# Build circuit
start = time.time()
c = buildCircuit(par)


if (par['statevector']):
    # Execute and get counts
    DM = qi.DensityMatrix.from_instruction(c)
    probs = DM.probabilities_dict()
    DM2 = qi.partial_trace(DM, [par['nQ']])
    probs = DM2.probabilities_dict()

print(start - time.time())

names = list(probs.keys())
values = list(probs.values())
print("this is it: ", list(probs.values())[0])
# plt.bar(probs)

plt.bar(range(len(probs)), values, tick_label=names)
plt.xticks(rotation=65)

# print(values)
# print(c)
ax = plt.gca()
ax.set_ylim([0, 1])
# name = "Plots/HardSAT20/{:0>3}.png".format(x)

name = "plots2/testPlot.png"
plt.savefig(name)

# dm.draw('latex', prefix='\\rho_{AB} = ')

b = circuit_drawer(c, output='mpl')
name = "plots2/testPlot2.png"
b.savefig(name)
