from QSAT import *
from Shonings import *
from StatSATtest import *

# Exhaustive search for solution

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
            
# Test and compare quantum and classical algorithms


def test3():
    #     # Initialize a GPU backend
#     # Note that the cloud instance for tutorials does not have a GPU
#     # so this will raise an exception.
#     try:
#         simulator_density_matrix = Aer.get_backend('aer_simulator_density_matrix')
#         simulator_density_matrix.set_options(device='GPU')
#     except AerError as e:
#         print(e)
        
#   Initialize parameters for SAT instance
    n = 10
    maxC = 12*n
    i = 10#n**2
    numIt = 100000
    ns = 2
    
    QSATpi = []
    QSAT2pi = []
    BPPSATpi = []
    Shoningpi = []
    Statpi = []
    
    for ncc in range(1,int(maxC/10)):
        nc = ncc*10
#         print(nc)
        par = setup(n,nc,ns)
        
        QSATavg = 0
        QSAT2avg = 0
        BPPSATavg = 0
        Shoningavg = 0
        Statavg = 0
        start = time.time()
        for SAT in par['satList']:
        #   Run statistical test of modified Shonings algorithm
            res = classicalSolve2(n,SAT,nc,i,numIt)
            cutoff = [x for x in res if x < i]
            values, base = np.histogram(cutoff, bins=n**2)
        #     values, base = np.histogram(res, bins=n**2)
            #evaluate the cumulative
            cumulative = np.cumsum(values)
            BPPSATavg += cumulative[-1]/numIt
            BPPSATpi.append(BPPSATavg/ns)
            print("BPP:",BPPSATpi)
            
        #   Run statistical test of regular Shonings algorithm
            resS = classicalSolveS(n,SAT,nc,i,numIt)
            cutoffS = [x for x in resS if x < i]
            valuesS, baseS = np.histogram(cutoffS, bins=n**2)
        #     valuesS, baseS = np.histogram(resS, bins=n**2)
            #evaluate the cumulative
            cumulativeS = np.cumsum(valuesS)
            Shoningavg += cumulativeS[-1]/numIt
            Shoningpi.append(Shoningavg/ns)
            print("SHO:",Shoningpi)
            
        #   Run statistical test of regular Ising algorithm
            resStat = classicalSolveStat(n,SAT,nc,i,numIt)
            cutoffStat = [x for x in resS if x < i]
            valuesStat, baseStat = np.histogram(cutoffStat, bins=n**2)
        #     valuesS, baseS = np.histogram(resS, bins=n**2)
            #evaluate the cumulative
            cumulativeStat = np.cumsum(valuesStat)
            Statavg += cumulativeStat[-1]/numIt
            Statpi.append(Statavg/ns)
            print("STAT:",Statpi)

#             par['numIt'] = i
#             shots = 10000
#             par['sat'] = SAT

#             # Build circuit
#             qc = buildCircuit(par,True)
#             circ = transpile(qc, simulator_density_matrix, optimization_level=0)
#             job_density_matrix = simulator_density_matrix.run(circ, shots=shots)
#             counts_density_matrix = job_density_matrix.result().get_counts(0)

#     #         start = time.time()
#     #         print("Iteration " + str(it) + " took " + str(start-time.time()))

#             total = 0
#             for state in counts_density_matrix:
#                 if(checkSAT(par['sat'], int(state, 2), len(par['sat']), par['nQ'])):
#                     total += counts_density_matrix[state]

#             QSATavg += total/shots
            
#             # Build circuit2
#             qc2 = buildCircuitNOT(par,True)#*Borat voice*
#             circ2 = transpile(qc2, simulator_density_matrix, optimization_level=0)
#             job_density_matrix2 = simulator_density_matrix.run(circ2, shots=shots)
#             counts_density_matrix2 = job_density_matrix2.result().get_counts(0)

#     #         start = time.time()
#     #         print("Iteration " + str(it) + " took " + str(start-time.time()))

#             total = 0
#             for state in counts_density_matrix2:
#                 if(checkSAT(par['sat'], int(state, 2), len(par['sat']), par['nQ'])):
#                     total += counts_density_matrix2[state]

#             QSAT2avg += total/shots

            
            
    #         QSATpi.append(QSATavg/ns)
    #         QSAT2pi.append(QSAT2avg/ns)
            print("Clauses " + str(nc) + " took " + str(start-time.time()))

            xAxis = list(range(1,ncc+1))

            plt.plot(xAxis, BPPSATpi, c='blue')
            plt.plot(xAxis, Shoningpi, c='green')
            plt.plot(xAxis, Statpi, c='red')
    #         plt.plot(xAxis, QSATpi, c='red')
    #         plt.plot(xAxis, QSAT2pi, c='magenta')

            name = "plots/firstTest/2Alg3SATtest20{:0>3}.png".format(n)
            plt.savefig(name)
            plt.clf()
        
def try250():
#   Initialize parameters for SAT instance
    
#     maxC = 12*n
    i = 1000000#int(n**3 * 4/100)
    numIt = 5
    ns = 1

    BPPSATpi = []
    Shoningpi = []
    Statpi = []    
    
    for y in range(10):
        BPPSATavg = 0
        Shoningavg = 0
        Statavg = 0
        start = time.time()
        
        SAT ,n ,nc = getSATS(y+1)
    
        #   Run statistical test of regular Ising algorithm
        resStat = classicalSolveStat(n,SAT,nc,i,numIt)
        cutoffStat = [x for x in resStat if x < i]
        valuesStat, baseStat = np.histogram(cutoffStat, bins=n**2)
    #     valuesS, baseS = np.histogram(resS, bins=n**2)
        #evaluate the cumulative
        cumulativeStat = np.cumsum(valuesStat)
        Statavg += cumulativeStat[-1]/numIt

    #   Run statistical test of modified Shonings algorithm
        res = classicalSolve2(n,SAT,nc,i,numIt)
        cutoff = [x for x in res if x < i]
        values, base = np.histogram(cutoff, bins=n**2)
    #     values, base = np.histogram(res, bins=n**2)
        #evaluate the cumulative
        cumulative = np.cumsum(values)
        BPPSATavg += cumulative[-1]/numIt

    #   Run statistical test of regular Shonings algorithm
        resS = classicalSolveS(n,SAT,nc,i,numIt)
        cutoffS = [x for x in resS if x < i]
        valuesS, baseS = np.histogram(cutoffS, bins=n**2)
    #     valuesS, baseS = np.histogram(resS, bins=n**2)
        #evaluate the cumulative
        cumulativeS = np.cumsum(valuesS)
        Shoningavg += cumulativeS[-1]/numIt

        BPPSATpi.append(BPPSATavg/ns)
        Shoningpi.append(Shoningavg/ns)
        Statpi.append(Statavg/ns)
    #         QSATpi.append(QSATavg/ns)
    #         QSAT2pi.append(QSAT2avg/ns)
        print("Clauses " + str(nc) + " took " + str(start-time.time()))

#         xAxis = list(range(1,y))

#         plt.plot(xAxis, BPPSATpi, c='blue')
#         plt.plot(xAxis, Shoningpi, c='green')
#         plt.plot(xAxis, Statpi, c='red')
    #         plt.plot(xAxis, QSATpi, c='red')
    #         plt.plot(xAxis, QSAT2pi, c='magenta')
        print(y)
        print("BPP:",BPPSATpi)
        print("SHO:",Shoningpi)
        print("STAT:",Statpi)

#         name = "plots/firstTest/First3Alg3SATtest75{:0>3}.png".format(n)
#         plt.savefig(name)
#         plt.clf()


def main():
    test3()
    
if __name__ == "__main__":
    main()

