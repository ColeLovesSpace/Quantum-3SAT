from QSAT import *
from Shonings import *
from StatSATtest import *
import multiprocessing
from multiprocessing import Process, Pool

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


def test3(n,ns,i,numIt,divs,begin,end):
    
    BPPSATpi = []
    Shoningpi = []
    Statpi = []
    
#     (maxC/n)
    plt.rcParams['text.usetex'] = True
    index = 0
    for ncc in np.linspace(begin*n,end*n,divs):
        index+=1
        nc = int(ncc)
        par = setup(n,nc,ns)
        
#         start = time.time()
        cores = multiprocessing.cpu_count()
        pool = Pool(cores)
        #   Run statistical test of modified Shonings algorithm
        start = time.time()
        BPPSATavg = 0
        results = []
        async_results = [pool.apply_async(classicalSolve,args=(n,SAT,nc,i,numIt)) for SAT in par['satList']]
        results = [ar.get() for ar in async_results]
        for res in results:
            cutoff = [x for x in res if x < i]
            values, base = np.histogram(cutoff, bins=n**2)
            #evaluate the cumulative
            cumulative = np.cumsum(values)
            BPPSATavg += cumulative[-1]/numIt
        BPPSATpi.append(BPPSATavg/ns)
        print("BPP:",BPPSATpi)
        
        print("BPP time to solve: " , start - time.time())
            
        #   Run statistical test of regular Shonings algorithm
        start = time.time()
        Shoningavg = 0
        results = []
        async_results = [pool.apply_async(classicalSolveS,args=(n,SAT,nc,i,numIt)) for SAT in par['satList']]
        results = [ar.get() for ar in async_results]
        for resS in results:
            cutoffS = [x for x in resS if x < i]
            valuesS, baseS = np.histogram(cutoffS, bins=n**2)
            #evaluate the cumulative
            cumulativeS = np.cumsum(valuesS)
            Shoningavg += cumulativeS[-1]/numIt
        Shoningpi.append(Shoningavg/ns)
        print("SHO:",Shoningpi)
        
        print("Shoning time to solve: " , start - time.time())
            
        #   Run statistical test of regular Ising algorithm
        start = time.time()
        Statavg = 0
        results = []
        async_results = [pool.apply_async(classicalSolveStat,args=(n,SAT,nc,i,numIt)) for SAT in par['satList']]
        results = [ar.get() for ar in async_results]
        for res in results:
            cutoffStat = [x for x in resS if x < i]
            valuesStat, baseStat = np.histogram(cutoffStat, bins=n**2)
            #evaluate the cumulative
            cumulativeStat = np.cumsum(valuesStat)
            Statavg += cumulativeStat[-1]/numIt
        Statpi.append(Statavg/ns)
        print("STAT:",Statpi)
        
        print("Stat time to solve: " , start - time.time())

        print("Clauses " + str(nc) + " took " + str(start-time.time()))
#         print("start",start)
#         print("index",index)
#         print()

#         xAxis = list(np.linspace(start,nc/n,index,endpoint=True))
        xAxis = np.linspace(begin,nc/n,index,endpoint=True)
                     
        name = "n{:0>3}".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)
#         xlabel = "$\langle P(\text{solved}|"+"i<{})$"
        fig = plt.figure(figsize=(9, 5), facecolor='lightskyblue')
        fig.suptitle(name)
        ax = fig.add_subplot(111)
        
        ax.set_title(r'$\langle P($solved$|i) \rangle_{\Phi}$', loc='left', fontstyle='oblique', fontsize='medium')
        
        ax.plot(xAxis, BPPSATpi, c='blue',label="ModShonings")
        ax.plot(xAxis, Shoningpi, c='green',label="Shonings")
        ax.plot(xAxis, Statpi, c='red',label="IsingSAT")
        ax.legend()
        ax.set_xlabel(r'$\frac{c}{n}$')
        name = "plots/clean/n{:0>3}".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)+".png"
        fig.savefig(name)
        fig.clf()
        
def try250():
#   Initialize parameters for SAT instance
    
#     maxC = 12*n
    i = 100000#int(n**3 * 4/100)
    numIt = 5
    ns = 1

    BPPSATpi = []
    Shoningpi = []
    Statpi = []    
    
#     satList = []
#     fo
#     SAT ,n ,nc = getSATS(y+1)
    
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
        res = classicalSolve(n,SAT,nc,i,numIt)
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
    
    print("Number of cpu : ", multiprocessing.cpu_count())
#     Initialize parameters for SAT instance
    n = 10
    ns = 80
#     maxC = 12*n
    i = 100000#n**2
    numIt = 100
    divs = 20
    start = 3
    end = 8
    test3(n,ns,i,numIt,divs,start,end)
#     try250()
    
if __name__ == "__main__":
    main()

