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
            
    
def setup(n,nc,nSAT):
    par = {}
    par['nQ'] = n
    par['nC'] = nc
    par['nSAT'] = nSAT
    par['sat'] = buildSatSingleSolution(par['nQ'])
    par['satList'] = SATset(n,nc,nSAT)
    par['numRepeats'] = 1000
    par['numIt'] = 5
    par['maxIt'] = n**2

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
        name = "plots/StatMult/n{:0>3}b".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)+".png"
        fig.savefig(name)
        fig.clf()
        
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

def test4(n,ns,i,numIt,divs,begin,end):
    
    BPPSATpi = []
    Shoningpi = []
    Statpi = []
    
    BPPSATpiStep = []
    ShoningpiStep = []
    StatpiStep = []
    
    cores = multiprocessing.cpu_count()
    pool = Pool(cores)
        
    
#     (maxC/n)
    plt.rcParams['text.usetex'] = True
    index = 0
    start = time.time()
    async_resultsC = []
    async_resultsS = []
    async_resultsStat = []
    for ncc in np.linspace(begin*n,end*n,divs):
        
        nc = int(ncc)
        print(nc)
        par = setup(n,nc,ns)
        
        async_resultsC.append([])
        async_resultsS.append([])
        async_resultsStat.append([])
     
        for SAT in par['satList']:
            #   Run statistical test of modified Shonings algorithm
            async_resultsC[index].append(pool.apply_async(classicalSolve,args=(n,SAT,nc,i,numIt)))
            #   Run statistical test of regular Shonings algorithm
            async_resultsS[index].append(pool.apply_async(classicalSolveS,args=(n,SAT,nc,i,numIt)))
            #   Run statistical test of regular Ising algorithm
            async_resultsStat[index].append(pool.apply_async(classicalSolveStat,args=(n,SAT,nc,i,numIt)))
        
      
        index+=1
    
    print("TimeToQueue: " , start - time.time())
    index = 0
    
    
#     start = time.time()
    for ncc in np.linspace(begin*n,end*n,divs):
        
        BPPSATavg = 0
        BPPSATavgStep = 0
        results = []
        nc = int(ncc)
        start = time.time()
        results = [ar.get() for ar in async_resultsC[index]]
        for res in results:
            cutoff = [x for x in res if x < i]
            values, base = np.histogram(cutoff, bins=n**2)
            #evaluate the cumulative
            tots = 0
            for x in res:
                tots += x/i
            if len(res)>0:
                BPPSATavgStep += tots/len(res)
            cumulative = np.cumsum(values)
            BPPSATavg += cumulative[-1]/numIt
        BPPSATpi.append(BPPSATavg/ns)
        BPPSATpiStep.append(BPPSATavgStep/ns)
        print("BPP:",BPPSATpi)
        
#         print("BPP time to solve: " , start - time.time())
            
        
        start = time.time()
        Shoningavg = 0
        ShoningavgStep = 0
        results = []
        results = [ar.get() for ar in async_resultsS[index]]
        for res in results:
            cutoffS = [x for x in res if x < i]
            valuesS, baseS = np.histogram(cutoffS, bins=n**2)
            #evaluate the cumulative
            tots = 0
            for x in res:
                tots += x/i
            if len(res)>0:
                ShoningavgStep += tots/len(res)
            cumulativeS = np.cumsum(valuesS)
            Shoningavg += cumulativeS[-1]/numIt
        Shoningpi.append(Shoningavg/ns)
        ShoningpiStep.append(ShoningavgStep/ns)
        print("SHO:",Shoningpi)
        
#         print("Shoning time to solve: " , start - time.time())
            
        
        start = time.time()
        Statavg = 0
        StatavgStep = 0
        results = []
        
        results = [ar.get() for ar in async_resultsStat[index]]
        for res in results:
            cutoffStat = [x for x in res if x < i]
            valuesStat, baseStat = np.histogram(cutoffStat, bins=n**2)
            #evaluate the cumulative
            tots = 0
            for x in res:
                tots += x/i
            if len(res)>0:
                StatavgStep += tots/len(res)
            cumulativeStat = np.cumsum(valuesStat)
            Statavg += cumulativeStat[-1]/numIt
        Statpi.append(Statavg/ns)
        StatpiStep.append(StatavgStep/ns)
        print("STAT:",Statpi)
        
#         print("Stat time to solve: " , start - time.time())
        index+=1
    
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
        name = "plots/AllCore/n{:0>3}b".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)+".png"
        fig.savefig(name)
        fig.clf()
        
        xAxis = np.linspace(begin,nc/n,index,endpoint=True)
                     
        name = "Steps n{:0>3}".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)
#         xlabel = "$\langle P(\text{solved}|"+"i<{})$"
        fig = plt.figure(figsize=(9, 5), facecolor='lightskyblue')
        fig.suptitle(name)
        ax = fig.add_subplot(111)
        
        ax.set_title(r'$\langle i \rangle_{\Phi}$', loc='left', fontstyle='oblique', fontsize='medium')
        
        ax.plot(xAxis, BPPSATpiStep, c='blue',label="ModShonings")
        ax.plot(xAxis, ShoningpiStep, c='green',label="Shonings")
        ax.plot(xAxis, StatpiStep, c='red',label="IsingSAT")
        ax.legend()
        ax.set_xlabel(r'$\frac{c}{n}$')
        name = "plots/AllCore/S-n{:0>3}b".format(n)+"i{}".format(i)+"nI{}".format(numIt)+"ns{}".format(ns)+".png"
        fig.savefig(name)
        fig.clf()

def main():
    
    print("Number of cpu : ", multiprocessing.cpu_count())
#     Initialize parameters for SAT instance
    n = 500
    ns = 20
#     maxC = 12*n
    i = 1000000#n**2
    numIt = 20
    divs = 10
    start = 3
    end = 8
    print(n)
    test4(n,ns,i,numIt,divs,start,end)
    
    
def main2():
    try250()
    
if __name__ == "__main__":
    main()

