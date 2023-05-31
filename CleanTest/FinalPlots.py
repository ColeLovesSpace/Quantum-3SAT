from QSAT import *
from Shonings import *
from StatSATtest import *
import multiprocessing
from multiprocessing import Process, Pool

# This program will statistically test how classical stochastic algorithms compare with 
# the quantum 3-SAT algorithm. First we will study how they compare when run on problems with a single solution. 

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
            
    
def setup(n):
    par = {}
    par['nQ'] = n
    par['sat'] = buildSatSingleSolution(par['nQ'])
    par['nC'] = len(par['sat'])

    #Plotting Conditionals
    par['latex'] = True
    par['statevector'] = True
    par['saveEnd'] = False
    
    par['measure'] = False
    par['had'] = True
    
    par['gpu'] = False
    
    return par
# Test and compare quantum and classical algorithms


def singleSol(n,imax,s):
    
    par['figName'] = "plots/singleSolution/n{n:0>3}-s{s:0>6}.png".format(n=par['nQ'],s=par['numRepeats'])
    par['xAsis'] = "$\frac{c}{n}$"
    par['yAsis'] = "$P(\text{solved})_{i<"+"+{}".format(par['maxIt'])+"}$"
    
    BPPSATpi = []
    Shoningpi = []
    QSATpi = []
    
    BPPSATpiStep = []
    ShoningpiStep = []
    QSATpiStep = []
    
    cores = multiprocessing.cpu_count()
    pool = Pool(cores)
        
    

    plt.rcParams['text.usetex'] = True
    
    async_resultsC = []
    async_resultsS = []
    async_resultsQSAT = []
    
    
    
    start = time.time()
    index = 0
    par = setup(n)
    its = range(1,imax+1)
    
    for i in its:
        async_resultsC.append([])
        async_resultsS.append([])
        async_resultsStat.append([])
     
        #   Run statistical test of modified Shonings algorithm
        async_resultsC[index].append(pool.apply_async(classicalSolve,args=(n,SAT,nc,i,s)))
        #   Run statistical test of regular Shonings algorithm
        async_resultsS[index].append(pool.apply_async(classicalSolveS,args=(n,SAT,nc,i,s)))
        #   Run statistical test of regular Ising algorithm
        async_resultsStat[index].append(pool.apply_async(classicalSolveStat,args=(n,SAT,nc,i,s)))
        
      
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
    n = 6
#     maxC = 12*n
    i = 1000000#n**2
    numIt = 20
    divs = 10
    start = 3
    end = 8
    print(n)
    test4(n,ns,i,numIt,divs,start,end)
    
if __name__ == "__main__":
    main()

