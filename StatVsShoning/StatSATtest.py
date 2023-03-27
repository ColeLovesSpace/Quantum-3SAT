import random
import math
import matplotlib.pyplot as plt
import numpy as np
import time

class Spin(object):
    def __init__(self,s):
        self.spin = s
        self.up = []
        self.down = []
        
    # connnect the spin to a clause. If orientation = 1, then add to up, else add to down connection list
    def connect(self,clause, orientation):
        if(orientation):
            self.up.append(clause)
        else:
            self.down.append(clause)
            
        if(orientation != self.spin):
            clause.red+=1
            
        clause.edges.append(self)
        clause.connections+=1
            
    def deltaH(self):
        rto = 0 # number of clauses flipped from red to orange
        otr = 0 # number of clauses flipped from orange to red
        #red means every connection is unsatisfied.
        #orange means one connection is satisfied.
        #green means more than one connection is satisfied
        if self.spin:
            for u in self.up:
                if u.red == u.connections-1:
                    otr += 1
            for d in self.down:
                if d.red == d.connections:
                    rto += 1
        else:
            for u in self.up:
                if u.red == u.connections:
                    rto += 1
            for d in self.down:
                if d.red == d.connections-1:
                    otr += 1
        
        return otr - rto
    
    def flip(self):
        self.spin = 1 - self.spin
        
        if self.spin:
            for u in self.up:
                u.red -= 1
            for d in self.down:
                d.red += 1
        else:
            for u in self.up:
                u.red += 1
            for d in self.down:
                d.red -= 1
        
class Clause:
    def __init__(self):
        self.red = 0
        self.connections = 0
        self.edges = []
        
        
class Graph:
    def __init__(self,n,SAT):
        self.SAT = SAT
        self.spins = []
        self.clauses = []
        self.c = len(SAT)
        self.H = 0
        for x in range(n):
            self.spins.append(Spin(random.choice([0,1])))
        
        for c in SAT:
            clause = Clause()
            for s in c:
                if s > 0:
                    self.spins[abs(s)-1].connect(clause,1)
                else:
                    self.spins[abs(s)-1].connect(clause,0)
                    
            self.clauses.append(clause)
            
            if clause.red == 3:
                self.H += 1
                
    def checkSAT(self):
        satisfied = True
        for c in self.SAT:
            numSat = 0
            for s in c:
                if s > 0:
                    if self.spins[abs(s)-1].spin == 1:
                        numSat += 1
                else:
                    if self.spins[abs(s)-1].spin == 0:
                        numSat += 1
            if numSat == 0:
                satisfied = False
            
        return satisfied
            
    def printState(self):
        state = ""
        for s in self.spins:
            state += str(s.spin)
        print(state)
                
    def solveSAT(self, iterations, Tb):
        i = 0
        solved = False
        flips = 0
        h = []
        h.append(self.H)
        kick = 0
        pdH = 0
        while i < iterations:
#             for s in self.spins:
            s = random.choice(self.spins)
            dH = s.deltaH()
            T = Tb + kick#*0.4#/(iterations/100000)#self.H/self.c
#             T = Tb + np.exp(kick)-1#/(iterations/100000)#self.H/self.c
#             T = 0.99
                
            pFlip = 1/(1+np.exp(dH/(T)))
            pdH = self.H
            
            if pFlip > random.random():
                s.flip()
                self.H += dH
                flips +=1
                h.append(self.H)
                i = i + 1
                
                if i%500000 == 0:
                    print("H",self.H)
                    print("T",T)
                    print("dH",dH)

            if self.H >= pdH:
                kick += 0.00001
            else:
                kick = 0
 
            if self.H == 0:
#                 self.printState()
#                 print("Number of iterations:",i)
                i = iterations
#                 print(self.checkSAT())
                break
            
        if self.H != 0: 
            print("Failed to find solution in " ,i, " steps")
        
#         print("Number of flips",flips)
#         print("Hamiltonian:",self.H)
#         plt.plot(h)
        return flips

    def solveSATclause(self, iterations, Tb):
        i = 0
        solved = False
        flips = 0
        h = []
        h.append(self.H)
        kick = 0.000001
        pdH = 0
        T = Tb
        while i < iterations:
#             for s in self.spins:
            i = i + 1
    
            T = (8*self.H/self.c)**3 + Tb + kick #*0.4#/(iterations/100000)# ((iterations-i)/iterations)
            pdH = self.H
#             for c in self.clauses:
            c = random.choice(self.clauses)
            if c.red == c.connections:
                pdH = self.H
                for s in c.edges:

                    dH = s.deltaH()
                    pFlip = 1/(1+np.exp(dH/(T)))
                    if pFlip > random.random():
                        s.flip()
                        self.H += dH
                        flips +=1
                        h.append(self.H)
                        
            if self.H >= pdH:
                kick += 0.00003
            else:
                kick = 0.000001
            T = (8*self.H/self.c)**3 + Tb + kick
            
            s = random.choice(self.spins)            
            dH = s.deltaH()
            pFlip = 1/(1+np.exp(dH/(T)))
            if pFlip > random.random():
                s.flip()
                self.H += dH
                flips +=1
                h.append(self.H)
                        
            if self.H >= pdH:
                kick += 0.00003
            else:
                kick = 0.000001

            if i%10000000 == 0:
                print("H",self.H)
                print("T",T)
#                 print("dH",dH)
                
            if self.H == 0:
#                 self.printState()
#                 print("Number of iterations:",i)
                i = iterations
#                 print(self.checkSAT())
                break
            
            
        if self.H != 0: 
            print("Failed to find solution in " ,i, " steps")
        
#         print("Number of flips",flips)
#         print("Hamiltonian:",self.H)
#         plt.plot(h)
        return np.array(h), flips

    def solveSATn(self, iterations):
        i = 0
        solved = False
        flips = 0
        h = []
        h.append(self.H)
        kick = 0
        pdH = 0
        while i < iterations:
#             for s in self.spins:
            toFlip = random.sample(self.spins,int(np.ceil(self.H*random.random()/2)))
            dH = 0
            for s in toFlip:
                dH += s.deltaH()
            
            if pdH == dH:
#                 kick += 0.01
                kick += 0.1
            else:
                pdH = dH
                kick = 0
                
            T = 0.2 + kick#*0.4#/(iterations/100000)#self.H/self.c
#             T = 0.2 + np.exp(kick)-1#/(iterations/100000)#self.H/self.c
#             T = 0.99
           
            pFlip = 1/(1+np.exp(dH/(T)))
        
            if pFlip > random.random():
                for s in toFlip:
                    s.flip()
                self.H += dH
                flips +=len(toFlip)
                h.append(self.H)
                i = i + 1
                
                if i%500000 == 0:
                    print("H",self.H)
                    print("T",T)
                    print("dH",dH)


            if self.H == 0:
#                 self.printState()
                print("Number of iterations:",i)
                i = iterations
                print(self.checkSAT())
                break
            
            
        if self.H != 0: 
            print("Failed to find solution in " ,i, " steps")
        
        print("Number of flips",flips)
#         print("Hamiltonian:",self.H)
#         plt.plot(h)
        return flips
        
    def solveSATclauseOnly(self, iterations, Tb):
        i = 0
        solved = False
        flips = 0
        h = []
        h.append(self.H)
        kick = 0.000001
        pdH = 0
        T = Tb
        itf = 0
        while i < iterations:
#             for s in self.spins:
            i = i + 1
            itf = i
        
            T = (8*self.H/self.c)**3 + Tb + kick #*0.4#/(iterations/100000)# ((iterations-i)/iterations)
            pdH = self.H
            for c in self.clauses:
#             c = random.choice(self.clauses)
                if c.red == c.connections:
                    pdH = self.H
                    for s in c.edges:

                        dH = s.deltaH()
                        pFlip = 1/(1+np.exp(dH/(T)))
                        if pFlip > random.random():
                            s.flip()
                            self.H += dH
                            flips +=1
                            h.append(self.H)

                if self.H >= pdH:
                    kick += 0.00003
                else:
                    kick = 0.000001
                T = (8*self.H/self.c)**3 + Tb + kick
            
                

#             if i%10000000 == 0:
#                 print("H",self.H)
#                 print("T",T)
# #                 print("dH",dH)
                
            if self.H == 0:
#                 self.printState()
#                 print("Number of iterations:",i)
                itf = i
                i = iterations
#                 print(self.checkSAT())
                break
            
            
#         if self.H != 0: 
#             print("Failed to find solution in " ,i, " steps")
        
#         print("Number of flips",flips)
#         print("Hamiltonian:",self.H)
#         plt.plot(h)
#         for s in self.spins:
#             s.spin = random.choice([0,1])

        return np.array(h), flips, itf
        
def buildSatSingleSolution(numVar):
    sat = [[-1,-2,-3],[-1,2,3],[1,-2,3],[1,2,-3],[-1,-2,3],[-1,2,-3],[1,-2,-3]]
    for n in range(4,numVar+1):
        sat.append([2,3,n])
        sat.append([-2,-3,-n])
        sat.append([-2,-3,n])
        
    return sat

def getSATS(s):
   ## 21 impossible
    filename = "../BPPSAT/SATS/UF250.1065.100/uf250-0" + str(s) + ".cnf"#"./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" ./SATS/UF75.325.100/uf75-024.cnf
#     filename = "../BPPSAT/SATS/UF75.325.100/uf75-0" + str(s + 75) + ".cnf"#"./SATS/uf20-91(1)/uf20-0"+str(s+1)+".cnf" ./SATS/UF75.325.100/uf75-024.cnf

    print(filename)

    with open(filename) as f:
        lines = f.readlines()

    dets = lines[7].strip().split(" ")
    n = int(dets[2])
    c = int(dets[4])
    lines = lines[8:-3]

    

    SAT = []
    for line in lines:
        l = np.array(line.strip().replace(" 0","").split(" ")).astype(int)
#         print(l)
        SAT.append(l)

#     SAT = np.array(SAT).astype('double')
#     for x in range(len(SAT)):
#         SAT[x] = -SAT[x]
    print(n, " ", len(SAT))
        
    return SAT, n, c

def getSATSrsa():
    filename = "./BPPSAT/SATS/rsa100.cnf"

    with open(filename) as f:
        lines = f.readlines()

    dets = lines[11].strip().split(" ")
    n = int(dets[2])
    c = int(dets[3])
    print(n)
    print(c)
    lines = lines[12:]

#     print(n, " ", c)

    SAT = []
    for line in lines:
        l = np.array(line.strip().replace(" 0","").split(" ")).astype(int)
#         print(l)
        SAT.append(l)

#     SAT = np.array(SAT).astype('double')
#     for x in range(len(SAT)):
#         SAT[x] = -SAT[x]
        
    return SAT, n, c

def test250():
    # n = 30
    # SAT = buildSatSingleSolution(n)
#     /Users/colecoughlin/Desktop/QSAT/SATS/UF75.325.100/uf75-04.cnf
    times = []
    nflips = []
    tot = time.time()
    for s in range(3):
        t = time.time()
        SAT, n, c = getSATS(s+1)
        # g = Graph(n,[[-1,2,3],[1,-2,4],[1,3,-4]])
        g = Graph(n,SAT)
    #     print("Hamiltonian: ", g.H)
    #     g.printState()
        num = g.solveSATclause(100000000,0.2)
        times.append(time.time()-t)
        nflips.append(num)
    #     g.printState()

        avg = 0
        for it in nflips:
            avg+=it
        avg/=len(nflips)
        print("average flips: ", avg)
        print("total time: ", time.time() - tot)
        plt.plot(times,c='red')
        plt.plot(nflips,c='blue')

        name = "plots2/StatTest/factorio{:0>3}.png".format(5)
        plt.savefig(name)
        
def testRSA():
    tot = time.time()
    
        
    SAT, n, c = getSATSrsa()
    
    g = Graph(n,SAT)
    #     print("Hamiltonian: ", g.H)
    #     g.printState()
#     Tb = 0.2
    Tb = 0.15
    h , flips = g.solveSAT(100000000000,Tb)
    
    g.printState()
    
    plt.plot(h,c='blue')
    name = "plots2/StatTest/rsa{:0>3}.png".format(5)
    plt.savefig(name)
    print(flips)
    print("total time: ", time.time() - tot)
    
    
def classicalSolveStat(n,SAT,nc,i,numIt):
    
    Tb = 0.1
    
    start = time.time()
    stepList = []
    
#     print("Stat time to build: " , start - time.time())
    start = time.time()
    for x in range(numIt):
        g = Graph(n,SAT)
        h , flips, its = g.solveSATclauseOnly(i,Tb)
        stepList.append(its)
#         print(stepList[-1])
    
#     print("Stat time to solve: " , start - time.time())
    
    return stepList
# def main():
#     testRSA()
    

# if __name__ == "__main__":
#     main()

    