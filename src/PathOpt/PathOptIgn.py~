import numpy
import random

class PathOpt:
    ns = 7
    CostMax = 30
    dist = numpy.zeros((ns, ns))
    cost = numpy.zeros((ns, ns))
    tours = numpy.zeros((ns+1, ns))
    tempi = [3,4,2,6,7,8,9]
    def __init__(self):
        for i in range(self.ns):
            for j in range(self.ns):
                self.dist[i,j] = 11.*random.random()
        for i in range(self.ns):
            self.dist[i,i] = 0.
    def riempi_costi(self):
        for i in range(self.ns):
            for j in range(self.ns):
                self.cost[i,j] = self.dist[i,j]*self.tempi[j-i]
    def riempi_tours(self):
        for i in range(self.ns):
            self.tours[1,i] = i+1
            
    def sol_iniziale(self):
        tInit = self.ns
        IsImproved = 1
        savings = numpy.zeros((self.ns, self.ns))
        estremi = numpy.zeros((2, self.ns))
        col = numpy.zeros((self.ns, 1))
#trova estremi
        for i in range(self.ns):
            estremi[0,i] = self.cost[i,i]
            for j in range(self.ns):
                col[j,0] = self.cost[j,i]
            k = 0
            for j in range(self.ns):
                if(self.tours[i,0] == 0):
                    k = self.tours[i-1,0]
            estremi[1,i]=k
#riempi costi            
        for i in range(self.ns):
            ext1 = estremi[0,i]
            for j in range(self.ns):
                ext2 = estremi[1,j]
                savings[i,j] = self.cost[0,ext1]+self.cost[ext2,0]-self.cost[ext2,ext1]
        for i in range(self.ns):
            savings[i,i] = 0

    def show(self):
        print self.cost



if __name__ == "__main__":
    po = PathOpt()
    po.riempi_costi()
    po.riempi_tours()
    po.sol_iniziale()
    po.show()
