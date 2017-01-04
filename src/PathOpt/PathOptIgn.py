import numpy
import random

ns = 7
CostMax = 30
dist = numpy.zeros((ns, ns))
cost = numpy.zeros((ns, ns))
tours = numpy.zeros((ns+1, ns))
tempi = [3,4,2,6,7,8,9]
for i in range(ns):
    for j in range(ns):
        dist[i,j] = 11.*random.random()
for i in range(ns):
    dist[i,i] = 0.
for i in range(ns):
    for j in range(ns):
        cost[i,j] = dist[i,j]*tempi[j-i]
for i in range(ns):
    tours[1,i] = i+1
tInit = ns
IsImproved = 1
savings = numpy.zeros((ns, ns))
estremi = numpy.zeros((2, ns))
col = numpy.zeros((ns, 1))
#trova estremi
for i in range(ns):
    estremi[0,i] = cost[i,i]
    for j in range(ns):
        col[j,0] = cost[j,i]
    k = 0
    for j in range(ns):
        if(tours[i,0] == 0):
            k = tours[i-1,0]
    estremi[1,i]=k
#riempi costi            
for i in range(ns):
    ext1 = estremi[0,i]
    for j in range(ns):
        ext2 = estremi[1,j]
        savings[i,j] = cost[0,ext1]+cost[ext2,0]-cost[ext2,ext1]
for i in range(ns):
    savings[i,i] = 0

print cost

