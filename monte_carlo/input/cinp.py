import os as os
import sys
import numpy as np
import matplotlib.pyplot as P
from importlib import reload
import h5py
import yaml
#from stripack import trmesh
import FuncMC2d as F
import dump_visit as dv
from FuncMC2d import MESH
import math
import glob
#
# P.style.use('matplotlibrc')
#


filei = sys.argv[1] 
rr = np.loadtxt(filei)
sv=F.SphVoronoi(rr)
cumlst, node_neighbour, bond_neighbour = F.neighbours(sv.points,sv._simplices)

with open ("positions.dat", "w") as f:
    for p in sv.points:
        f.write("%lf %lf %lf\n" %(p[0], p[1], p[2]))

# def WritePos(rr,fname='input.h5'):
hf = h5py.File('input.h5','w')
hf.create_dataset('pos',data=sv.points)

with open ("triangles.dat", "w") as f:
    for p in sv._simplices:
        f.write("%d %d %d\n" %(p[0], p[1], p[2]))

hf.create_dataset('triangles',data=sv._simplices)

np.savetxt("cumlist.dat", cumlst, fmt="%d")

hf.create_dataset('cumu_list',data=cumlst)
np.savetxt("node_neighbour.dat", node_neighbour, fmt="%d")

hf.create_dataset('node_nbr',data=node_neighbour)
nbn = np.zeros([len(node_neighbour),2], dtype=int)
with open ("bond_neighbour.dat", "w") as f:
    for i,bn in enumerate(bond_neighbour):
        f.write("%d %d\n" %(bn[0], bn[1]))
        nbn[i,0] = bn[0]
        nbn[i,1] = bn[1]

print(np.shape(nbn))

hf.create_dataset('bond_nbr',data=nbn)


hf.close()
