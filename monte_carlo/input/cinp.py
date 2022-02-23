import sys
sys.path.append("../tools/")
from MESH import *
import functiontools as ft
import numpy as np
from importlib import reload
import h5py
import dump_visit as dv
from scipy.spatial import SphericalVoronoi
#################################### FUNCTIONS ###################################
def SphVoronoi(rr,R=1,lplot=False):
    Np = np.shape(rr)[0]
    xyz = np.zeros([Np,3])
    for ip in range(Np):
        tht = rr[ip,0]
        phi = rr[ip,1]
        x = R*np.sin(tht)*np.cos(phi)
        y = R*np.sin(tht)*np.sin(phi)
        z = R*np.cos(tht)
        xyz[ip] = np.array([x,y,z])
    sv = SphericalVoronoi(xyz,radius=R)
    if lplot:
        plot_voronoi(xyz,sv)
    return sv
#################################### SCRIPT ###################################
filei = sys.argv[1]
rr = np.loadtxt(filei)
sv = SphVoronoi(rr)
# cumlst, node_neighbour, bond_neighbour = ft.neighbours(sv.points,sv._simplices)
msh = Mesh(sv.points,sv._simplices)
msh.assign_nbrs()
hf = h5py.File('input.h5','w')
hf.create_dataset('pos',data=sv.points)
hf.create_dataset('triangles',data=sv._simplices)
hf.create_dataset('cumu_list',data=msh.cmlst)
hf.create_dataset('node_nbr',data=msh.node_nbr)
nbn = np.zeros([len(msh.node_nbr),2], dtype=int)
for i,bn in enumerate(msh.bond_nbr):
    nbn[i,0] = bn[0]
    nbn[i,1] = bn[1]
hf.create_dataset('bond_nbr',data=nbn)
hf.close()