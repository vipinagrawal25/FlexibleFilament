import os as os
# import pencil_old as pc
import numpy as np
import matplotlib.pyplot as P
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
#%matplotlib notebook
from importlib import reload
import h5py
import yaml
#from stripack import trmesh
import FuncMC2d as F
import dump_visit as dv
from FuncMC2d import MESH
#
P.style.use('matplotlibrc')
#
## moved all the styling to matplotlibrc
#---------------------------------------------#
with open(r'input.yaml') as file:
    """input for the run"""
    inp = yaml.full_load(file)
#F.triang_moebius()
#F.triang_sph()
Np=inp['num_particles']
#
debug = inp['debug']
if debug:
    rrini = np.loadtxt('initial_rrini.dat')
    np.savetxt('initial_rrini.dat', rrini)
else:
    rrini = F.rand_sph(Np)
    np.savetxt('initial_rrini.dat', rrini)
#
svini = F.SphVoronoi(rrini)
dv.dump_visit('initial_rand_points.vtk', svini.points, svini._simplices)
#
if(inp['surface_montecarlo']):
    rr = F.MC_surf(Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=10000,kBT=1.,
                                 dfac=Np,interactive=False)
    np.savetxt('final_rrini.dat', rr)
else:
    rr = rrini
#
if(inp['read_ini_particle']):
    hf=h5py.File("fin_pos.h5","r")
    rr=np.array(hf.get('rr'))
    hf.close()
sv = F.SphVoronoi(rr)
cmlst,node_neighbour,bond_neighbour=F.neighbours(sv)
#
mesh=MESH(Np=Np,
        R=sv.points,
        BB=1,
        HH=1,
        cmlst=cmlst,
        node_nbr=node_neighbour,
        bond_nbr=bond_neighbour)
#
dv.dump_visit('initial_points.vtk', mesh.R, sv._simplices)
#
mesh=F.MC_mesh(mesh,maxiter=10000,kBT=1.,interactive=True,dfac=64)
#
dv.dump_visit('final_points.vtk', mesh.R, sv._simplices)