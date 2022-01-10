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
from FuncMC2d import MESH
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

if(inp['do_montecarlo']):
    rr = F.MC_surf(Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=1000,kBT=1.,
                                 dfac=Np,interactive=False)
else:
    rr = rrini
#
if(inp['read_ini_particle']):
    hf=h5py.File("fin_pos.h5","r")
    rr=np.array(hf.get('rr'))
    hf.close()
sv = F.SphVoronoi(rr)
cumlst,node_neighbour,bond_neighbour=F.neighbours(sv)
#
mesh=MESH(Np,sv.points,cumlst,node_neighbour,bond_neighbour)
print(mesh.dis(1,2))