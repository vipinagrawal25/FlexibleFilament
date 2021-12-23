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
P.style.use('matplotlibrc')

## moved all the styling to matplotlibrc
#---------------------------------------------#
with open(r'input.yaml') as file:
    """input for the run"""
    inp = yaml.load(file, Loader=yaml.FullLoader)
#F.triang_moebius() 
#F.triang_sph()

Np=inp['num_particles']

rrini = F.rand_sph(Np)

debug = inp['debug']

if(debug):
    F.WritePos(rrini,fname='ini_pos')

if(inp['do_montecarlo']):
    rr = F.MC_surf(Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=1000,kBT=1.,
                                 dfac=Np,interactive=False)
else:
    rr = rrini

if(inp['read_ini_particle']):
    hf=h5py.File("fin_pos.h5","r")
    rr=np.array(hf.get('rr'))
    hf.close()


sv = F.SphVoronoi(rr)
F.assign_newmems(sv);
if(debug):
    print(type(sv))
    print(sv.regions)
    print(sv.tris)

# print(F.normal(sv,1))
# cn_NNL,NNL = F.nearest_neighbour(sv)
# print("Nearst neighbour list: ")
# print(NNL)
# print("Area = "+str(np.sum(sv.calculate_areas())))
#lats = rr[:,0]-np.pi/2
#lons = rr[:,1]
#print('max longitude=',(np.abs(lons)).max())
#print('min latitude=',(np.abs(lats)).max()-np.pi/2)
#P.plot(lons, lats,'*')
#P.grid(True)
#P.show()
#print('lat(min,max)',np.min(lats),np.max(lats))
#print(np.pi/2)
#print("MC done, now doing triangulation")
#tri = trmesh(lons, lats)

#print(dir(tri))
#F.lat_lon_list(tri)
#print('writing x,y,z coordinates to a hdf5 file')
#F.print_xyz(lats,lons,fname='ini_sph',radius=1.)
#fig = P.figure()
#ax = fig.add_subplot(111)
#ax = F.plot_pos(ax,rr,rr)
#P.show()
