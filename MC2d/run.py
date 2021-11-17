import os as os
import pencil_old as pc
import numpy as np
import matplotlib.pyplot as P
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
#%matplotlib notebook
from importlib import reload
import h5py
#from stripack import trmesh
import FuncMC2d as F
#-----------------------------
P.rc('font',size=22)
P.rc('figure',figsize=(8,6))
P.rc('figure.constrained_layout',use=True)
P.rc('xtick',direction='in')
P.rc('ytick',direction='in')
P.rc('xtick.major',size=8,width=2)
P.rc('xtick.minor',visible=True,size=4,width=1)
P.rc('ytick.major',size=8,width=2)
P.rc('ytick.minor',visible=True,size=4,width=1)
#---------------------------------------------#
#F.triang_moebius()
#F.triang_sph()
Np=16
rrini = F.rand_sph(Np)
print(rrini)
#F.MC_surf(Np,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',maxiter=1000,kBT=1.,
#          dfac=Np,interactive=True)
rr = F.MC_surf(Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=1000,kBT=1.,
                     dfac=Np,interactive=False)
sv = F.SphVoronoi(rr)
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

