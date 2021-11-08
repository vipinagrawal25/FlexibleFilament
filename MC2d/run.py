import os as os
import pencil_old as pc
import numpy as np
import matplotlib.pyplot as P
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
#%matplotlib notebook
from importlib import reload
import h5py
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
Np=64
#rrini = F.rand_sph(Np)
#print(rrini)
#F.MC_surf(Np,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',maxiter=1000,kBT=1.,
#          dfac=Np,interactive=True)
F.MC_surf(Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=1000,kBT=1.,
          dfac=Np,interactive=True)
#fig = P.figure()
#ax = fig.add_subplot(111)
#ax = F.plot_pos(ax,rr,rr)
#P.show()

