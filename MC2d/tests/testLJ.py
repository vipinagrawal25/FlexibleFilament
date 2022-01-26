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
P.style.use('../matplotlibrc')
#---------------------------------------------#
rad=1
Np=4096
lopt=np.sqrt(8*np.pi*rad**2/(2*Np-4))
sigma=lopt/(2**(1/6))
print(lopt)
#
for zz in np.linspace(0.015,0.05,1000):
	P.plot(zz,F.LJ(zz,sigma=sigma,epsilon=1), 'ok')
P.grid(True)
P.show()
#F.LJ(zz=np.range(0,4,100),sigma=1,zcut=4,eps_bot=4)