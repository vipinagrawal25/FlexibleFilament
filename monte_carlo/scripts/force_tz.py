import numpy as np
import matplotlib.pyplot as plt
import sys
import os
Home=os.getenv('Home')
ES=os.getenv('ES')
sys.path.append(ES+'/monte_carlo/tools/')
import functiontools as ft
plt.style.use(Home+'/.matplotlibrc')
plt.rcParams['text.usetex'] = True
import glob

tzall=sorted(glob.glob("*/mc_log"))
folders=[tz.replace("/mc_log","/") for tz in tzall]
folders.reverse()
start=0
dvert,d_err=ft.dvert(folders=folders,start=start,nopush='nopull/',subfol="./")
dvert,force,baseE,f_error=ft.avg_quantity_tz(folders=folders,index=-4,datadir="",subfol="./",
                        nopush="nopull/",index2=-5,start=start)
# np.savetxt('scripts/dvert_force.txt',"# dvert d_err force f_error\n")
np.savetxt('scripts/dvert_force.txt',np.transpose([dvert,d_err,force,f_error]))
# Ks = ft.Ks_spring(start=start,datadir=exo65+"1_BB2/Zn1.0/")
# print(Ks)
f,ax=plt.subplots()
# plt.plot(dvert,force,'o')
plt.errorbar(dvert,force,yerr=f_error,fmt='o')
ax.set(xlabel=r'1-$\frac{\langle z \rangle}{\langle z_0 \rangle}$',ylabel=r'Force*radius/KbT')
plt.show()