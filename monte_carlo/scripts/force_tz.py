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
dvert,force,baseE,error=ft.avg_quantity_tz(folders=folders,index=-4,datadir="./",subfol="./",
										   nopush='Zn1.0/')
f,ax=plt.subplots()
plt.plot(dvert,force,'o-')
plt.errorbar(dvert,force,yerr=error,fmt='o')
ax.set(xlabel=r'1-$\frac{\langle z \rangle}{\langle z_0 \rangle}$',ylabel=r'Force*radius/KbT')
plt.show()