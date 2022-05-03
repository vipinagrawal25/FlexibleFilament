from header import *
tzall=sorted(glob.glob("*/mc_log"))
folders=[tz.replace("/mc_log","/") for tz in tzall]
folders.reverse()
start=800000
#
mc_log = ft.read_mc_log(folders=folders,start=start)
dvert,d_err=ft.dvert(mc_log=mc_log[1:][:],nopush='nopull/',subfol="./")
force,f_error=ft.avg_quantity_tz(mc_log=mc_log[1:][:],index=-6,index2=-7)
Ks,Ks_err=ft.Ks_spring(mc_log=mc_log,datadir="",subfol="./",start=start)
np.savetxt('scripts/dvert_force.txt',np.transpose([dvert,d_err,force,f_error]))
np.savetxt('scripts/Ks.txt',np.transpose([Ks,Ks_err]))
f,ax=plt.subplots()
plt.errorbar(dvert,force,yerr=f_error,fmt='o')
ax.set(xlabel=r'1-$\frac{\langle z \rangle}{\langle z_0 \rangle}$',ylabel=r'Force*radius/KbT')
plt.show()