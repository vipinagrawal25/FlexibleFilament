from header import *
tzall=sorted(glob.glob("*/mc_log"))
folders=[tz.replace("/mc_log","/") for tz in tzall]
folders.reverse()
start=0
# dvert,d_err=ft.dvert(folders=folders,start=start,nopush='nopull/',subfol="./")
mc_log = ft.read_mc_log(folders=folders)
dvert,force,baseE,f_error=ft.avg_quantity_tz(mc_log=mc_log,index=-4,datadir="",subfol="./",
                        nopush="nopull/",index2=-5,start=start)
Ks,Ks_err=ft.Ks_spring(mc_log=mc_log,datadir="",subfol="./",nopush="nopull/",start=start)
# np.savetxt('scripts/dvert_force.txt',"# dvert d_err force f_error\n")
np.savetxt('scripts/dvert_force.txt',np.transpose([dvert,d_err,force,f_error]))
np.savetxt('scripts/Ks.txt',np.transpose([Ks,Ks_err]))
# Ks = ft.Ks_spring(start=start,datadir=exo65+"1_BB2/Zn1.0/")
# print(Ks)
f,ax=plt.subplots()
# plt.plot(dvert,force,'o')
plt.errorbar(dvert,force,yerr=f_error,fmt='o')
ax.set(xlabel=r'1-$\frac{\langle z \rangle}{\langle z_0 \rangle}$',ylabel=r'Force*radius/KbT')
plt.show()