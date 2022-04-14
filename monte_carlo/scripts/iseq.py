from header import *
#------------------------------------------------------------------------------------#
index1 = -5
index2 = -4
datadir="Zn0.5/"
f,ax=plt.subplots()
mc_log=np.loadtxt(datadir+"/mc_log")
plt.plot(0.5*(mc_log[:,index1]-mc_log[:,index2]))
ft.isgaussian(datadir=datadir,start=10000)
plt.show()