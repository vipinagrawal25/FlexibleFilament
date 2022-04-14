import header
#------------------------------------------------------------------------------------#
index = -2
datadir="Zn0.75/"
f,ax=plt.subplots()
mc_log=np.loadtxt("./"+"Zn0.75/"+"/mc_log")
plt.plot(mc_log[:,index])
plt.show()