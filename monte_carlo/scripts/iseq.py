from header import *
#------------------------------------------------------------------------------------#
index1 = -3
index2 = -2
datadir="nopull/"
f,ax=plt.subplots()
mc_log=np.loadtxt(datadir+"/mc_log")
plt.plot(0.5*(mc_log[:,index1]-mc_log[:,index2]))
ft.isgaussian(datadir=datadir,start=200000)
plt.show()
##### PDF of height flucutuations #################
# files=sorted(glob.glob(datadir+"/*.vtk"))
# hist = ft.height_field(files,Np=None,nbin=50)
# cen = hist[1][1:] + hist[1][:-1]
# plt.plot(cen*0.5,hist[0],'.-')