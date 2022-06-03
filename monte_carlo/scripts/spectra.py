# from header import *
# lmax=22
# hlm=np.zeros(lmax)
# istart=10	
# iend=3000
# for foln in range(istart,iend,10):
# 	fol1='output_fig2_yellow/part_'+str(foln).zfill(5)+'.vtk'
# 	hlm=hlm+ft.spectra(fol1,lmax=lmax,Np=23042)
# 	print(fol1)
# plt.plot(hlm*hlm,'o-')
# plt.show()
from header import *
lmax=25
istart=100
iend=3000
for iarg in range(1,len(sys.argv)):
	hlm=np.zeros(lmax)
	for foln in range(istart,iend,10):
		fol=sys.argv[iarg]+'/part_'+str(foln).zfill(5)+'.vtk'
		hlm=ft.spectra(fol,lmax=lmax,Np=23042)
		print(fol)
		xx = np.arange(lmax)
		np.savetxt(sys.argv[iarg]+"spectra/spec_"+str(foln).zfill(5),np.vstack([xx,hlm]).T)
		# hlm2=hlm*hlm
		# plt.plot(np.arange(1,lmax),hlm2[1:]/hlm2[0],'o-',label=fol)
# plt.legend()
# plt.show()