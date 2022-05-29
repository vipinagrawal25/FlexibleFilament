from header import *
lmax=22
hlm=np.zeros(lmax)
istart=2000
iend=2010
for foln in range(istart,iend,10):
	fol1='output_fig2_yellow/part_'+str(foln).zfill(5)+'.vtk'
	hlm=hlm+ft.spectra(fol1,lmax=lmax,Np=23042)
	print(fol1)
plt.plot(hlm*hlm,'o-')
plt.show()