from header import *
lmax=25
istart=10
iend=0
for iarg in range(1,len(sys.argv)):
	if iend>istart:
		fols=[sys.argv[iarg]+'/part_'+str(foln).zfill(5)+'.vtk' for foln in range(istart,iend)]
	else:
		fols=sorted(glob.glob(sys.argv[iarg]+"/*.vtk"))
	hlm=np.zeros([len(fols),lmax])
	for i,fol in enumerate(fols):
		print(fol)
		hlm[i]=ft.spectra(fol,lmax=lmax,Np=5120)
	np.savetxt(sys.argv[iarg]+"/spec_l0.txt",hlm)