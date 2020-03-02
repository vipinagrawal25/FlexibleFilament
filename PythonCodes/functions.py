from pylab import *
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.fftpack import dst

def curvatureplot(FILE='output',omega=3):
	dd = loadtxt(FILE+'/curvature.txt')
	dd_SS = loadtxt(FILE+'/material_point.txt')
	Y = dd_SS
	nrowcol=dd_SS.shape
	nsnap=nrowcol[0]
	Np = nrowcol[1]
	# Now normalize the material point coordinate
	for isnap in range(0,nsnap-1,1):
		Y[isnap,:] = Y[isnap,:]/Y[isnap,-1]

	timeAxis = loadtxt(FILE+'/time.txt')
	IndexAxis = arange(Np)
	heightAxis = dd[:,1:Np+1]

	# heightAxis = ndimage.rotate(heightAxis,90)
	# imshow(heightAxis,cmap=plt.cm.gray)
	# colorbar()
	(IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
	fig = figure()
	ax = fig.add_subplot(111)
	# surf = ax.contourf(X,Y,heightAxis,vmin=10**(-5),vmax=10**5,cmap='plasma')
	surf = ax.contourf(X,Y,heightAxis,cmap='plasma')
	cbar = colorbar(surf)
	ax.set_aspect(20)
	ax.set_xlabel('time')
	ax.set_ylabel('material_point')
	plt.show()
	plt.savefig('curvatureplot.eps')
	close()

def GPUcurvatureplot(FILE='data',omega=3):
	heightAxis = loadtxt(FILE+'/curvature.txt')
	dd_SS = loadtxt(FILE+'/material_point.txt')
	nrowcol=dd_SS.shape
	nsnap=nrowcol[0]
	Np = nrowcol[1]

	for isnap in range(0,nsnap):
		for ip in range(1,Np):
			dd_SS[isnap,ip] = dd_SS[isnap,ip-1]+dd_SS[isnap,ip]

	Y = dd_SS
	# Now normalize the material point coordinate
	for isnap in range(0,nsnap):
		Y[isnap,:] = Y[isnap,:]/Y[isnap,-1]

	PSI = loadtxt(FILE+'/PSI')
	timeAxis = PSI[:,0]

	IndexAxis = arange(1,Np+1)

	# heightAxis = ndimage.rotate(heightAxis,90)
	# imshow(heightAxis,cmap=plt.cm.gray)
	# colorbar()
	(IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
	fig = figure()
	ax = fig.add_subplot(111)
	# surf = ax.contourf(X,Y,heightAxis,vmin=10**(-5),vmax=10**5,cmap='plasma')
	surf = ax.contourf(X,Y,heightAxis,cmap='plasma')
	cbar = colorbar(surf)
	ax.set_aspect(20)
	ax.set_xlabel('time')
	ax.set_ylabel('material_point')
	plt.show()
	plt.savefig('curvatureplot.eps')
	close()

def SineCurvature(FILE='output',aspect_ratio=1):
	dd = loadtxt(FILE+'/curvature.txt')
	dd_SS = loadtxt(FILE+'/material_point.txt')
	Y = dd_SS
	nrowcol=dd_SS.shape
	nsnap=nrowcol[0]
	Np=nrowcol[1]
	# Now normalize the material point coordinate
	for isnap in range(0,nsnap-1,1):
		Y[isnap,:] = Y[isnap,:]/Y[isnap,-1]

	Z = dd[:,1:Np+1]		#Removed the first entry
	Zsine=zeros((nsnap,Np))
	#
	for isnap in range(0,nsnap):
		Zsine[isnap,:] = dst(Z[isnap,:],type=1)

	timeAxis = loadtxt(FILE+'/time.txt')
	IndexAxis = arange(Np)
	# heightAxis = ndimage.rotate(heightAxis,90)
	# imshow(heightAxis,cmap=plt.cm.gray)
	# colorbar()
	(IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
	fig = figure()
	ax = fig.add_subplot(111)
	# surf = ax.contourf(X,Y,heightAxis,vmin=10**(-5),vmax=10**5,cmap='plasma')
	surf = ax.contourf(X,Y,Zsine)
	cbar = colorbar(surf)
	# ax.set_aspect(aspect_ratio)
	ax.set_xlabel('time')
	ax.set_ylabel('material_point')
	plt.show()
	plt.savefig('sinecurvature.eps')
	close()

def MSD_plot(FILE='MSD.txt'):
	file = loadtxt(FILE)
	dips = np.where((file[1:-1] < file[0:-2]) * (file[1:-1] < file[2:]))[0] + 1
	# plt.plot(file)
	plt.plot(dips,file[dips],'o')
	# plt.ylim((0, 0.1))
	plt.savefig('MSD_cycle.eps')
	plt.show()

# This function calculates Mean square displacement and substracts the translation of the rod after every cycle if there is any.
# For this, basically I substract the co-ordinate of middle point.
def MSD_no_trans(FILE=''):
	file = loadtxt(FILE+'MSD.txt')
	dips = np.where((file[1:-1] < file[0:-2]) * (file[1:-1] < file[2:]))[0] + 1		

	dd_ini = loadtxt(FILE+'output/position0.txt')
	nrowcol=dd_ini.shape
	Np = nrowcol[0]
	ndips = len(dips)
	MSDnew=zeros(ndips)

	if Np%2==0:
		MidP=int(Np/2)
	else:
		MidP=(Np+1)/2
		MidP=int(MidP)

	for idip in range(0,ndips):
		isnap=dips[idip]
		dd=loadtxt(FILE+'output/position'+str(isnap)+'.txt')
		dd[:,0]=dd[:,0]-dd[MidP,0]-dd_ini[:,0]
		dd[:,1]=dd[:,1]-dd[MidP,1]-dd_ini[:,1]
		dd[:,2]=dd[:,2]-dd[MidP,2]-dd_ini[:,2]
		#
		for ip in range(0,Np-1):
			MSDnew[idip]=MSDnew[idip]+dd[ip,0]*dd[ip,0]+dd[ip,1]*dd[ip,1]+dd[ip,2]*dd[ip,2]
			MSDnew[idip]=sqrt(MSDnew[idip])

	# plt.plot(dips,file[dips],'o')
	# plt.legend(['First Line'])
	plt.plot(dips,MSDnew,'.-')
	plt.ylim([1.2, 1.4])
	plt.legend(['MSD cycle','MSD after substracting the translation'])
	plt.savefig('MSD_without_translation.eps')
	close()

def Energy(AA,HH,FILE='output'):
	ddBE=loadtxt(FILE+'/curvature.txt')
	ddSE = loadtxt(FILE+'/material_point.txt')

	nrowcol=ddBE.shape
	nrows=nrowcol[0]
	ncols=nrowcol[1]
	#
	aa = 1./(ncols-2) 			# aa is the equlibrium distance
	BE=	zeros(nrows)			# Bending energy
	SE=zeros(nrows)				# Stretching energy
	TE=zeros(nrows)				# Total Energy
	time = zeros(nrows)
	#
	for i in range(0,nrows-1):
		time[i] = ddBE[i,0]
		for j in range(1,ncols):
			BE[i] = BE[i]+AA/aa*ddBE[i,j]
			SE[i] = SE[i]+HH/(2*aa)*(ddSE[i,j-1]-aa*(j-1))**2
			# TE[i] = BE[i]+SE[i]
	#
	# fig,ax = plt.subplots()
	# plt.plot(time,BE)
	# ax.plot(SE,label='Stretching Energy')
	# ax.plot(TE,label='Total Energy')
	
	# plt.xlabel('time')
	# plt.title('Bending Energy')
	# plt.savefig('Energy.eps')
	# plt.show()
	return time,BE,SE

def LeebyL(Np,Folder='output'):
	ddMP = loadtxt(Folder+'/material_point.txt')
	time = loadtxt(Folder+'/time.txt')
	nrowcol = ddMP.shape
	nsnap=nrowcol[0]

	LeebyL = zeros(nsnap)
	for isnap in range(0,nsnap,1):
		dd = loadtxt(Folder+'/var'+str(isnap)+'.txt')
		Lee = (dd[0,0]-dd[-1,0])**2 + (dd[0,1]-dd[-1,1])**2 + (dd[0,2]-dd[-1,2])**2
		Lee = sqrt(Lee)
		#
		L = 0
		for ip in range(0,Np-1):
			L = L+(dd[ip+1,0]-dd[ip,0])**2 + (dd[ip+1,1]-dd[ip,1])**2 + (dd[ip,2]-dd[ip+1,2])**2
		LeebyL[isnap] = Lee/L

	plt.plot(time,LeebyL,'o-')
	plt.plot(time,Lee,'o-')
	# plt.ylim((0, 0.1))
	plt.savefig('LeebyL.eps')
	plt.show()

	return LeebyL

# 	for idip in range(0,ndip-1):
# 		dd = loadtxt(FILE+'output/position'+str(dips(idip))+'.txt')
# 		dd = 
# 		MSD[idip] = 
	# plt.plot(file)	
	# plt.plot(dips,file[dips],'o')
	# plt.ylim((0, 0.1))
	# plt.savefig('MSD_cycle.eps')
	# plt.show()


# def BendingEnergy(AA,FILE='output/curvature.txt'):
# 	dd=loadtxt(FILE)
# 	nrowcol=dd.shape
# 	nrows=nrowcol[0]
# 	ncols=nrowcol[1]
# 	#
# 	aa = 1./(ncols-2) 					# aa is the equlibrium distance
# 	BE=zeros([nrows,1])
# 	time = zeros([nrows,1])
# 	#
# 	for i in range(0,nrows-1):
# 		time[i] = dd[i,0]
# 		for j in range(1,ncols):
# 			BE[i] = BE[i]+AA/aa*dd[i,j]
# 	#		
# 	plt.plot(time,BE)
# 	plt.title('Bending Energy')
# 	plt.xlabel('time')
# 	plt.ylabel('Bending Energy')
# 	plt.show()
# 	plt.savefig('BendingEnergy.eps')
# 	return time,BE

# def StretchEnergy(HH,FILE='output/material_point.txt'):
# 	dd = NP.loadtxt(FILE)
# 	nrowcol=dd.shape
# 	nrows=nrowcol[0]
# 	ncols=nrowcol[1]

# 	aa=1./(ncols-2)
# 	SE=NP.zeros([nrows,1])

# 	for i in range(0,nrows-1):
# 		for j in range(1,ncols):
# 			SE[i] = SE[i] + HH/(2*aa)*((dd[i,j]-aa)**2)

# 	plt.plot(SE)
# 	return SE

# chitra = plt.figure();
# def MSD_plot(FILE='MSD.txt',step=1, fig=chitra):
# 	step = int(around(step))
# 	# MSD = loadtxt(FILE+'/MSD.txt',delimiter=';',usecols=range(2000))
# 	MSD = loadtxt(FILE)
# 	# fig = figure()
# 	plt.plot(1/100*MSD[0::step],'o')
# 	# index = 0;
# 	plt.show()
# 	return fig
