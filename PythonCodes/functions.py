from pylab import *
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.fftpack import dst
from scipy.interpolate import interp1d

def curvatureplot(FILE='output',omega=3,length=1,tmin=0,tmax=-1):
	dd = loadtxt(FILE+'/kappa.txt')
	dd_SS = loadtxt(FILE+'/material_point.txt')
	Y = dd_SS[tmin:tmax]
	nrowcol=dd_SS.shape
	nsnap=nrowcol[0]
	Np = nrowcol[1]
	# Now normalize the material point coordinate
	for isnap in range(0,nsnap-1,1):
		Y[isnap,:] = Y[isnap,:]*length/Y[isnap,-1]

	timeAxis = loadtxt(FILE+'/time.txt')
	timeAxis = timeAxis[tmin:tmax]
	IndexAxis = arange(Np)
	heightAxis = dd[tmin:tmax,1:Np+1]
	
	# heightAxis = ndimage.rotate(heightAxis,90)
	# imshow(heightAxis,cmap=plt.cm.gray)
	# colorbar()
	(IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
	fig = figure()
	ax = fig.add_subplot(111)
	surf = ax.contourf(X,Y,heightAxis)
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

def MSD_plot(FILE='',omega=1.5):
	file = loadtxt(FILE+'MSD.txt')
	tt = loadtxt(FILE+'output/time.txt')

	ttsclf = (tt[2]-tt[1])*omega/(2*pi)
	icy = int(1/ttsclf)

	cycle_index = arange(0,tt.size,icy)
	plt.plot(cycle_index,file[cycle_index],'.-') 
	# plt.plot(file)
	# plt.ylim((0, 0.1))
	plt.savefig('MSD_cycle.eps')
	plt.show()

# This function calculates Mean square displacement and substracts the translation of the rod after every cycle if there is any.
# For this, basically I substract the co-ordinate of middle point.
def MSD_no_trans(FILE='',omega=3):
	tt = loadtxt(FILE+'output/time.txt')
	dd_ini = loadtxt(FILE+'output/var0.txt')
	nrowcol=dd_ini.shape
	Np = nrowcol[0]
	
	ttsclf = (tt[2]-tt[1])*omega/(2*pi)
	dcy = int(1/sclf)
	cycle_index = arange(0,tt.size,icy)
	
	MSDnew=zeros(floor(tt.size/icy))

	if Np%2==0:
		MidP=int(Np/2)
	else:
		MidP=(Np+1)/2
		MidP=int(MidP)

	for idip in range(0,ndips):
		isnap=dips[idip]
		dd=loadtxt(FILE+'output/var'+str(isnap)+'.txt')
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

def PowerSpec(hh,Delta,deali=True):
    hh=hh-NP.mean(hh)
    NN = int(hh.size)
    NNBy2 = int(NN/2)
    # Check whether hh is one dimensional, or power of 2
    ff = (1/(NN*Delta))*arange(0,NNBy2+1,dtype=int)
    HH = fft.rfft(hh)
    Pxx = abs(HH)
    
    if(deali):
        ff = ff[0:int(0.65*ff.size)] 		#De-aliasing the data
        Pxx = Pxx[0:int(0.65*Pxx.size)]
        HH = HH[0:int(0.65*HH.size)]
    
    return ff, HH, Pxx

def FindPeaks(function,facTh=0.1,neighbours=2):
	if neighbours==2:
		pks = where((function[2:-2] > function[0:-4]) * (function[2:-2] > function[1:-3]) 
	              * (function[2:-2] > function[3:-1])  * (function[2:-2] > function[4:]))[0] + 2
	else:
		pks = where((function[1:-1] > function[0:-2]) * (function[1:-1] > function[2:]) )[0] + 1
	
	MaxFunc = max(function)
	threshold = facTh*MaxFunc
	peaks = []
	if ( (function[1]>function[0])*(function[1]>function[2]) ):
	    peaks.append(1)

	for ipk in pks:
	    if (function[ipk]> threshold):
	        peaks.append(ipk)
	        
	return peaks

def CurvatureSign(kappasqr,yy,zz,eps):
	NN=yy.size
	sign=zeros(NN)
	znew=zeros(NN)
	# First shift everything by the linear line joining the first and last line.
	yy=yy-(yy[0]-yy[-1])/(zz[0]-zz[-1])*zz
	# Now interpolate things
	for iN in range(0,NN):
		znew[iN] = zz[0]+(zz[-1]-zz[0])*iN/NN
	ff = interp1d(zz,yy)
	ynew=ff(znew)
	# Assign sign based on the slope in starting. Remember first two points have zero curvature so just go two more points
	# to make sure slope calculation is right and implement it.
	ynewdiff=diff(ynew)
	ynewdiff=ynewdiff[0:3]
	ynewdiff[ynewdiff<0]=-1
	ynewdiff[ynewdiff>0]=1

	sign[1:4]=ynewdiff
	# Now initial sign has been defined. we can change the sign whenever double derivative i.e. kappasqr->0
	dips=where((kappasqr[4:-1] < kappasqr[3:-2]) * (kappasqr[4:-1] < kappasqr[5:]) )[0] + 4
	# Sign should continue till the first minimum is hit
	if dips.size>0:
		sign[4:dips[0]+1]=sign[3]
	else:
		sign[4:]=sign[3]
	# Now let's change the sign alternatively if minimum is below some epsilon value
	for index in range(1,dips.size):
		if kappasqr[dips[index]]<eps:
			sign[dips[index-1]+1:dips[index]+1]=-sign[dips[index-1]]
		else:
			sign[dips[index-1]+1:dips[index]+1]=sign[dips[index-1]]
	# Sign should just continue after the last minimum
	if dips.size>0:
		sign[dips[-1]+1:]=-sign[dips[-1]]
	# Multiply signs to the absolute value of kappasqr
	return sign

def GetCurv(Folder='output/',code='CPU',eps=0.04):
	# This function will take the square root of curvature. Sign of the final thing would be decided
	# by double derivative of the position vector. If the function is convex, curvature can be negative
	# else the curvature would be positive.	
	dd = loadtxt(Folder+'curvature.txt')
	dd[dd<0]=0
	nrowcol = dd.shape
	NN = nrowcol[1] -1
	nsnap = nrowcol[0]
	time=dd[:,0]
	curvsqr = dd[:,1:NN+1]
	kappa=zeros([nsnap,NN+1])
	kappa[:,0]=time
	if (code=='CPU'):
		for isnap in range(1,nsnap):
			dd = loadtxt(Folder+'var'+str(isnap)+'.txt')
			zz=dd[:,2]
			yy=dd[:,1]
			kappa[isnap,1:NN+1]=sqrt(curvsqr[isnap,:])*CurvatureSign(curvsqr[isnap,:],yy,zz,eps)
	elif (code == 'GPU'):
		dd=loadtxt(Folder+"PSI")
		zz=zeros(NN)
		yy=zeros(NN)
		for isnap in range(1,nsnap):
			for iN in range(0,NN):
				yy[iN] = dd[isnap,3*iN+1]
				zz[iN] = dd[isnap,3*iN+3]
			kappa[isnap,1:NN+1]=sqrt(curvsqr[isnap,:])*CurvatureSign(curvsqr[isnap,:],yy,zz,eps)
	
	savetxt(Folder+'kappa.txt',kappa,fmt='%.5e')
