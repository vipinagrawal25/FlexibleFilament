from numpy import *
import matplotlib
matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import sys


FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
# <<<<<<< HEAD
# writer = FFMpegWriter(fps=10, metadata=metadata)

# fig = plt.figure()
# # ax1=fig.add_subplot(1,1,1,projection='3d')
# ax2=fig.add_subplot(1,1,1)

# FILE = 'output'
# height = 2

# def MakePlot3D(ax,Xaxis,Yaxis,Zaxis,isnap):
# =======
writer = FFMpegWriter(fps=30, metadata=metadata)
fig=plt.figure()
height=1
#####################################################
def MakePlot3D(ax,Xaxis,Yaxis,Zaxis,tsnap,symbol='.',colors=[],colormap='gnuplot'):
# >>>>>>> NewtonKrylov
	# ax = fig.add_subplot(2,1,1, projection='3d')
	if colors==[]:
		colors = ones(Xaxis.shape[0])
	ax.plot(Yaxis,Xaxis,Zaxis,symbol)
	# ax.set_xlim(-1, 1)
	ax.set_xlim(-0.2, 0.2)
	ax.set_zlim(0,height)
	plt.title(str(tsnap))
	ax.set_xlabel('Y direction',fontsize=12)
	ax.set_ylabel('X direction',fontsize=12)
	ax.set_zlabel('Z direction',fontsize=12)
	ax.grid(True)
#####################################################
def MakePlot2D(ax,Xaxis,Yaxis,tsnap,symbol='.'):
	# ax.subplot(2,1,2)
# <<<<<<< HEAD
# 	# ax.set_xlim(-1,2*height)
# 	ax.set_ylim(-1,1.3)

# 	ax.plot(Xaxis,Yaxis,'o-')
# 	plt.title(str(time[isnap]))
# =======
	ax.set_xlim(-2,2)
	ax.set_ylim(-2,2)
	ax.plot(Xaxis,Yaxis,symbol)
	plt.title(str(tsnap))
# >>>>>>> NewtonKrylov
	# plt.hold(False)
	# ax.grid(True)
#####################################################
def GetCoordinate(arr,dim,nparticles=256):
	row,col = arr.shape
	xx = zeros(row)
	if col==6 or col==3:						#It means that the simulation was run in 3 dimension
		xx=arr[0:nparticles,0]
		yy=arr[0:nparticles,1]
		zz=arr[0:nparticles,2]
	elif col==4:
		yy = arr[0:nparticles,0]
		zz = arr[0:nparticles,1]
	return xx,yy,zz
#####################################################
# def GetCoordinate2(arr):
# 	row = arr.shape
# 	xx = zeros(row)
# 	yy = zeros(row)
# 	zz = zeros(row)
# 	for ir in range(row):
# 		yy[ir] = 
# 		zz[ir] = 
# 	return xx,yy,zz
#####################################################
def MultiFileMovie(FILE='output',dim=3):
	if dim==2:
		ax=fig.add_subplot(1,1,1)
	elif dim==3:
		ax=fig.add_subplot(1,1,1,projection='3d')
	else:
		sys.exit('The dimension does not exist.')
	time = loadtxt(FILE+'/time.txt')
	nrowcol=time.shape
	nsnap=nrowcol[0]
	with writer.saving(fig,"movie.mp4", 100):
		for isnap in range(0,nsnap,1):
	 		dd=loadtxt(FILE+'/var'+str(isnap)+'.txt')
	 		xx,yy,zz=GetCoordinate(dd,dim)
	 		if dim==2:
	 			MakePlot2D(ax,yy,zz,time[isnap],'o')
	 		else:
	 			MakePlot3D(ax,xx,yy,zz,time[isnap],'o')
	 		writer.grab_frame()
	 		ax.clear()
	 		print('plot = ' + str(isnap)+ 'Done')
#####################################################
def TracerMovie(FILE='output',dim=3,symbol='o'):
	if dim==2:
		ax=fig.add_subplot(1,1,1)
	elif dim==3:
		ax=fig.add_subplot(1,1,1,projection='3d')
	else:
		sys.exit('The dimension does not exist.')
	time = loadtxt(FILE+'/time.txt')
	nrowcol=time.shape
	nsnap=nrowcol[0]
	with writer.saving(fig,"movie_tracer.mp4", 100):
		for isnap in range(1,nsnap,1):
	 		dd=loadtxt(FILE+'/tracer'+str(isnap)+'.txt')
	 		xx,yy,zz=GetCoordinate(dd,dim)
	 		if dim==2:
	 			MakePlot2D(ax,yy,zz,time[isnap],symbol)
	 			# MakePlot2D_text(ax,yy,zz,time[isnap])
	 		else:
	 			MakePlot3D(ax,xx,yy,zz,time[isnap])
	 		writer.grab_frame()
	 		ax.clear()
	 		print('plot = ' + str(isnap)+ 'Done')
#####################################################
def TracerMovieCycle(omega,FILE='output',dim=3,symbol='o',nparticles=256):
	if dim==2:
		ax=fig.add_subplot(1,1,1)
	elif dim==3:
		ax=fig.add_subplot(1,1,1,projection='3d')
	else:
		sys.exit('The dimension does not exist.')
	time = loadtxt(FILE+'/time.txt')
	nrowcol=time.shape
	nsnap=nrowcol[0]
	TT = 2*pi/omega
	ncy = int(nsnap/(TT*50))
	with writer.saving(fig,"movie_tracer_cy.mp4", 100):
		for icy in range(0,ncy,1):
	 		isnap = int(icy*TT*50)
	 		dd=loadtxt(FILE+'/tracer'+str(isnap)+'.txt')
	 		xx,yy,zz=GetCoordinate(dd,dim,nparticles)
	 		if dim==2:
	 			MakePlot2D(ax,yy,zz,time[isnap],symbol)
	 			# MakePlot2D_text(ax,yy,zz,time[isnap])
	 		else:
	 			MakePlot3D(ax,xx,yy,zz,time[isnap])
	 		writer.grab_frame()
	 		ax.clear()
	 		print('plot = ' + str(isnap)+ 'Done')
# #####################################################
# def TracerMovie(cy1,cy2,omega,FILE='output',dim=3,symbol='o'):
# 	if dim==2:
# 		ax=fig.add_subplot(1,1,1)
# 	elif dim==3:
# 		ax=fig.add_subplot(1,1,1,projection='3d')
# 	else:
# 		sys.exit('The dimension does not exist.')
# 	time = loadtxt(FILE+'/time.txt')
# 	TT = 2*pi/omega;
# 	snap1 = int(50*TT*cy1)
# 	snap2 = int(50*TT*cy2)
# <<<<<<< HEAD
# with writer.saving(fig,"movie.mp4", 100):
# 	for isnap in range(1,itn,1):
# 		dd = loadtxt(FILE+'/position'+str(isnap)+'.txt')
# 		xx = dd[:,0]
# 		yy = dd[:,1]
# 		zz = dd[:,2]
# 		# MakePlot3D(ax1,xx,yy,zz,isnap)
# 		MakePlot2D(ax2,yy,zz,isnap)
# 		writer.grab_frame()
# 		ax2.clear()
# 		# ax1.clear()
# 		print('plot = ' + str(isnap)+ 'Done')
# =======
# 	with writer.saving(fig,"movie_tracer"+str(cy2-cy1)+".mp4", 100):
# 		for isnap in range(snap1,snap2,1):
# 	 		dd=loadtxt(FILE+'/tracer'+str(isnap)+'.txt')
# 	 		xx,yy,zz=GetCoordinate(dd,dim)
# 	 		if dim==2:
# 	 			MakePlot2D(ax,yy,zz,time[isnap],symbol)
# 	 			# MakePlot2D_text(ax,yy,zz,time[isnap])
# 	 		else:
# 	 			MakePlot3D(ax,xx,yy,zz,time[isnap])
# 	 		writer.grab_frame()
# 	 		ax.clear()
# 	 		print('plot = ' + str(isnap)+ 'Done')
#####################################################
def SingleFileMovie(FILE='data/PSI',dim=2, output="data/movie.mp4"):
	dd=loadtxt(FILE)
	nrowcol=dd.shape
	nrows=nrowcol[0]
	# nrows=400
	ncols=nrowcol[1]
	# print(ncols)
	NN=int((ncols-1)/3)
	#
	if dim==2:
		ax=fig.add_subplot(1,1,1)
	elif dim==3:
		ax=fig.add_subplot(1,1,1,projection='3d')
	else:
		sys.exit('The dimension does not exist.')
	#
	with writer.saving(fig,output,100):
		for irow in range(0,nrows):
			xx=zeros(NN)
			yy=zeros(NN)
			zz=zeros(NN)
			#
			time=dd[irow,0]
			for icol in range(0,NN):
				xx[icol]=dd[irow,3*icol+1]
				yy[icol]=dd[irow,3*icol+2]
				zz[icol]=dd[irow,3*icol+3]

			if dim==2:
				MakePlot2D(ax,xx,zz,time)
			else:
				MakePlot3D(ax,yy,xx,zz,time)

			writer.grab_frame()
			ax.clear()
			print('plot = ' + str(irow)+ 'Done')
# >>>>>>> NewtonKrylov