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
writer = FFMpegWriter(fps=30, metadata=metadata)
fig=plt.figure()
height=1

#####################################################
def MakePlot3D(ax,Xaxis,Yaxis,Zaxis,tsnap):
	# ax = fig.add_subplot(2,1,1, projection='3d')
	ax.plot(Yaxis,Xaxis,Zaxis,'o-')
	ax.set_xlim(-1, 1)
	ax.set_ylim(-1, 1)
	ax.set_zlim(0,height)
	plt.title(str(tsnap))
	ax.set_xlabel('Y direction',fontsize=12)
	ax.set_ylabel('X direction',fontsize=12)
	ax.set_zlabel('Z direction',fontsize=12)
	ax.grid(True)
#####################################################
def MakePlot2D(ax,Xaxis,Yaxis,tsnap):
	# ax.subplot(2,1,2)
	ax.set_xlim(-0.5,1.5)
	ax.set_ylim(-0.3,1.3)
	ax.plot(Xaxis,Yaxis,'o-')
	plt.title(str(tsnap))
	# plt.hold(False)
	# ax.grid(True)
#####################################################
def GetCoordinate(arr,dim):
	row,col = arr.shape
	xx = zeros(row)

	if col==6 or col==3:	# It means that the simulation was run in 3 dimension
		xx=arr[:,0]
		yy=arr[:,1]
		zz=arr[:,2]
	elif col==4:
		yy = arr[:,0]
		zz = arr[:,1]
	return xx,yy,zz
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
	# nsnap=600
	with writer.saving(fig,"movie.mp4", 100):
		for isnap in range(1,nsnap,1):
	 		dd=loadtxt(FILE+'/var'+str(isnap)+'.txt')
	 		xx,yy,zz=GetCoordinate(dd,dim)
	 		# xx=dd[:,0]
	 		# yy=dd[:,1]
	 		# zz=dd[:,2]

	 		if dim==2:
	 			MakePlot2D(ax,yy,zz,time[isnap])
	 		else:
	 			MakePlot3D(ax,xx,yy,zz,time[isnap])

	 		writer.grab_frame()
	 		ax.clear()
	 		print('plot = ' + str(isnap)+ 'Done')

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