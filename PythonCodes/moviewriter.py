from numpy import *
import matplotlib
matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

itn = 2000
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=10, metadata=metadata)

fig = plt.figure()
# ax1=fig.add_subplot(1,1,1,projection='3d')
ax2=fig.add_subplot(1,1,1)

FILE = 'output'
height = 2

def MakePlot3D(ax,Xaxis,Yaxis,Zaxis,isnap):
	# ax = fig.add_subplot(2,1,1, projection='3d')
	ax.plot(Xaxis,Yaxis,Zaxis,'.-')
	ax.set_xlim(-1, 1)
	ax.set_ylim(-0.5, 0.5)
	ax.set_zlim(0,height)
	plt.title(str(time[isnap]))

def MakePlot2D(ax,Xaxis,Yaxis,isnap):	
	# ax.subplot(2,1,2)
	# ax.set_xlim(-1,2*height)
	ax.set_ylim(-1,1.3)

	ax.plot(Xaxis,Yaxis,'o-')
	plt.title(str(time[isnap]))
	# plt.hold(False)
	ax.grid(True)

time = loadtxt(FILE+'/time.txt')

with writer.saving(fig,"movie.mp4", 100):
	for isnap in range(1,itn,1):
		dd = loadtxt(FILE+'/position'+str(isnap)+'.txt')
		xx = dd[:,0]
		yy = dd[:,1]
		zz = dd[:,2]
		# MakePlot3D(ax1,xx,yy,zz,isnap)
		MakePlot2D(ax2,yy,zz,isnap)
		writer.grab_frame()
		ax2.clear()
		# ax1.clear()
		print('plot = ' + str(isnap)+ 'Done')
