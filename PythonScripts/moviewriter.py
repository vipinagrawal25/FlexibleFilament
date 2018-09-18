from numpy import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# from pylab import *
# import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation

itn = 3000
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

# l, = ax.plot([], [], [], 'o-' )
# m, = plt.title([])

#     
#     ax = fig.gca(projection='3d')
fig = plt.figure()
ax1=fig.add_subplot(1,1,1,projection='3d')
# ax2=fig.add_subplot(1,2,2)

def MakePlot3D(ax,Xaxis,Yaxis,Zaxis,isnap):
	# ax = fig.add_subplot(2,1,1, projection='3d')
	ax.plot(Xaxis,Yaxis,Zaxis,'o-')
	ax.set_xlim(-1, 1)
	ax.set_ylim(-0.5, 0.5)
	ax.set_zlim(0,1.3)
	plt.title(str(time[isnap]))

def MakePlot2D(ax,Xaxis,Yaxis,isnap):
	# ax.subplot(2,1,2)
	ax.plot(Xaxis,Yaxis,'o-')
	ax.set_xlim(-1,1.3)
	ax.set_ylim(-1,1.3)
	plt.title(str(time[isnap]))
	# plt.hold(False)
	ax.grid(True)

time = loadtxt('output/time.txt')

with writer.saving(fig,"output/movie.mp4", 100):
	for isnap in range(1,itn,1):
		dd = loadtxt('output/position'+str(isnap)+'.txt')
		xx = dd[:,0]
		yy = dd[:,1]
		zz = dd[:,2]
		MakePlot3D(ax1,xx,yy,zz,isnap)
		ax1.hold(False)
		# MakePlot2D(ax2,xx,yy,isnap)
		# ax2.hold(False)
		writer.grab_frame()
		print('plot = ' + str(isnap)+ 'Done')

# with writer.saving(fig, "output/movie.mp4", 100):	
#     for isnap in range(1,itn,1):
#     	# start = time.time()
#     	dd = loadtxt('output/position'+str(isnap)+'.txt')
#     	zz = dd[:,2]
#     	yy = dd[:,1]
#     	xx = dd[:,0]
#     	# plt.text(0.5, 0.5, str(time[isnap]))
#     	# l.set_data(xx,yy,zz)
#     	ax.plot(xx,yy,zz,'o-')
#     	ax.set_xlim(-1, 1)
#     	ax.set_ylim(-0.5, 0.5)
#     	ax.set_zlim(0,1.3)
#     	plt.title(str(time[isnap]))
#     	# print(time[isnap])
#     	writer.grab_frame()
#     	ax.hold(False)
#     	# end = time.time()
#     	# print(end-start)
#     	print('plot = ' + str(isnap)+ ' Done' )

