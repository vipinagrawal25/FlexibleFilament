from numpy import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# from pylab import *
# import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation
# FFMpegWriter = manimation.writers['ffmpeg']
# metadata = dict(title='Movie Test', artist='Matplotlib',
#                 comment='Movie support!')
# writer = FFMpegWriter(fps=5, metadata=metadata)
FILE = 'output'
itn = 2000
Np = 100
dd = loadtxt(FILE+'/curvature.txt')
dd_SS = loadtxt(FILE+'/material_point.txt')
time = dd[:,0]
fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
for isnap in range(1600,2000):
	plt.plot(dd_SS[isn
		ap],dd[isnap,1:Np+1])
	# ax.hold(True)
plt.savefig('temp.eps')
plt.show()
# with writer.saving(fig, FILE+"/movieCurvature.mp4", 100):
# 	for isnap in range(1,itn-1,1):
# 		curvature = log10(dd[isnap,1:Np+1]*(1/Np))
# 		SS = dd_SS[isnap]/1.21
# 		ax.plot(SS,curvature)
# 		plt.xlabel('arc length')
# 		plt.ylabel('curvature')
# 		ax.hold(False)
# 		ax.grid(True)
# 		ax.set_ylim(-13,0)
# 		ax.set_xlim(0,1)
# 		# ax.set_yscale('log')
# 		# double const aa =  1.0/(double)(Np); 	// distance between two nodes.
# 		writer.grab_frame()
# 		print('plot = ' + str(isnap)+ 'Done')