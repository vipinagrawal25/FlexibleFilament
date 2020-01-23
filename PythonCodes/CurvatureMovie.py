from numpy import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=5, metadata=metadata)

FILE = 'output'
dd = loadtxt(FILE+'/curvature.txt')
nrowcol=dd.shape
nsnap=nrowcol[0]
Np = nrowcol[1]
Np=Np-1
dd_SS = loadtxt(FILE+'/material_point.txt')
time = dd[:,0]

fig=plt.figure()
with writer.saving(fig, FILE+"/movieCurvature.mp4", 100):
	for isnap in range(0,nsnap-1,1):
		curvature=dd[isnap,1:Np+1]*(1/(Np*Np))
		SS = dd_SS[isnap,:]/dd_SS[isnap,-1]
		ax=fig.add_subplot(1,1,1)
		ax.plot(SS,curvature)
		plt.xlabel('arc length')
		plt.ylabel('curvature')
		ax.grid(True)
		ax.set_ylim(0,4)
		ax.set_xlim(0,1)
		# ax.set_yscale('log')
		# double const aa =  1.0/(double)(Np); 	// distance between two nodes.
		writer.grab_frame()
		ax.clear()
		print('plot = ' + str(isnap)+ 'Done')

