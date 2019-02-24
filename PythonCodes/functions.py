from pylab import *
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def curvatureplot( FILE='output', Np=100, omega=3 ):
	dd = loadtxt(FILE+'/curvature.txt')
	dd_SS = loadtxt(FILE+'/material_point.txt')
	dd = dd
	Y = dd_SS
	timeAxis = loadtxt(FILE+'/time.txt')*omega/3.14
	IndexAxis = arange(Np)
	heightAxis = dd[:,1:Np+1]
	# heightAxis = ndimage.rotate(heightAxis,90)
	# imshow(heightAxis,cmap=plt.cm.gray)
	# colorbar()

	(IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
	fig = figure()
	ax = fig.add_subplot(111)
	surf = ax.contourf(X,Y,heightAxis)
	cbar = colorbar(surf)
	ax.set_xlabel('time')
	ax.set_ylabel('material_point')

	plt.savefig('curvatureplot.eps')
	plt.savefig('curvatureplot.png',dpi=300)
	# show()
	close()
	return 

def MSD_plot(FILE='.',step=1):
	step = int(step)
	print (step)
	MSD = loadtxt(FILE+'/MSD.txt',delimiter=';',usecols=range(2000))
	# MSD = loadtxt(FILE+'/MSD.txt')
	fig = figure()
	fig=plt.plot(MSD[1::step],'o')
	return fig


