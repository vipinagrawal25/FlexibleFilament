from pylab import *
from scipy import ndimage
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D

# This section is intended for calculating the surface plot of curvature.
Np = 20

FILE = 'output'

dd = loadtxt(FILE+'/curvature.txt')
dd_SS = loadtxt(FILE+'/material_point.txt')
dd = dd/400
Y = dd_SS
timeAxis = loadtxt(FILE+'/time.txt')/3.14


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
# ax.set_zlabel('Curvature')
# ax.set_zlabel('Curvature')
show()
