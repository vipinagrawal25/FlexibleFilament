from pylab import *
from scipy import ndimage
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D

# This section is intended for calculating the surface plot of curvature.
Np = 20
dd = loadtxt('10/curvature.txt')
dd_SS = loadtxt('10/material_point.txt')
dd = dd
Y = dd_SS[1500,2000,:]
timeAxis = 0.02*dd[1500:2000,0]


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
show()