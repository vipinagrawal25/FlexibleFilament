import numpy as np
import matplotlib.pyplot as plt

Np = 256
iext = int(np.sqrt(Np))
jext = int(np.sqrt(Np))
npoly = (iext-1)*(iext-1);
extent = np.asarray([-0.1, 0.1, -0.1, 0.1])
dx = (extent[1] - extent[0])/jext;
rad = 0.2;
z_bot = 1.05
with open("afm_tip.vtk", "w") as f:
   f.write('# vtk DataFile Version 2.0 \n')
   f.write('grid, time 110\n')
   f.write('ASCII \n')
   f.write('DATASET POLYDATA \n')
   f.write('POINTS  '+str(Np)+'  float\n')
   for j in range (0, jext):
       for i in range (0, iext):
            x0 = extent[0] + dx*(i + 0.5);
            y0 = extent[2] + dx*(j + 0.5);
            z0 = (x0/rad)**2 + (y0/rad)**2 + z_bot; 
            f.write("%16.8f %16.8f %16.8f\n" %(x0, y0, z0))

   f.write('POLYGONS  '+str(npoly)+'  ' +str(5*npoly) + '\n')
   for j in range (0, jext-1):
       for i in range (0, iext-1):
           here = j*iext + i;
           est = (i+1)%iext + ((iext + j)%iext) *iext;
           nth = ((iext + i)%iext) + (j + 1 + iext)%iext *iext;
           n_est = ((i + 1)%iext) + ((j + 1 + iext)%iext)*iext;
           # print(here, est, nth, n_est)
           f.write("%d %d %d %d %d\n" %(4, here, est, n_est, nth))


x = np.asarray([-1.2,1.2])
y = np.asarray([-1.2,1.2])
z = np.asarray([-1.05,-0.85])

with open("bot_wall.vtk", "w") as f:
    f.write('# vtk DataFile Version 2.0 \n')
    f.write('grid, time 110\n')
    f.write('ASCII \n')
    f.write('DATASET POLYDATA \n')
    f.write('POINTS  '+str(8)+'  float\n')
    f.write("%16.8f %16.8f %16.8f\n" %(x[0], y[0], z[0]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[1], y[0], z[0]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[1], y[1], z[0]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[0], y[1], z[0]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[0], y[0], z[1]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[1], y[0], z[1]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[1], y[1], z[1]))
    f.write("%16.8f %16.8f %16.8f\n" %(x[0], y[1], z[1]))
    f.write('POLYGONS  '+str(6)+'  30\n')
    f.write('4 0 1 2 3\n')
    f.write('4 4 5 6 7\n')
    f.write('4 0 1 5 4\n')
    f.write('4 2 3 7 6\n')
    f.write('4 0 4 7 3\n')
    f.write('4 1 2 6 5\n')


