from functions import *

frequency = 0.6;
Np = 100;
TMAX= 200*frequency/3.14;

curvatureplot(omega=frequency)
MSD_plot(step=2000/TMAX)
