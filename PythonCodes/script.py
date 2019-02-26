from functions import *
import matplotlib.pyplot as plt

frequency = 3;
Np = 100;
TMAX= 15;

curvatureplot(omega=frequency)

figure=MSD_plot(step=1)
# show(figure)
plt.savefig('MSD_complete.png', dpi=400)
close()

figure=MSD_plot(step=2000*3.14*2/(frequency*TMAX))
# show(figure)
plt.savefig('MSD_cycle.png', dpi=400)
close()
