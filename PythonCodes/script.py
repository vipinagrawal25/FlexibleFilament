from functions import *
import matplotlib.pyplot as plt

frequency = 1.5;
Np = 101;
TMAX= 80;
Totalfiles=4000;
file = 'MSD.txt';


curvatureplot(omega=frequency, Np=Np)
# figure=MSD_plot(step=1)
# # show(figure)
# plt.savefig('MSD_complete.png', dpi=400)
# close()

# chitra = plt.figure()
# chitra=MSD_plot(fig=chitra,step=Totalfiles*3.14*2/(frequency*TMAX))
# step=Totalfiles*3.14*2/(frequency*TMAX)
# print(step)
# show(figure)
	# plt.savefig('MSD_cycle.png', dpi=400)
	# close()
	