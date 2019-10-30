from functions import *
import matplotlib.pyplot as plt
import moviewriter as MOVIE

frequency = 3;
Np = 100;
TMAX= 15;
Totalfiles=2000;

#This is an array of things that you want to calculate.
Thing1 = "Movie"
Thing2 = "MSD"

Things = [Thing2]

curvatureplot(omega=frequency, Np=Np)
figure=MSD_plot(step=1)
show(figure)
plt.savefig('MSD_complete.png', dpi=400)
close()

for i in range(0,len(Things)):
	if Things[i] == "movie" or Things[i] == "Movie":
		MOVIE.MultifileMovie(nsnap=2000)
	elif Things[i] == "MSD":
		chitra = plt.figure()
		step=Totalfiles*3.14*2/(frequency*TMAX)
		chitra=MSD_plot(fig=chitra,step=step)
		plt.savefig('MSD_cycle.png', dpi=400)
		close()


# BendE = BendingEnergy(AA=0.00006,FILE='data/curvature.txt')
# plt.plot(BendE)
# plt.savefig('data/BendE.png',dpi=300)

# StretchE = StretchEnergy(HH=39.4,FILE='data/material_point.txt')
# plt.figure()
# plt.plot(StretchE)
# plt.savefig('data/StretchE.png',dpi=300)

# TotalE=BendE+StretchE
# plt.figure()
# plt.plot(TotalE)
# plt.savefig('data/TotalE.png',dpi=300)

