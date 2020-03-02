from functions import *
import matplotlib.pyplot as plt
import moviewriter as MOVIE

AA = 1.5*0.00001*1
frequency = 0;

HH = 64*AA*(40**2)
Np = 100;
aspect_ratio=1
# TMAX= 15;
# Totalfiles=2000;

#This is an array of things that you want to calculate.
Thing1 = "Movie"
Thing2 = "MSD"
Thing3 = "Curvature"
Thing4 = "Energy"
Thing5 = "LeebyL"
Thing6 = "SineCurvature"
Thing7 = "MSD_no_trans"

Things = [Thing1]

# figure=MSD_plot(step=1)
# show(figure)
# plt.savefig('MSD_complete.png', dpi=400)
# close()
for i in range(0,len(Things)):
	if Things[i] == "movie" or Things[i] == "Movie":
		MOVIE.MultifileMovie(dim=2)
	elif Things[i] == "MSD":
		MSD_plot()
		close()
	elif Things[i] == "Curvature" or Things[i] == "curvatureplot" or Things[i] == "curvature":
		curvatureplot(omega=frequency)
	elif Things[i] == "Energy":
		time,BE,SE = Energy(AA,HH)
		plt.plot(time[0:-1],BE[0:-1])
		peaks = np.where((BE[1:-1] > BE[0:-2]) * (BE[1:-1] > BE[2:]))[0] + 1
		plt.plot(time[peaks],BE[peaks],'o-')
		plt.savefig('Energy.eps')
		close()
	elif Things[i] == "LeebyL":
		LeebyL(Np=100)
	elif Things[i]=="SineCurvature" :
		SineCurvature(aspect_ratio=aspect_ratio)
		close()
	elif Things[i]=="MSD_no_trans":
		MSD_no_trans()
# BendE = BendingEnergy(A=0.00006,FILE='data/curvature.txt')
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