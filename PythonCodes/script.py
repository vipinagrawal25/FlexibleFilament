from functions import *
import matplotlib.pyplot as plt
import moviewriter as MOVIE
import numpy as NP
import sys
print(sys.argv[1])
# exit(1)
#
AA = 1.5*0.00001*1*0.4**4
Np = 41
HH = 16*AA*(0.005**2)
aspect_ratio=1
length=1.28
sigma=0.75
ShearRate=2
omega = ShearRate*sigma

# This is an array of things that you want to calculate.
Thing1 = "Movie"
Thing2 = "MSD"
Thing3 = "Curvature"
Thing4 = "Energy"
Thing5 = "LeebyL"
Thing6 = "SineCurvature"
Thing7 = "MSD_no_trans"
Thing8 = "SqrtCurvature"
Thing9 = "TracerVelocity"

Things = [Thing8]
Xtracer2 = [0.01,0,0.5]
Xtracer9 = [0.01,1.2,0.5]

for i in range(0,len(Things)):
	if Things[i] == "movie" or Things[i] == "Movie":
		MOVIE.MultiFileMovie(FILE='output',dim=2)
	elif Things[i] == "MSD":
		MSD_plot()
		close()
	elif Things[i] == "Curvature" or Things[i] == "curvatureplot" or Things[i] == "curvature":
		curvatureplot(omega=omega,length=length)
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
	elif Things[i]=="SqrtCurvature":
		GetCurv(Folder=sys.argv[1],dim=2,wDataMeth=1)
	elif Things[i]=="TracerVelocity":
		Vtracer=VtraceTimeSer(sigma=0.75,Xtracer=Xtracer9)
		NP.save('output/Vtracer9.npy',Vtracer)