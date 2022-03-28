# from MESH import *
import os
import visitio as vio
import sys
import functiontools as f
import paramio as pio
import matplotlib.pyplot as plt
import numpy as np
from MESH import *
# import matplotlib
# matplotlib.use('TkAgg')
################# SCRIPT ################
# if len(sys.argv < 2):
# 	print("Please pass the file name for diagnostics.")
# 	exit(1)
# What variables do you want to compute?
argc=len(sys.argv)
compute={"obtuse":0,"curvature":0,
		"free_energy":0,"freeenergy":0,"free-energy":0,
		"mean_energy":0,"meanE":0,
		"pdf":0
		}
fnames=[]
for arg in sys.argv[1:]:
	if arg in compute.keys():
		compute[arg]=1;
	else:
		fnames.append(arg)
# foutname = "0.75_curv_fixed"
# ------------------ Read data -------------------- #
# nrun = len(fnames)
# mc_log=np.empty(nrun,dtype=object)
# tz=np.zeros(nrun)
# for ifol,folder in enumerate(fnames):
# 	mc_log[ifol]=np.loadtxt(folder+"/mc_log")
# 	tz[ifol]=pio.read_param(fname=folder+'/para_file.in')['tip_pos_z']
# N=mc_log[0].shape[0]
# ------------------ Free energy -------------------- #
if compute['free_energy'] ==1 or compute['freeenergy'] ==1 or compute['free-energy']==1:
	FFall = np.zeros(nrun)
	for ifol in range(nrun):
		tot_ener=mc_log[ifol][:,2]
		Eminis,FF=f.free_energy(tot_ener)
		FFall[ifol]=FF[-1]
		plt.plot(FF,'.-')
		print(FF[-1])
	plt.figure()
	plt.plot(FFall,'.-')
# ------------------ Mean energy -------------------- #
	# plt.grid(True)
# ------------------ PDF -------------------- #
if compute['pdf']==1:
	for ifol in range(nrun):
		plt.figure()
		data=mc_log[ifol][:,2]
		hist=np.histogram(data,bins=25,density=True)
		plt.plot(hist[1][1:],hist[0],'.-')
plt.show()
for file in fnames:
	if compute['obtuse']==1 or compute['curvature']==1:
		folder=f.foldername(file)
		Np=pio.read_param(folder+'/para_file.in')['N']
		points,cells=vio.vtk_to_data(infile=file,Np=Np)
		mesh = Mesh(points,cells)
		# foutname=file
		foutname=''.join(file.split(".")[0:-1])
	if compute['obtuse']==1:
		isobtuse=mesh.compute_obtuse()
		foutname+="_obtuse"
	if compute['curvature']==1:
		mesh.assign_nbrs()
		curv=mesh.curvature()
		foutname=foutname+"_curv"
		print(curv)
	vio.dump_visit(foutname+".vtk", points, cells)
	vio.dump_visit_points_scalar(foutname+".vtk", points, curv, name_scalar='curvature')