from MESH import *
import os
import visitio as vio
import sys
import functiontools as f
import paramio as pio
import meshzoo
################# SCRIPT ################
# if len(sys.argv < 2):
# 	print("Please pass the file name for diagnostics.")
# 	exit(1)
# What variables do you want to compute?
argc=len(sys.argv)
compute={"obtuse":0,"curvature":0,
		"free_energy":0,"freeenergy":0,"free-energy":0}
fnames=[]
for arg in sys.argv[1:]:
	if arg in compute.keys():
		compute[arg]=1;
	else:
		fnames.append(arg)
# ---------------------------- Now call required functions -------------------------------- #
for file in fnames:
	folder=f.foldername(file)
	Np=pio.read_param(folder+'/para_file.in')['N']
	points,cells=vio.vtk_to_data(infile=file,Np=Np)
	mesh = MESH(points,cells)
	foutname=''.join(file.split(".")[0:-1])
	if compute['obtuse']==1:
		isobtuse=mesh.compute_obtuse()
		foutname+="_obtuse"
	if compute['curvature']==1:
		mesh.assign_nbrs()
		curv=mesh.curvature()
		foutname=foutname+"_curv"
	if compute['free_energy'] ==1 or compute['freeenergy'] ==1 or compute['free-energy']==1:
		tot_ener=np.loadtxt(folder+"/mc_log")[:,2]
		FF=f.free_energy(tot_ener)
		plt.plot(FF)
	vio.dump_visit(foutname+".vtk", points, cells)
	vio.dump_visit_points_scalar(foutname+".vtk", points, curv, name_scalar='curvature')