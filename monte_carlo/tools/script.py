from functiontools import MESH
import os
import visit_IO as vio
import sys
import functiontools as f
#################FUNCTIONS################
def foldername(file):
	'''The function return the folder name for a given filename'''
	x=file.split("/")
	name=""
	for n in x[0:-1]:
		name=name+n
	return name
#################SCRIPT################
# if len(sys.argv < 2):
# 	print("Please pass the file name for diagnostics.")
# 	exit(1)
# What variables do you want to compute?
argc=len(sys.argv)
compute={"obtuse":0,"curvature":0}
fnames=[]
for arg in sys.argv[1:]:
	if arg in compute.keys():
		compute[arg]=1;
	else:
		fnames.append(arg)
for file in fnames:
	folder=foldername(file)
	Np=f.read_param(folder+'/para_file.in')['N']
	points,cells=vio.vtk_to_data(infile=file,Np=Np)
	mesh = MESH(points,cells)
	foutname="something1.vtk"
	if compute['obtuse']==1:
		isobtuse=mesh.compute_obtuse()
		# foutname+="_obtuse"
	if compute['curvature']==1:
		mesh.assign_nbrs()
		curv=mesh.curvature()
		vio.dump_visit(foutname, points, cells)
		vio.dump_visit_points_scalar(foutname, points, curv, name_scalar='curvature')
		# foutname+="_curv"
	# if compute['obtuse']==1:
	# 	vio.dump_visit_points_scalar(foutname, points, isobtuse, name_scalar='obtuse_data')
	# if compute['curvature']==1:
	# 	vio.dump_visit_points_scalar(foutname, points, curv, name_scalar='curvature')