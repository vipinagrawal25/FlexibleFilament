import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.special import legendre
from scipy.special import sph_harm
# from vtk import *
# from vtk.util.numpy_support import vtk_to_numpy
############# PARAMETERS #########################
Np=23042
############# FUNCTIONS #########################
def vtk_to_data(infile,Np=Np,skiprows=5):
	# Read data from vtk file.
	points=np.loadtxt(infile,skiprows=skiprows,max_rows=Np)
	cells=np.loadtxt(infile,skiprows=skiprows+Np+1,max_rows=2*Np-4)
	return points, cells
# ---------------------------------------#
def 
############# SCRIPT #########################
len_arg = len(sys.argv)
if len_arg>1:
	folder = sys.argv[1]
	filen = sys.argv[2].zfill(5)
else:
	folder="../output_fig2_yellow/"
	filen = "02990"
infile = folder+"part_"+filen+".vtk"
infile0 = folder+"part_00000.vtk"
# ----------get the data ----------------#
points,cells = vtk_to_data(infile)
points0,cells0 = vtk_to_data(infile0)
# ---------------------------------------#
height_field = points-points0