import sys
sys.path.append("../tools/")
from MESH import *
import functiontools as ft
import numpy as np
import meshzoo
# ############################# FUNCTIONS ###################################
# def SphVoronoi(rr,R=1,lplot=False):
#     Np = np.shape(rr)[0]
#     xyz = np.zeros([Np,3])
#     for ip in range(Np):
#         tht = rr[ip,0]
#         phi = rr[ip,1]
#         x = R*np.sin(tht)*np.cos(phi)
#         y = R*np.sin(tht)*np.sin(phi)
#         z = R*np.cos(tht)
#         xyz[ip] = np.array([x,y,z])
#     sv = SphericalVoronoi(xyz,radius=R)
#     if lplot:
#         plot_voronoi(xyz,sv)
#     return sv
# # ------------------------ save datasets ------------------------------------ #
# def dump_mesh_hdf5(msh,file):
#     if file.split(".")[-1]=="h5":
#         pass
#     else:
#         file=file+".h5"
#     hf = h5py.File(file,'w')
#     hf.create_dataset('pos',data=msh.points)
#     hf.create_dataset('triangles',data=msh.cells)
#     hf.create_dataset('cumu_list',data=msh.cmlst)
#     hf.create_dataset('node_nbr',data=msh.node_nbr)
#     nbn = np.zeros([len(msh.node_nbr),2], dtype=int)
#     for i,bn in enumerate(msh.bond_nbr):
#         nbn[i,0] = bn[0]
#         nbn[i,1] = bn[1]
#     hf.create_dataset('bond_nbr',data=nbn)
#     hf.close()
# ###################### read rr #############################
# filei = sys.argv[1]
# if len(sys.argv)==2:
#     fileo = filei.split("/")[-1].split(".")[0]
# rr = np.loadtxt(filei)
# sv = ft.SphVoronoi(rr)
# # cumlst, node_neighbour, bond_neighbour = ft.neighbours(sv.points,sv._simplices)
# msh = Mesh(sv.points,sv._simplices)
# msh.assign_nbrs()
# msh.dump_hdf5(fileo)
# ###################### Meshzoo ###############################
points,cells=meshzoo.icosa_sphere(65)
fileo="icosa"+"_N"+str(points.shape[0])
# print(fileo)
# exit()
msh = Mesh(points,cells)
msh.assign_nbrs()
msh.dump_hdf5(fileo)