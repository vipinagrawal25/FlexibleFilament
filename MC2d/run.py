import os as os
# import pencil_old as pc
import numpy as np
import matplotlib.pyplot as P
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
#%matplotlib notebook
from importlib import reload
import h5py
import yaml
#from stripack import trmesh
import FuncMC2d as F
import dump_visit as dv
from FuncMC2d import MESH
import math
import meshzoo
#
# P.style.use('matplotlibrc')
#
## moved all the styling to matplotlibrc
#---------------------------------------------#
with open(r'input.yaml') as file:
    """input for the run"""
    inp = yaml.full_load(file)
#F.triang_moebius()
#F.triang_sph()
Np=inp['num_particles']
#
debug = inp['debug']
if debug:
    rrini = np.loadtxt('final_rrini.dat')
    # np.savetxt('initial_rrini.dat', rrini)
else:
    if (inp['regular_lattice']):
        points,cells=meshzoo.icosa_sphere(48)
        # print("Intersection = ",F.mesh_intersect(points,cells))
        Np = points.shape[0]
        print("Np =",Np)
    else:
        rrini = F.rand_sph(Np)
        np.savetxt('initial_rrini.dat', rrini)
#
if(inp['do_montecarlo']):
    rr=rrini
    rad=1
    lopt=np.sqrt(8*np.pi*rad**2/(2*Np-4))
    sigma=lopt/(2**(1/6))
    print("sigma=",sigma)
    sv=F.SphVoronoi(rr)
    dv.dump_visit('output/rrini'+str(0).zfill(4)+'.vtk', sv.points, sv._simplices)
    with open('energy.dat', "w") as file:
        for i in range(1,1000):
            rr,E = F.MC_surf(rr,Np,Lone=np.pi,Ltwo=2*np.pi,metric='sph',maxiter=10**(int(math.log10(Np))+1),kBT=1.,
                                 dfac=Np,interactive=False,sigma=sigma)
            sv=F.SphVoronoi(rr)
            obtuse = F.check_obtuse(sv.points, sv._simplices)
            mesh=MESH(BB=2.5,HH=1,sv=sv)
            print("|avlij0-lopt|/lopt =",np.abs(np.mean(mesh.lij0)-lopt)/lopt, obtuse.sum())
            dv.dump_visit('output/rrini'+str(i).zfill(4)+'.vtk', sv.points, sv._simplices)
            np.savetxt('output/rrini'+str(i).zfill(4)+'.dat', rr)
            file.write("%d %15.8f %15.8f\n" %(i, E, np.sum(obtuse)))
            file.flush()
    points=sv.points
    cells=sv.cells
else:
    if inp['regular_lattice']:
        pass
    else:    
        rr=rrini
        sv = F.SphVoronoi(rr)
        points=sv.points
        cells=sv._simplices
#
# if(inp['read_ini_particle']):
#     hf=h5py.File("fin_pos.h5","r")
#     rr=np.array(hf.get('rr'))
#     hf.close()
# sv = F.SphVoronoi(rr)
# cmlst,node_neighbour,bond_neighbour=F.neighbours(sv)
#
eps=0.05
mesh=MESH(BB=50*eps,HH=1,R=points,cells=cells)
avlij0 = np.mean(mesh.lij0)
print(avlij0)
exit()
YY=inp['facH']*eps/(avlij0*avlij0)
mesh.HH=YY*np.sqrt(3)/2
print("mesh.HH =",mesh.HH)
F.WriteMesh(points,mesh.cmlst,mesh.node_nbr,mesh.bond_nbr,cells,
          fname='../monte_carlo/input/input_icosa')
#
with open('energy.dat', "w") as file:
    dv.dump_visit('output/var'+str(0).zfill(4)+'.vtk', mesh.R, cells)
    obtuse = F.check_obtuse(mesh.R, cells)
    dv.dump_visit_cells_scalar('output/var'+str(0).zfill(4)+'.vtk', cells, obtuse, 
            name_scalar='is90')
    file.write("%d %15.8f %15.8f\n" %(0, 0.0, np.sum(obtuse)))
    file.flush()
    for i in range(1,10**(int(math.log10(Np)))):
        mesh,E=F.MC_mesh(mesh,kBT=1.,interactive=True,dfac=32,rr=lopt)
        obtuse = F.check_obtuse(mesh.R, cells)
        dv.dump_visit('output/var'+str(i).zfill(4)+'.vtk', mesh.R, cells)
        dv.dump_visit_cells_scalar('output/var'+str(i).zfill(4)+'.vtk', cells, obtuse, 
            name_scalar='is90')
        file.write("%d %15.8f %15.8f\n" %(i, E, np.sum(obtuse)))
        file.flush()