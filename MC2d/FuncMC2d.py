from __future__ import division
import os as os
import numpy as np
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
from scipy.interpolate import interp1d
from numpy import linalg as LA
import h5py
import matplotlib.cm as cm
from scipy.spatial import Delaunay
from scipy.spatial import SphericalVoronoi, geometric_slerp
from mpl_toolkits.mplot3d import proj3d
from dataclasses import dataclass
import math
import ttisect as ttis
#-----------------------------------------
class MESH:
    def __init__(self,BB,HH,R,cells):
            # if (R is None) or (cmlst is None) or (node_nbr is None) or (bond_nbr is None):
            #     cmlst,node_nbr,bond_nbr=neighbours(sv.points,sv._simplices)
            #     R=sv.points
            cmlst,node_nbr,bond_nbr=neighbours(R,cells)
            self.BB=BB
            self.HH=HH
            self.Np=R.shape[0]
            self.R=R
            self.cmlst=np.array(cmlst).astype(np.int)
            self.node_nbr=np.array(node_nbr).astype(np.int)
            self.bond_nbr = np.array(bond_nbr).astype(tuple)
            self.lij0=self.lengths()
            self.radius=1
            self.sp_curv=2/self.radius
#
    def lengths(self):
        R=self.R
        Np=self.Np
        lengs=np.zeros(self.cmlst[-1],dtype=np.double)
        for i in range(Np):
            start=self.cmlst[i]
            end=self.cmlst[i+1]
            count=0
            for j in self.node_nbr[start:end]:
                lengs[start+count]=LA.norm(R[i]-R[j])
                count=count+1
        return lengs
    #
    def cot(self,i,j,k):
        #TODO: Compute cot angle by actual size of bonds.
        # angle theta for triangle ijk
        R=self.R
        rik=R[i]-R[k]
        rjk=R[j]-R[k]
        # dotP=np.inner(rik,rjk)/(LA.norm(rjk)*LA.norm(rik))
        # return dotP/np.sqrt(1-dotP*dotP)
        return np.inner(rik,rjk)/LA.norm(np.cross(rik,rjk))
    #
    def bend_energy(self,i):
        BB=self.BB
        R=self.R
        start=self.cmlst[i]
        end=self.cmlst[i+1]
        #
        sigma_i=0
        cot_times_rij=np.zeros(3)
        count=0
        for j in self.node_nbr[start:end]:
            rij=R[i]-R[j]
            k=self.bond_nbr[start+count][0]
            kprime=self.bond_nbr[start+count][1]
            cot_sum=0.5*(self.cot(i,j,k)+self.cot(i,j,kprime))
            sigma_i=sigma_i+np.inner(rij,rij)*cot_sum
            cot_times_rij=cot_times_rij+cot_sum*rij
            #
            count=count+1
        sigma_i=sigma_i/4
        curv=(1/sigma_i)*LA.norm(cot_times_rij)
        curv0=self.sp_curv
        return 0.5*BB*sigma_i*(curv-curv0)**2
    #
    def stretch_energy(self,i):
        HH=self.HH
        R=self.R
        lij0=self.lij0
        start=self.cmlst[i]
        end=self.cmlst[i+1]
        # print(type(self.node_nbr[start:end]))
        count=0
        stretchE=0
        for j in self.node_nbr[start:end]:
            rij=R[i]-R[j]
            stretchE=stretchE+(LA.norm(rij)-lij0[start+count])**2
            count=count+1
        return 0.5*HH*stretchE
    #
    def bend_energy_nbr(self,i):
        BB=self.BB
        R=self.R
        start=self.cmlst[i]
        end=self.cmlst[i+1]
        #
        bendE=self.bend_energy(i)
        for j in self.node_nbr[start:end]:
            bendE=bendE+self.bend_energy(j)
        return bendE
    #
    def tot_energy(self):
        beE=0
        stE=0
        for i in range(self.Np):
            beE=beE+self.bend_energy(i)
            stE=stE++0.5*self.stretch_energy(i)
        # print("beE=",beE)
        # print("stE=",stE)
        return beE+stE
#-----------------------------------------
def mesh_intersect(points,cells):
    """
    Right now I am writing O(n^2) algorithm.
    points=[N,3] as numpy double array
    cells=[2N-4,3] as numpy integer array
    """
    N=points.shape[0]
    Ntr=2*N-4
    Nisect=0
    for i in range(Ntr):
        for j in range(i+1,Ntr):
            Nisect+=ttis.tri_tri_isect(points[cells[i,0]],points[cells[i,1]],points[cells[i,2]],
                                       points[cells[j,0]],points[cells[j,1]],points[cells[j,2]])
    return Nisect
#-----------------------------------------
def MC_step_mesh(mesh,maxiter=100,kBT=1.,dfac=64):
    Np=mesh.Np
    move = 0
    for iter in range(maxiter):
        """ Choose a random integer between 0 and N-1 """
        kp = np.random.randint(low=0,high=Np)
        Eini,do_mc=denergy_kp(mesh,kp,lwall=False)
        # Eini = En(rr,kp,metric=metric)
        Rrem = mesh.R[kp].copy()
        rand_increment_mesh(mesh,kp)
        Efin,do_mc=denergy_kp(mesh,kp,lwall=False)
        if do_mc:
            dE = Efin-Eini
            Accept = Metropolis(dE,kBT)
        else:
            Accept=False
        if Accept:
            move = move+1
        else:
            mesh.R[kp]=Rrem.copy()
    E = mesh.tot_energy()
    # print(rr[0,0],rr[1,0])
    return move,E
#------------------------------
def check_obtuse(points, triangles):
    """
    Takes the triangles and points as input and returns 
    an array of shape number of triangles 
    #
    The output is 1 if one of the triangle is greater than
    90, and 0 if all the angle is acute (less than 90)
    """
    # angle theta for triangle ijk
    isgt90 = np.zeros(np.shape(triangles)[0], 
            dtype=np.float64)
    piby2 = 0.5*np.pi
    a_ij_ik=np.zeros(triangles.shape[0])
    a_ji_jk=np.zeros(triangles.shape[0])
    a_ki_kj=np.zeros(triangles.shape[0])
    for i, tris in enumerate(triangles):
        ri = points[tris[0]]
        rj = points[tris[1]]
        rk = points[tris[2]]
        rij = ri - rj
        rij = rij/LA.norm(rij)
        #
        rik = ri - rk
        rik = rik/LA.norm(rik)
        #
        rkj = rk - rj
        rkj = rkj/LA.norm(rkj)
        #
        a_ij_ik[i] = np.arccos(np.inner(rij,rik))
        a_ji_jk[i] = np.arccos(np.inner(-rij,-rkj))
        a_ki_kj[i] = np.arccos(np.inner(-rik,rkj))
        #
        logic = ((a_ij_ik[i] > piby2) or 
                (a_ji_jk[i] > piby2) or 
                (a_ki_kj[i] > piby2))
        if(logic):
            isgt90[i] = 1e0
        # print (a_ij_ik, a_ji_jk, a_ki_kj)
        # print (a_ij_ik + a_ji_jk + a_ki_kj)
    return isgt90
#-----------------------------------------
def denergy_kp(mesh,kp,lwall=True,zwall=-1-1/16):
    eval_other_energies=True
    dE=0.0
    LJener = 0.0
    if lwall:
        zz=mesh.R[kp,2]
        dz = zz - zwall
        if(dz < 0e0):
            # I can go ahead with evaluating the lj energy
            LJener = LJ(dz)
        else:
            eval_other_energies = False
    if(eval_other_energies):
        dE = (mesh.bend_energy_nbr(kp)
            +mesh.stretch_energy(kp))
    return dE+LJener,eval_other_energies
#-----------------------------------------
def LJ(zz,sigma=1/32,zcut=None,epsilon=4):
    if zcut is None:
        zcut=4*sigma
    if zz>zcut:
        pot=0
    else:
        sbyz6 = (sigma/zz)**6
        pot = 4*epsilon*(sbyz6**2-sbyz6)
    return pot
#-----------------------------------------
def rand_increment_mesh(mesh,kp,rr=1,dfac=64):
    mesh.R[kp] = mesh.R[kp] + (rr/dfac)*np.random.uniform(low=-1,high=1,size=3)
#-----------------------------------------
def MC_mesh(mesh,maxiter=None,kBT=1.,interactive=True,dfac=64):
    Np=mesh.Np
    if maxiter is None:
        maxiter=10**(int(math.log10(Np))+1)
    # cont=True
    Eini = mesh.tot_energy()
    print('Eini/Np=' + str(Eini/Np))
    # if interactive:
    #     fig = P.figure()
    #     ax = fig.add_subplot(111)
    #     ax.plot(rr[:,0],rr[:,1],'o',color='C0')
    move,Efin = MC_step_mesh(mesh,maxiter=maxiter,kBT=kBT,dfac=dfac)
    print('MC steps', maxiter)
    print('accepted moves',move)
    print('Efin/Np=',Efin/Np)
        # if interactive:
        #     ax.plot(rr[:,0],rr[:,1],'*',color='C3')
        #     P.show()
        #     answer=input('continue(Y) ? ')
        #     if answer != 'Y':
        #         cont = False
        #         print('Run over, go home')
        #     else:
        #         print('another', maxiter, 'MC steps')
        #         P.close()
        # else:
        #     print('Run over, go home')
        #     break
    return mesh,Efin
#-----------------------------------------
# @dataclass
# class MESH:
#     Np: int
#     R: np.double
#     BB: np.double
#     cmlst: np.int
#     node_nbr: np.int
#     bond_nbr: np.int
#     #
#     def cot(self,i,j,k):
#         # angle theta for triangle ijk
#         R=self.R
#         rik=R[i]-R[k]
#         rjk=R[j]-R[k]
#         return np.inner(rik,rjk)/LA.norm(np.cross(rik,rjk))
#     def dis(self,i,j):
#         return LA.norm(self.R[i]-self.R[j])
#-----------------------------------------
def neighbours(points,cells):
    simpl=sorted_simplices(cells)
    r1=simpl[:,0]
    r2=simpl[:,1]
    r3=simpl[:,2]
    Np=len(points)
    lst=np.zeros(Np,dtype=int)
    cumlst=np.zeros(Np+1,dtype=int)
    for i in range(0, Np):
        lst[i]=len(r1[r1==i])/2
    cumlst[1:] = np.cumsum(lst)
    node_neighbour = np.zeros(cumlst[-1],dtype=int)
    bond_neighbour = np.zeros(cumlst[-1],dtype=tuple)
    for i in range(0, cumlst[-1], 1):
        node_neighbour[i]=r2[2*i]
        bond_neighbour[i]=(r3[2*i],r3[2*i+1])
    return cumlst,node_neighbour,bond_neighbour
#-----------------------------------------#
def sorted_simplices(cells):
    lsimples = len(cells)
    nsimplices = np.asarray([], dtype=np.int32)
    for scles in cells:
        nscles = np.sort(scles)
        nsimplices = np.hstack([nsimplices, nscles])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[2], nscles[0]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[0], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[0], nscles[2], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[0], nscles[2]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[1], nscles[0]]])
    nsimpl = nsimplices.reshape(lsimples*6, 3)
    nsimpl = np.asarray(sorted(nsimpl, key=lambda x: (x[0], x[1])))
    return nsimpl
# #-----------------------------------------#
# def print_neighbors_(sv):
#     lsimples = len(sv._simplices)
#     nsimplices = np.asarray([], dtype=np.int32)
#     for scles in sv._simplices:
#         nscles = np.sort(scles)
#         nsimplices = np.hstack([nsimplices, nscles])
#         nsimplices = np.hstack([nsimplices, [nscles[1], nscles[2], nscles[0]]])
#         nsimplices = np.hstack([nsimplices, [nscles[2], nscles[0], nscles[1]]])
#         nsimplices = np.hstack([nsimplices, [nscles[0], nscles[2], nscles[1]]])
#         nsimplices = np.hstack([nsimplices, [nscles[1], nscles[0], nscles[2]]])
#         nsimplices = np.hstack([nsimplices, [nscles[2], nscles[1], nscles[0]]])
#     nsimpl = nsimplices.reshape(lsimples*6, 3)
#     nsimpl = np.asarray(sorted(nsimpl, key=lambda x: (x[0], x[1])))
#     r1 = nsimpl[:,0]
#     r2 = nsimpl[:,1]
#     r3 = nsimpl[:,2]
#     print("Num neighbors")
#     for i in range(0, len(sv.points)):
#         print(int(len(r1[r1==i])/2))

#     print("neighbors index")
#     for i in range(0, len(r2), 2):
#         print(r2[i])

#     print("adjecent neighbors index")
#     for i in range(0, len(r3), 1):
#         print(r3[i])
# REQUIREMENTS: sv.tris is assigned and defined if not call make_tris function first
def get_tri(sv,point):
    if 'tris' in dir(sv):
        return sv.tris[sv.cntris[point]:sv.cntris[point+1]]
    else:
        print("# tris(triangles for a given point) is not a member of SphericalVoronoi object. Assign it by calling calc_tris function.")
        print("# Anyway I am doing it for you now.")
        sv.tris,sv.cntris=list_to_arr(sv.regions)
        return sv.tris[sv.cntris[point]:sv.cntris[point+1]]
#-----------------------------------------#
def list_to_arr(lst):
    Np=len(lst)
    nelement = np.zeros(Np,dtype=int);
    for ip in range(Np):
        nelement[ip] = len(lst[ip])
    cnelement=np.zeros(Np+1,dtype=int)
    cnelement[1:]=np.cumsum(nelement)
    arr=np.zeros(cnelement[-1],dtype=int)
    trcount=0
    for ip in range(Np):
        arr[cnelement[ip]:cnelement[ip+1]]=lst[ip]
    return arr,cnelement
#-----------------------------------------#
def calc_tris(sv):
    sces = sv._simplices
    Np = len(sv.points)
    # it stores the number of triangle for a given point.
    ntris = np.zeros(Np,dtype=int);
    for ices in sces:
        for ipnt in ices:
            ntris[ipnt]=ntris[ipnt]+1
    cntris=np.zeros(Np+1,dtype=int)
    cntris[1:]=np.cumsum(ntris)
    tris=np.zeros(cntris[-1],dtype=int)
    ind_tris=np.zeros(Np,dtype=int)
    trcount=0;
    for ices in sces:
        for ipnt in ices:
            tris[cntris[ipnt]+ind_tris[ipnt]]=trcount
            ind_tris[ipnt]=ind_tris[ipnt]+1
        trcount=trcount+1
    sv.tris=tris
    sv.cntris=cntris
#-----------------------------------------#
# the function computes the normal for all the triangles and stores them in sv object.
def calc_trinrmls(sv):
    sces=sv._simplices
    points=sv.points
    trinrmls=np.zeros([len(sces),3])
    trcount=0
    for ices in sces:
        trinrmls[trcount,:]=np.cross(points[ices[0]]-points[ices[1]],points[ices[1]]-points[ices[2]])
        trinrmls[trcount,:]=trinrmls[trcount,:]/LA.norm(trinrmls[trcount,:])
        trcount=trcount+1
    sv.trinrmls=trinrmls
    # return np.cross(sv.)
#-----------------------------------------#
# Given coordinates of a point, what is the normal?
# Ans: it's sum of all the normals neighbouring the point
# REQUIREMENTS: sv.tris,sv.trinrmls is assigned and defined if not call make_tris,calc_trinrmls functions first.
def normal(sv,point):
    tris_point=get_tri(sv,point)
    print(tris_point)
    nrml = np.zeros(3)
    for itri in tris_point:
        nrml = nrml+sv.trinrmls[itri]
    nrml=nrml/LA.norm(nrml)
    return nrml
#-----------------------------------------#
def nearest_neighbour(sv,Np=0):
    if (~Np):
        Np = len(sv.regions)
    n_NNL = [len(region) for region in sv.regions]
    cn_NNL = np.zeros(Np+1,dtype=int)
    cn_NNL[1:] = np.cumsum(n_NNL)
    NNL = np.zeros(cn_NNL[-1],dtype=int)
    for ip in range(Np):
        NNL[cn_NNL[ip]:cn_NNL[ip+1]] = np.array(sv.regions[ip])
    return cn_NNL,NNL
#-----------------------------------------#
def SphVoronoi(rr,R=1,lplot=False):
    Np = np.shape(rr)[0]
    xyz = np.zeros([Np,3])
    for ip in range(Np):
        tht = rr[ip,0]
        phi = rr[ip,1]
        x = R*np.sin(tht)*np.cos(phi)
        y = R*np.sin(tht)*np.sin(phi)
        z = R*np.cos(tht)
        xyz[ip] = np.array([x,y,z])
    sv = SphericalVoronoi(xyz,radius=R)
    if lplot:
        plot_voronoi(xyz,sv)
    return sv
#-----------------------------------------#
def plot_voronoi(points,sv):
    t_vals = np.linspace(0, 1, 2000)
    fig = P.figure()
    ax = fig.add_subplot(111, projection='3d')
    # plot the unit sphere for reference (optional)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='y', alpha=0.05)
    # plot generator points
    ax.scatter(points[2:, 0], points[2:, 1], points[2:, 2], c='r')
    # plot Voronoi vertices
    #ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')
    #ax.axhline(0.5, ls=':')
    ax.scatter(points[0:2, 0], points[0:2, 1], points[0:2, 2], c='k')
    # ax.scatter(sv.points[:, 0],sv.points[:, 1], points[:, 2], c='g')
    for region in sv._simplices:
        n = len(region)
        for i in range(n):
            start = sv.points[region][i]
            end = sv.points[region][(i + 1) % n]
            # print(i)
            # print(start)
            # print(end)
            result = geometric_slerp(start, end, t_vals)
            ax.plot(result[..., 0],
                    result[..., 1],
                    result[..., 2],
                    c='k')
    ax.azim = 10
    ax.elev = 40
    _ = ax.set_xticks([])
    _ = ax.set_yticks([])
    _ = ax.set_zticks([])
    #fig.set_size_inches(4, 4)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    P.show()
#-----------------------------------------#
def MC_surf(rr,N,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',maxiter=100,kBT=1.,
            interactive=True,dfac=64,sigma=None):
    # rrini = inidist(N,Lone,Ltwo,metric=metric)
    # rr = rrini
    # print(rr[0,:])
    # WritePos(rr,fname='ini_pos')
    if sigma is None:
        lopt=np.sqrt(8*np.pi*rad**2/(2*Np-4))
        sigma=lopt/(2**(1/6))
    cont=True
    while cont:
        Eini = tot_energy(rr,metric=metric,sigma=sigma)
        print('Eini/N=',Eini/N)
        if interactive:
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(rr[:,0],rr[:,1],'o',color='C0')
        rr,move,Efin = MC_step(rr,Lone=Lone,Ltwo=Ltwo,metric=metric,
                               maxiter=maxiter,kBT=kBT,dfac=dfac,sigma=sigma)
        print('MC steps',maxiter)
        print('accepted moves',move)
        print('Efin/N=',Efin/N)
        # WritePos(rr,fname='fin_pos')
        if interactive:
            ax.plot(rr[:,0],rr[:,1],'*',color='C3')
            P.show()
            answer=input('continue(Y) ? ')
            if answer != 'Y':
                cont = False
                print('Run over, go home')
            else:
                print('another', maxiter, 'MC steps')
                P.close()
        else:
            print('# Run over, go home')
            break
    return rr,Efin/N
#-----------------------------------------#
def inidist(N,Lone=2*np.pi,Ltwo=np.pi,metric='cart'):
    if metric == 'cart':
        rrini = uni_random(N,Lone,Ltwo)
    elif metric == 'sph':
        rrini = rand_sph(N)
    else:
        print(metric,'not coded, quiting')
        quit()
    return rrini
#-----------------------------------------#
def rand_sph(N):
    ran = np.zeros([N,2],dtype='double')
    # the first two points are always at the two poles
    ran[0,0]=0
    ran[0,1]=0
    ran[1,0]=np.pi
    ran[1,1]=0
    for ip in range(2,N):
        rth = np.arccos(np.random.uniform(low=-1.0,high=1.0))
        if rth==0:
            rth = rth+ np.arccos(np.random.uniform(low=-1.0,high=1.0)) 
            #np.random.uniform(low=0.0,high=np.pi)
        if rth==np.pi:
            rth = rth- np.arccos(np.random.uniform(low=-1.0,high=1.0)) 
            #np.random.uniform(low=0.0,high=np.pi)
        ran[ip,0] = rth
        ran[ip,1] = np.random.uniform(low=0.0,high=2*np.pi)
    return ran
#-----------------------------------------#
def plot_pos(ax,rr,rrini):
    N = np.shape(rr)[0]
    print(rrini[0,:])
    print(rr[0,:])
    ax.plot(rrini[:,0],rrini[:,1],'o',color='C0')
    ax.plot(rr[:,0],rr[:,1],'*',color='C3')
    return ax
#-----------------------------------------#
def MC_step(rr,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',maxiter=100,kBT=1.,
            dfac=64,sigma=None):
    N = np.shape(rr)[0]
    move = 0
    for iter in range(maxiter):
        """ Choose a random integer between 0 and N-1 """
        kp = np.random.randint(low=0,high=N)
        if metric=='sph':
            kp = np.random.randint(low=2,high=N)
        Eini = En(rr,kp,metric=metric,sigma=sigma)
        xrem = rr[kp,0]
        yrem = rr[kp,1]
        rr = rand_increment(rr,kp,Lone=Lone,Ltwo=Ltwo,dfac=dfac,metric=metric)
        Efin = En(rr,kp,metric=metric,sigma=sigma)
        dE = Efin-Eini
        Accept = Metropolis(dE,kBT)
        if Accept:
            move = move+1
        else:
            rr[kp,0]=xrem
            rr[kp,1]=yrem
    E = tot_energy(rr,metric=metric,sigma=sigma)
    # print(rr[0,0],rr[1,0])
    return rr,move,E
#-----------------------------------------#
def rand_increment(rr,kp,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',dfac=64):
    if metric == 'cart':
        rr = increment_cart(rr,kp,Lone=Lone,Ltwo=Ltwo,dfac=dfac)
    elif metric == 'sph':
        rr = increment_sph(rr,kp,dfac=dfac)
    else:
        print('metric',metric,'not coded, exiting')
        quit()
    return rr
#-----------------------------------------#
def increment_cart(rr,kp,Lone=2*np.pi, Ltwo=2*np.pi,dfac=64):
    x0 = rr[kp,0]
    dx = (Lone/dfac)*np.random.uniform(low=-1,high=1.)
    x1 = pbc(x0+dx,Lone)
    y0 = rr[kp,1]
    dy = (Ltwo/dfac)*np.random.uniform(low=-1,high=1.)
    y1 = pbc(y0+dy,Ltwo)
    rr[kp,0] = x1
    rr[kp,1]= y1
    return rr
#-----------------------------------------#
def increment_sph(rr,kp,ntry=10,dfac=64):
#        print(kp)
    th0 = rr[kp,0]
    th1 = get_rand_th1(th0,dfac=dfac)
    phi0 = rr[kp,1]
    dphi = (2*np.pi/dfac)*np.random.uniform(low=-1,high=1.)
    phi1 = pbc(phi0+dphi,2*np.pi)
    rr[kp,0] = th1
    rr[kp,1]= phi1
    return rr
#-----------------------------------------#
def get_rand_th1(th0,dfac=32,ntry=10):
    dth = (np.pi/dfac)*np.random.uniform(low=-1,high=1.)
    th1 = th0+dth
    if th1 > np.pi or th1 < 0 :
        #print(th0,dth,th1)
        th1 = get_rand_th1(th0)
    return th1
#-----------------------------------------#
def pbc(x,Lmax):
    if x > Lmax :
        x = x-Lmax
    if x < 0 :
        x = Lmax-x
    return x
#-----------------------------------------#
def Metropolis(dE,kBT=1.):
    if dE < 0:
        Accept=True
    else:
        rand = np.random.uniform(low=0.0,high=1.)
        if rand < np.exp(-dE/kBT):
            Accept=True
        else:
            Accept=False
    if dE == 0:
        Accept=False
    return Accept
#-----------------------------------------#
def uni_random(N=100,Lone=2*np.pi,Ltwo=2*np.pi):
    ran_array = np.zeros([N,2])
    for ip in range(N):
        ran_array[ip,0] = np.random.uniform(low=0.0,high=Lone)
        ran_array[ip,1] = np.random.uniform(low=0.0,high=Ltwo)
    return ran_array
#-----------------------------------------#
def cont():
    answer=input('continue(Y) ? ')
    if answer == 'Y':
        cont()
    else:
        print('bye')
#-----------------------------------------#
def WritePos(rr,fname='pos'):
    hf = h5py.File(fname+'.h5','w')
    hf.create_dataset('rr',data=rr)
    hf.close()
#-----------------------------------------#
def distance2d(x1,y1,x2,y2,metric='cart'):
    if metric == 'cart':
        ds = cartesian_distance(x1,y2,x2,y2)
    elif metric == 'sph':
        ds = sph_distance(x1,y2,x2,y2)
    else:
        print('metric',metric,'not coded, exiting')
        quit()
    if ds == 0:
        print('distance is zero..')
        print(x1,y2,x2,y2)
        print('quiting')
        quit()
    return ds
#------------------------------------------#
def sph_distance(th1,phi1,th2,phi2):
    """ This is not real distance along geodesics. This
is the 3d distance between the points on the surface of a unit sphere
which is embedded in three dimensions."""
    x1,y1,z1 = cart3d(th1,phi1)
    x2,y2,z2 = cart3d(th2,phi2)
    ds = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 ) 
    return ds
#------------------------------------------#
def cart3d(th,phi):
    x = np.sin(th)*np.cos(phi)
    y = np.sin(th)*np.sin(phi)
    z = np.cos(th)
    return x,y,z
#------------------------------------------#
def cartesian_distance(x1,y1,x2,y2):
    ds = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return ds
#------------------------------------------#
def tot_energy(rr,metric='cart',sigma=None):
    """ The input rr is an array of size (2,N) where 
        N is the number of particle"""
    N=np.shape(rr)[0]
    E = 0
    for i in range(N):
        #print(i)
        Ei =  En(rr,i,metric='cart',sigma=sigma)
        E = E + Ei
        E = E/2. # because every pair is counted twice
    return E
#---------------------------------------------#
def En(rr,kp,metric='cart',sigma=None):
    N = np.shape(rr)[0]
    E = 0
    xk = rr[kp,0]
    yk = rr[kp,1]
    for j in range(N):
        if j != kp:
            xj = rr[j,0]
            yj = rr[j,1]
            ds = distance2d(xk,yk,xj,yj,metric=metric)
            ee = LJ(ds,epsilon=1e-20,sigma=sigma)
            #print(j,ds,ee)
            E = E + ee
    return E
#---------------------------------------------#
def LenardJonesRep(R,epsilon=1e-8,sigma=0.2):
    if R == 0:
        print('zero relative distance betn two particles! quiting')
        quit()
    else:
        if R > sigma:
            V = 0
        else:
            V = epsilon*(sigma/R)**12
    return V
#---------------------------------------------#
def print_xyz(lats,lons,fname='ini_sph',radius=1.):
    phi = lons
    theta = lats+np.pi/2
    x=radius*np.sin(theta)*np.cos(phi)
    y=radius*np.sin(theta)*np.sin(phi)
    z=radius*np.cos(theta)
    hf = h5py.File(fname+'.h5','w')
    hf.create_dataset('x',data=x)
    hf.create_dataset('y',data=y)
    hf.create_dataset('z',data=z)
    hf.close()
#------------------------------------
def lat_lon_list(tri):
    print(tri.lst)
    print(tri.lptr)
    print(tri.lend)
    # we go over the number of particles
    """ We assume the algorithm works in the following way:
    First start randomly with one node. Then go to all the neighbours 
    of that node. The pointer lend tells you when to stop. Once you
    stop assume that the place you stop is your next node and go around
    listing its neighbours. But in this way, every pair of particles 
    should appear only once """
    node_ptr=1
    for ip in range(tri.npts):
        end_ptr = tri.lend[ip]
        kp = tri.lst[node_ptr-1]
        neigh_ptr = tri.lptr[node_ptr-1]
        print("=================================")
        print('particle,end_ptr',kp,end_ptr)
        print("=================================")
        print('Its neighbours are:')
        iter1 = 0
        for iter1 in range(8):
            jp = tri.lst[neigh_ptr-1]
            print('particle,iter1',jp,iter1)
            if neigh_ptr == end_ptr:
                node_ptr = neigh_ptr
                break
            neigh_ptr = tri.lptr[neigh_ptr-1]
#------------------------------------
def lat_lon_neighbour(tri,k):
    # print the latitude and longitude of the kth node and its 6 neighbours
    print(tri.lats[k],tri.lons[k])
    n1 = tri.lst[tri.lptr[k]]
#--------------------------------
def func_moebius():
    u=np.linspace(0,2*np.pi, 24)
    v=np.linspace(-1,1, 8)
    u,v=np.meshgrid(u,v)
    u=u.flatten()
    v=v.flatten()

    #evaluate the parameterization at the flattened u and v
    tp=1+0.5*v*np.cos(u/2.)
    x=tp*np.cos(u)
    y=tp*np.sin(u)
    z=0.5*v*np.sin(u/2.)

    #define 2D points, as input data for the Delaunay triangulation of U
    points2D=np.vstack([u,v]).T
    tri = Delaunay(points2D)#triangulate the rectangle U
    return x,y,z,tri
#------------------------------------------
def map_z2color(zval, colormap, vmin, vmax):
    #map the normalized value zval to a corresponding color in the colormap

    if vmin>vmax:
        raise ValueError('incorrect relation between vmin and vmax')
    t=(zval-vmin)/float((vmax-vmin))#normalize val
    R, G, B, alpha=colormap(t)
    return 'rgb('+'{:d}'.format(int(R*255+0.5))+','+'{:d}'.format(int(G*255+0.5))+\
           ','+'{:d}'.format(int(B*255+0.5))+')'
#-----------------------------------
def map_z2color(zval, colormap, vmin, vmax):
    #map the normalized value zval to a corresponding color in the colormap

    if vmin>vmax:
        raise ValueError('incorrect relation between vmin and vmax')
    t=(zval-vmin)/float((vmax-vmin))#normalize val
    R, G, B, alpha=colormap(t)
    return 'rgb('+'{:d}'.format(int(R*255+0.5))+','+'{:d}'.format(int(G*255+0.5))+\
           ','+'{:d}'.format(int(B*255+0.5))+')'
#---------------------------------
#To plot the triangles on a surface, we set in Plotly Mesh3d the lists of x, y, respectively z- coordinates of the 
#vertices, and the lists of indices, i, j, k, for x, y, z coordinates of all vertices:
#-------------------------------------------------------
def tri_indices(simplices):
    #simplices is a numpy array defining the simplices of the triangularization
    #returns the lists of indices i, j, k

    return ([triplet[c] for triplet in simplices] for c in range(3))
#--------------------------------------------------------------
def plotly_trisurf(x, y, z, simplices, colormap=cm.RdBu, plot_edges=None):
    #x, y, z are lists of coordinates of the triangle vertices
    #simplices are the simplices that define the triangularization;
    #simplices  is a numpy array of shape (no_triangles, 3)
    #insert here the  type check for input data

    points3D=np.vstack((x,y,z)).T
    tri_vertices=map(lambda index: points3D[index], simplices)# vertices of the surface triangles     
    zmean=[np.mean(tri[:,2]) for tri in tri_vertices ]# mean values of z-coordinates of 
                                                      #triangle vertices
    min_zmean=np.min(zmean)
    max_zmean=np.max(zmean)
    facecolor=[map_z2color(zz,  colormap, min_zmean, max_zmean) for zz in zmean]
    I,J,K=tri_indices(simplices)

    triangles=go.Mesh3d(x=x,
                     y=y,
                     z=z,
                     facecolor=facecolor,
                     i=I,
                     j=J,
                     k=K,
                     name=''
                    )

    if plot_edges is None:# the triangle sides are not plotted 
        return [triangles]
    else:
        #define the lists Xe, Ye, Ze, of x, y, resp z coordinates of edge end points for each triangle
        #None separates data corresponding to two consecutive triangles
        lists_coord=[[[T[k%3][c] for k in range(4)]+[ None]   for T in tri_vertices]  for c in range(3)]
        Xe, Ye, Ze=[np.ufunc.reduce(lambda x,y: x+y, lists_coord[k]) for k in range(3)]

        #define the lines to be plotted
        lines=go.Scatter3d(x=Xe,
                        y=Ye,
                        z=Ze,
                        mode='lines',
                        line=dict(color= 'rgb(50,50,50)', width=1.5)
               )
        return [triangles, lines]
#------------------------------------
def triang_moebius():
    x,y,z,tri=func_moebius()
    data1=plotly_trisurf(x,y,z, tri.simplices, colormap=cm.RdBu, plot_edges=True)
    axis = dict(
        showbackground=True,
        backgroundcolor="rgb(230, 230,230)",
        gridcolor="rgb(255, 255, 255)",
        zerolinecolor="rgb(255, 255, 255)",
    )
    layout = go.Layout(
        title='Moebius band triangulation',
        width=800,
        height=800,
        scene=dict(
            xaxis=dict(axis),
            yaxis=dict(axis),
            zaxis=dict(axis),
            aspectratio=dict(
                x=1,
                y=1,
                z=0.5
            ),
        )
    )
    fig1 = go.Figure(data=data1, layout=layout)

    py.iplot(fig1, filename='Moebius-band-trisurf')
#------------------------------------------------------
def sphere(theta, phi): 
    x=np.cos(phi)*np.cos(theta)
    y=np.cos(phi)*np.sin(theta)
    z=np.sin(phi)  
    return x,y,z
#------------------------------------------------------
def barycentric(n):
    return [(i/n,j/n, 1-(i+j)/n) for i in range(n,-1, -1) for j in range(n-i, -1, -1)]
#------------------------------------------------------
def triangles(n, point):
    triplets=[]
    i=0
    j=1
    for nr in range(1,n+1):
        for k in range(nr):
            triplets.append([point[i], point[j], point[j+1]])
            i+=1
            j+=1
        j+=1
    return triplets
#------------------------------------------------------
def divide_sides(t):
    pts=[t[k] for k in range(3)]+[t[0]]
    for k in range(1, 7, 2):
        m=int((k-1)/2)
        pts.insert(k, (t[m]+t[(m+1)%3])/2)
    return pts
#------------------------------------------------------
def line_pts(tri, surface, color='#0000A0'):
    # tri is the list of triangulation sub-triangles
    # surface is the function implementing a surface parameterization
    lines = []
    for t in tri:
        pts=divide_sides(t)
        coords=zip(*[surface(pts[k][0], pts[k][1]) for k in range(7)])
        lines.append(Scatter3d(x=coords[0], 
                               y=coords[1],
                               z=coords[2], 
                               mode='lines', 
                               line=Line(color=color,  
                                         width=5)))
    return lines
#------------------------------------------------------
def triang_sph():
    theta=np.linspace(0,2*np.pi,40)
    phi=np.linspace(-np.pi/2, np.pi/2, 30)
    theta,phi=np.meshgrid(theta,phi)
    x,y,z=sphere(theta,phi)
    A=np.array([0, np.pi/12])
    B=np.array([np.pi/3, np.pi/12])
    C=np.array([np.pi/4, 7*np.pi/16])
    T=[A,B,C]# Parametric triangle of vertices A,B,C
    cartesian_coords=lambda  T, w: w[0]*T[0]+w[1]*T[1]+w[2]*T[2]
    n=8
    bar=barycentric(n)
    pts_tri=[cartesian_coords(T, w) for w in bar]#list of triangulation points
    tri=triangles(n, pts_tri)#list of sub-triangles in T
    trace = Surface(
        z=z,
        x=x,
        y=y,
        colorscale='Viridis',
    )
    axis = dict(
        showbackground=True, 
        backgroundcolor="rgb(230, 230,230)",
        gridcolor="rgb(255, 255, 255)",      
        zerolinecolor="rgb(255, 255, 255)",  
    )
    data = Data([trace]+line_pts(tri, sphere))
    # plot both the sphere as a Surface object and the triangular wireframe
    layout = Layout(
        title='Triangular wireframe on  sphere',
        width=700,
        height=700,
        scene=Scene(  
            xaxis=XAxis(axis),
            yaxis=YAxis(axis), 
            zaxis=ZAxis(axis), 
        ),
        showlegend=False
    )
    
    fig1 = Figure(data=data, layout=layout)
    py.plot(fig1, filename='triangular-wireplot-2')
#-----------------------------------------