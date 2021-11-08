import os as os
import pencil_old as pc
import numpy as np
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
#import matplotlib.ticker as mticker
#from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.interpolate import interp1d
from numpy import linalg as LA
import h5py
#-----------------------------------------#
def MC_surf(N,Lone=2*np.pi,Ltwo=2*np.pi,metric='cart',maxiter=100,kBT=1.,
            interactive=True,dfac=64):
    rrini = inidist(N,Lone,Ltwo,metric=metric)
    rr = rrini
    print(rr[0,:])
    WritePos(rr,fname='ini_pos')
    cont=True
    while cont:
        Eini = tot_energy(rr,metric=metric)
        print('Eini/N=',Eini/N)
        if interactive:
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(rr[:,0],rr[:,1],'o',color='C0')
        rr,move,Efin = MC_step(rr,Lone=Lone,Ltwo=Ltwo,metric=metric,
                               maxiter=maxiter,kBT=kBT,dfac=dfac)
        print('MC steps',maxiter)
        print('accepted moves',move)
        print('Efin/N=',Efin/N)
        WritePos(rr,fname='fin_pos')
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
            print('Run over, go home')
            break
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
    ran = np.zeros([N,2])
    # the first two points are always at the two poles
    ran[0,0]=0
    ran[0,1]=0
    ran[1,0]=np.pi
    ran[1,1]=0
    for ip in range(2,N):
        rth = np.random.uniform(low=0.0,high=np.pi)
        if rth==0:
            rth = rth+ np.random.uniform(low=0.0,high=np.pi)
        if rth==np.pi:
            rth = rth- np.random.uniform(low=0.0,high=np.pi)
        ran[ip,0] = rth
        ran[ip,1]=np.random.uniform(low=0.0,high=2*np.pi)
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
            dfac=64):
    N = np.shape(rr)[0]
    move = 0
    for iter in range(maxiter):
        """ Choose a random integer between 0 and N-1 """
        kp = np.random.randint(low=0,high=N)
        if metric=='sph':
            kp = np.random.randint(low=2,high=N)
        Eini = En(rr,kp,metric=metric)
        xrem = rr[kp,0]
        yrem = rr[kp,1]
        rr = rand_increment(rr,kp,Lone=Lone,Ltwo=Ltwo,dfac=dfac,metric=metric)
        Efin = En(rr,kp,metric=metric)
        dE = Efin-Eini
        Accept = Metropolis(dE,kBT)
        if Accept:
            move = move+1
        else:
            rr[kp,0]=xrem
            rr[kp,1]=yrem
    E = tot_energy(rr,metric=metric)
    print(rr[0,0],rr[1,0])
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
    for itry in range(ntry) :
        dth = (np.pi/dfac)*np.random.uniform(low=-1,high=1.)
        th1 = th0+dth
        if th1 < np.pi or th1 > 0 :
            break
        else:
            if itry == ntry-1:
                print(ntry,' tries failed to move a particle..quiting')
                quit()
    phi0 = rr[kp,1]
    dphi = (2*np.pi/dfac)*np.random.uniform(low=-1,high=1.)
    phi1 = pbc(phi0+dphi,2*np.pi)
    rr[kp,0] = th1
    rr[kp,1]= phi1
    return rr
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
def tot_energy(rr,metric='cart'):
    """ The input rr is an array of size (2,N) where 
N is the number of particle"""
    N=np.shape(rr)[0]
    E = 0
    for i in range(N):
        #print(i)
        Ei =  En(rr,i,metric='cart')
        E = E + Ei
        E = E/2. # because every pair is counted twice
    return E
#---------------------------------------------#
def En(rr,kp,metric='cart'):
    N = np.shape(rr)[0]
    E = 0
    xk = rr[kp,0]
    yk = rr[kp,1]
    for j in range(N):
        if j != kp:
            xj = rr[j,0]
            yj = rr[j,1]
            ds =  distance2d(xk,yk,xj,yj,metric=metric)
            ee = LenardJonesRep(ds)
            #print(j,ds,ee)
            E = E + ee
    return E
#---------------------------------------------#
def LenardJonesRep(R,epsilon=1e-10,sigma=1.):
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
