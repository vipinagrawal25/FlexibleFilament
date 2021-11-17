from __future__ import division
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
import plotly as pl
import plotly.graph_objs as go
import matplotlib.cm as cm
from scipy.spatial import Delaunay
from scipy.spatial import SphericalVoronoi, geometric_slerp
from mpl_toolkits.mplot3d import proj3d
#import plotly as py
#from plotly.graph_objs import *
#import plotly.tools as tls
#-----------------------------------------#
def SphVoronoi(rr,R=1,lplot=True):
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
    ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2],
                   c='g')
    #ax.axhline(0.5, ls=':')
    ax.scatter(points[0:2, 0], points[0:2, 1], points[0:2, 2], c='k')
    for region in sv.regions:
        n = len(region)
        for i in range(n):
            start = sv.vertices[region][i]
            end = sv.vertices[region][(i + 1) % n]
            print(i)
            print(start)
            print(end)
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
    return rr
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
#To plot the triangles on a surface, we set in Plotly Mesh3d the lists of x, y, respectively z- coordinates of the vertices, and the lists of indices, i, j, k, for x, y, z coordinates of all vertices:
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
    
    data = Data([trace]+line_pts(tri, sphere))#plot both the sphere as a Surface object and the triangular wireframe
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
