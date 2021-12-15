import numpy as NP
import matplotlib.pyplot as plt
import matplotlib
from scipy.fftpack import *
matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D
#---------------------------------------#
def length(arr):
    ll = 0.0
    for iN in range(1,arr.shape[0]):
        ll = ll + NP.sqrt( (arr[iN,0]-arr[iN-1,0])**2 + (arr[iN,1] - arr[iN-1,1])**2 + (arr[iN,2]-arr[iN-1,2])**2)
    return ll

def snap_one(ax,dirname,omega,snaptime,symbol='o-',txt=False,legend=True):
    TT=2*NP.pi/omega
    isnap=int(NP.around(50*TT*snaptime))
    dd = NP.loadtxt(dirname+"/output/var"+str(isnap)+".txt")
    ax.plot(dd[:,1],dd[:,2],symbol,label =str(snaptime)+'T')
    if(txt):
        ax.text(dd[0,1]-0.1,dd[0,2]-0.1, str((snaptime))+'T',fontsize=20,weight='bold')
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    if(legend):
        ax.legend(fontsize=15)
    ax.set_ylim()
    ax.tick_params(axis='x',labelsize=17)
    ax.tick_params(axis='y',labelsize=17)
    return ax

def snap(ax,dirname,omega,snaptime,symbol='o-',txt=False,legend=True):
    TT=2*NP.pi/omega
    nfigs = NP.shape(snaptime)[0]
    for ifig in range(nfigs):
        isnap=int(NP.around(50*TT*snaptime[ifig]))
        dd = NP.loadtxt(dirname+"/output/var"+str(isnap)+".txt")
        ax.plot(dd[:,1],dd[:,2],symbol,label =str(snaptime[ifig])+'T')
        if(txt):
            ax.text(dd[0,1]-0.1,dd[0,2]-0.1, str(int(snaptime[ifig]))+'T',fontsize=20,weight='bold')
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    if(legend):
        ax.legend(fontsize=15)
    ax.set_ylim()
    ax.tick_params(axis='x',labelsize=17)
    ax.tick_params(axis='y',labelsize=17)
    return ax

def curv_compare(ax,dirname,omega,snaptime):
    nfigs = NP.shape(snaptime)[0]
    TT = 2*NP.pi/omega
    ss = 1.28*NP.loadtxt(dirname+"/output/material_point.txt")
    curv = 1/NP.sqrt(2)*NP.loadtxt(dirname+"/output/kappa.txt")
    for ifig in range(nfigs):
        isnap = int(NP.around(50*TT*snaptime[ifig]))
        ax.plot(ss[isnap,:],curv[isnap,1:],'.-',label =str(snaptime[ifig])+'T')
    ax.legend(fontsize=15)
    plt.rcParams["mathtext.fontset"] = "cm"
    ax.set_ylabel(r'$a\kappa$', fontsize = 23)
    ax.set_xlabel('$s$',fontsize = 23)
    ax.tick_params(axis='x',labelsize=17)
    ax.tick_params(axis='y',labelsize=17)
    return ax

def VtrPlot(ax,dirname,omega,snapmin,snapmax,VtrPoint=2,ind1=0,ind2=2,marker='s',label=" "):
    labelname = ['$V_z/V_{max}$','$V_x/V_{max}$','$V_y/V_{max}$']
    TT = 2*NP.pi/omega
    tmin=int(TT*snapmin*50)
    tmax=int(TT*snapmax*50)
    Vtr = 1/2.56*NP.load(dirname+"/output/Vtracer" + str(VtrPoint) +".npy")
    if label == " ":
        ax.plot(Vtr[tmin:tmax,ind1],Vtr[tmin:tmax,ind2],marker,ms=3)
    else:
        ax.plot(Vtr[tmin:tmax,ind1],Vtr[tmin:tmax,ind2],marker,ms=3,label=label)
    ax.set_xlabel(labelname[ind1],fontsize=20)
    ax.set_ylabel(labelname[ind2],fontsize=20)
    ax.tick_params(axis='x',labelsize=17)
    ax.tick_params(axis='y',labelsize=17)
    return ax

def VtrPlot3D(ax,dirname,omega,snapmin,snapmax,marker='s-'):
    labelname = ['$V_x/V_{max}$','$V_y/V_{max}$','$V_z/V_{max}$']
    TT = 2*NP.pi/omega
    tmin=int(TT*snapmin*50)
    tmax=int(TT*snapmax*50)
    Vtr = 1/2.56*NP.load(dirname+"/output/Vtracer" + str(VtrPoint) +".npy")
    ax.plot(Vtr[tmin:tmax,0],Vtr[tmin:tmax,1],Vtr[tmin:tmax,2],marker,ms=2)
    ax.set_xlabel(labelname[0],fontsize=23)
    ax.set_ylabel(labelname[1],fontsize=23)
    ax.set_zlabel(labelname[2],fontsize=23)
    ax.grid(True)
    return ax

def LangVTrPlot(ax,dirname,omega,icy):
    TT = 2*NP.pi/omega
    dd0= NP.loadtxt(dirname+'output/tracer0.txt')[:,0:3] 
    dd = NP.loadtxt(dirname+'output/tracer'+str(int(icy*TT*50))+'.txt')[:,0:3]
    nsnap = int(dd0.shape[0]/2)
    colors = NP.zeros(nsnap)
    for isnap in range(nsnap):
        colors[isnap] = isnap
    ax.scatter(dd[0:256,0],dd[0:256,1],dd[0:256,2],'o',c=colors,cmap='gnuplot2')
    ax.tick_params(axis='x',labelsize=17,rotation=30)
    ax.tick_params(axis='y',labelsize=17)
    ax.tick_params(axis='z',labelsize=17)
    ax.xaxis.labelpad = 30
    ax.yaxis.labelpad = 15
    ax.zaxis.labelpad = 15
    #
    ax.set_xlabel('X',fontsize=23)
    ax.set_ylabel('Y',fontsize=23)
    ax.set_zlabel('Z',fontsize=23)
    return ax

def SineTransform(ax,dirname,omega,snaptime,length=1.28,dia=0.005):
    TT = 2*NP.pi/omega
    nfigs = NP.shape(snaptime)[0]
    curv = NP.loadtxt(dirname+'/output/kappa.txt')
    for ifig in range(nfigs):
        isnap = int(NP.around(50*TT*snaptime[ifig]))
        dd = curv[isnap,1:]
        NN=dd.size
        curvsine=1/(NN)*dst(dd,type=1)
        wavenrs = length/(dia*NN)*NP.linspace(0,NN,NN)
        ax.plot(wavenrs,curvsine,'.-',label= str(snaptime[ifig])+'T')
        
    plt.rcParams["mathtext.fontset"] = "cm"
    ax.set_ylabel('Amplitude', fontsize = 23)
    ax.set_xlabel('L/D X Wave Numbers',fontsize = 23)
    ax.tick_params(axis='x',labelsize=17)
    ax.tick_params(axis='y',labelsize=17)
    ax.legend(fontsize=17)
    return ax

def tracerPlot(ax,dirname,omega,cy,ntr=-1):
    TT = 2*NP.pi/omega
    dd0= NP.loadtxt(dirname+'output/tracer0.txt')[:,0:3]
    dd1 = NP.loadtxt(dirname+'output/tracer'+str(int(cy*TT*50))+'.txt')[:,0:3]
    #
    nsnap = int(dd0.shape[0])
    colors = NP.zeros(nsnap)
    for isnap in range(nsnap):
        colors[isnap] = NP.linalg.norm(dd0[isnap,:]-NP.array([0.02,0.,0.64]))
    if ntr==-1:
        ax.scatter(dd1[ind,0],dd1[ind,1],dd1[ind,2],c=colors,marker='.',s=40,cmap='rainbow')
    else:
        ind = minTracerDisplacementIndex(dirname,ntr,cy1=0,cy2=cy)
        ax.scatter(dd1[ind,1],dd1[ind,2],dd1[ind,0],c=colors[0:ntr],s=40,marker='.',cmap='rainbow')
    #
    ax.locator_params(axis='x', nbins=3)
    ax.locator_params(axis='y', nbins=4)
    ax.locator_params(axis='z', nbins=3)
    ax.set_xlabel('x',fontsize=size1)
    ax.set_ylabel('y',fontsize=size1)
    ax.set_zlabel('z',fontsize=size1)
    return ax

def minTracerDisplacementIndex(dirname,omega,maxtracer,cy1=0,cy2=-1):
    TT = 2*NP.pi/omega
    snap1 = int(cy1*TT*50)
    dd1 = NP.loadtxt(dirname+'output/tracer'+str(snap1)+'.txt')[:,0:3]
    if cy2==-1:
        nsnap = NP.loadtxt(dirname+'output/time.txt').shape[0]
        cy2 = int(nsnap/(50*TT))
    snap2 = int(cy2*TT*50)
    dd2 = NP.loadtxt(dirname+'output/tracer'+str(snap2)+'.txt')[:,0:3]
    Np = dd1.shape[0]
    diff = NP.zeros([Np])
    maxdiff = NP.zeros(maxtracer)
    for ip in range(Np):
        diff[ip] = NP.linalg.norm(dd1[ip,:]-dd2[ip,:])
    ind = NP.sort(NP.argsort(diff)[0:maxtracer])
    return ind

def DiffCoeff(dirname,omega,cy1=0,cy2=-1,ind=[],nparticles=256):
    nsnap = NP.loadtxt(dirname+'output/time.txt').shape[0]
    TT = 2*NP.pi/omega
    if cy2==-1:
        ncy = int(nsnap/(50*TT))
    else:
        ncy = cy2
    MSDTr = NP.zeros([ncy-cy1])
    if ind==[]:
        ind = NP.linspace(0,nparticles-1,nparticles)
        ind = ind.astype(int)
    dd0 = NP.loadtxt(dirname+'output/tracer'+str(int(50*TT*cy1))+'.txt')[ind,1:3]
    for icy in range(cy1+1,ncy):
        ddcy = NP.loadtxt(dirname+'output/tracer'+str(int(50*TT*icy))+'.txt')
        ddcy=ddcy[ind,1:3]
        MSDTr[icy-cy1] = NP.linalg.norm(ddcy-dd0)**2/nparticles
    return MSDTr

def posneg(array):
    bol = array > 0
    pos = array[bol]
    bol = array < 0
    neg = array[bol]
    return pos,neg

def rank_order(x):
    N = NP.size(x)
    xout = NP.sort(x)
    cpdf = NP.arange(N,0,-1)/float(N)
    return xout,cpdf

def PDF_disp(dirname,omega,ind=[],cy1=0,cy2=-1): 
    TT = 2*NP.pi/omega
    if cy2==-1:
        nsnap = NP.loadtxt(dirname+'output/time.txt').shape[0]
        ncy = int(nsnap/(50*TT))
    else:
        ncy=cy2
    if ind==[]:
        nparticles = NP.loadtxt(dirname+'output/tracer0.txt').shape[0]
        ind = NP.linspace(0,nparticles-1,nparticles)
        ind = ind.astype(int)
    else:
        nparticles = ind.shape[0]
    alldisp = NP.zeros([nparticles*(ncy-1-cy1),3])
    ix=0
    dd1 = NP.loadtxt(dirname+'output/tracer'+str(int(50*TT*(cy1+1)))+'.txt')[ind,0:3]
    for icy in range(cy1,ncy-1):
        dd2 = NP.loadtxt(dirname+'output/tracer'+str(int(50*TT*(icy+1)))+'.txt')[ind,0:3]
        disp = (dd2-dd1)
        dd1 = NP.copy(dd2)
        for ip in range(nparticles):
            alldisp[ix] = disp[ip]
            ix=ix+1
    dia = 0.005
    alldisp=alldisp/dia
    dxp,dxn=posneg(alldisp[:,1])
    xp,pdfxp = rank_order(dxp)
    xn,pdfxn = rank_order(-dxn)
    dyp,dyn=posneg(alldisp[:,2])
    yp,pdfyp = rank_order(dyp)
    yn,pdfyn = rank_order(-dyn)
    dzp,dzn=posneg(alldisp[:,0])
    zp,pdfzp = rank_order(dzp)
    zn,pdfzn = rank_order(-dzn)
    return xp,pdfxp,xn,pdfxn,yp,pdfyp,yn,pdfyn,zp,pdfzp,zn,pdfzn,alldisp