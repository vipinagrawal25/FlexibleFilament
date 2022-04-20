import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time
import paramio as pio
import os
import psutil
from scipy.spatial import SphericalVoronoi
#------------------------------------------------------------------------------------#
def voronoi_area(cotJ,cotK,jsq,ksq,area):
    '''Q. Given two cotangent angles, it returns either the area due to perpendicular 
        bisector, or the barycenter.'''
    sigma=0
    if cotJ>0 and cotK>0:
        if cotJ*cotK <1:
            sigma = 0.125*(cotJ*jsq+cotK*ksq)
        else:
            sigma = 0.5*area
    else:
        sigma = 0.25*area
    return sigma;
#------------------------------------------------------------------------------------#
def foldername(file):
    '''The function return the folder name for a given filename'''
    x=file.split("/")
    name=""
    for n in x[0:-1]:
        name=name+n+"/"
    return name
#------------------------------------------------------------------------------------#
def partition(energy,KbT=1):
    '''Returns the running average of partition function i.e. Z=<exp(-beta*E_tot)>'''
    beta=1.0/KbT
    Nensemble=energy.shape[0]
    # O(n^2) algorithm
    Emin = np.min(energy)
    ZZ=np.zeros(Nensemble)
    ZZ = np.cumsum(np.exp(-beta*(energy - Emin)))
    return Emin,ZZ
#------------------------------------------------------------------------------------#
def free_energy(energy,KbT=1):
    '''Returns the free energy i.e. Z=<exp(-beta*E_tot)>, F= -KbT*np.log10(Z)'''
    beta = 1.0/KbT
    Nens = energy.shape[0]
    # FF = Emin
    Emini,ZZ=partition(energy,KbT=KbT)
    FF = Emini - KbT*np.log(ZZ)
    return Emini, FF
#------------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------------#
def isrunning(procname='run'):
    for p in psutil.process_iter(['username','name']):
        if p.info['name']==procname:
            running = 1
        else:
            running = 0
    return running
#---------------------------------------------------------------- #
def wait(procname='run',timedelay=10):
    running=1
    while running==1:
        time.sleep(timedelay)
        running=isrunning(procname=procname)
#---------------------------------------------------------------- #
def movetip(tz_start,tz_end,step=-0.01,timedelay=10,restart=None):
    tz_all=np.arange(tz_start,tz_end,step)
    for tz in tz_all:
        pio.change_param(tip_pos_z=tz)
        print("# tz = ", tz)
        g = str(float("{0:.3f}".format(tz)))
        os.system("mkdir "+g)
        if restart is None:
            pio.change_param(is_restart=0)
        else:
            pio.change_param(is_restart=1)
            os.system("cp "+restart+"/restart.h5 "+g)
        os.system("./run para_file.in "+ g + "> "+g+"/terminal.txt")
        restart=g
        wait()
#---------------------------------------------------------------- #
def avg_quantity_tz(folders,index=2,datadir="./",subfol="rerun/",start=1000,
                    nch=10,error=True,nopush="noafm/",index2=None,Ks=False):
    if index2 is None:
        index2=index
    nruns=len(folders)
    mc_log=np.empty(nruns,dtype=object)
    for i,fol in enumerate(folders):
        mc_log[i]=np.loadtxt(datadir+fol+subfol+"/mc_log")
    # ------------------ compute average -------------------- #
    dvert=np.zeros(nruns)
    mc_nopush = np.loadtxt(datadir+nopush+subfol+"/mc_log")
    zoavg=np.mean(mc_nopush[start:,-2])-np.mean(mc_nopush[start:,-1])
    baseE = np.mean(mc_nopush[start:,index])
    avgE=np.zeros(nruns)
    std=np.zeros(nruns)
    err=np.zeros(nruns)
    for ifol in range(nruns):
        tot_ener=np.abs(mc_log[ifol][start:,index])+np.abs(mc_log[ifol][start:,index2])
        tot_ener=0.5*tot_ener
        zavg=np.mean(mc_log[ifol][start:,-2])-np.mean(mc_log[ifol][start:,-1])
        dvert[ifol]=1-zavg/zoavg;
        #
        Nmc=tot_ener.shape[0]
        chsz=int(Nmc/nch)
        Ech=np.zeros(nch)
        for ich in range(nch):
            Ech[ich] = np.mean(tot_ener[ich*chsz:(ich+1)*chsz])
        err[ifol]=np.std(Ech)
        avgE[ifol]=np.mean(Ech)
    if error:
        return dvert,avgE,baseE,err
    else:
        return dvert,avgE,baseE
# # ---------------------------------------------------------------- #
def Ks_vs_Y3d(datadir="./",subsubfol="rerun/",start=1000,
                nch=10,error=True,nopush="noafm/",fit_end=10):
    fols=sorted(glob.glob(datadir+"/run"))[0:-1]
    fols = [ifol.replace("/run","/") for ifol in fols]
    Y3d = [float(ifol.replace("/run",""))*(10**9)/3 for ifol in fols]
    mm=np.zeros(len(fols))
    bb=np.zeros(len(fols))
    dvert_all=np.empty(len(fols),dtype=object)
    for i,ifol in enumerate(fols):
        os.chdir(ifol)
        tzall=sorted(glob.glob("*/mc_log"))
        step=float(tzall[1].replace("/mc_log",""))-float(tzall[0].replace("/mc_log",""))
        tz_start=float(tzall[-4].replace("/mc_log",""))
        tz_end=float(tzall[0].replace("/mc_log",""))
        nruns=int(abs(tz_start-tz_end)/step+2)
        folders=np.linspace(tz_start,tz_end,nruns)
        dvert,force,baseE,error=avg_quantity_tz(folders=folders,index=-4,datadir="./",
                                                subfol=subsubfol)
        m,b = np.polyfit(dvert[0:fit_end], force[0:fit_end], 1)
        # plt.plot(dvert,m*dvert+b,'-')
        print(m)
        mm[i]=m
        bb[i]=b
        dvert_all[i]=dvert
        # plt.legend()
        os.chdir("../")
    return Y3d, mm, bb, dvert_all
#-------------------------------------------------------------------#
def dvert(folders,datadir="./",subfol="rerun/",start=1000,end=-1,
          nch=10,error=True,nopush="noafm/"):
    nruns=len(folders)
    mc_nopush = np.loadtxt(datadir+nopush+subfol+"/mc_log")
    zoavg=np.mean(mc_nopush[start:end,-2])-np.mean(mc_nopush[start:end,-1])
    dvert=np.zeros(nruns)
    err=np.zeros(nruns)
    mc_log=np.empty(nruns,dtype=object)
    for i,fol in enumerate(folders):
        mc_log[i]=np.loadtxt(datadir+fol+subfol+"/mc_log")
    for ifol in range(nruns):
        zz=mc_log[ifol][start:end,-2]-mc_log[ifol][start:end,-1]
        Nmc=zz.shape[0]
        chsz=int(Nmc/nch)
        zch=np.zeros(nch)
        for ich in range(nch):
            zch[ich] = np.mean(zz[ich*chsz:(ich+1)*chsz])
        print(np.mean(zch))
        err[ifol]=1-np.std(zch)/zoavg;
        dvert[ifol]=1-np.mean(zch)/zoavg;
    if error:
        return dvert,err
    else:
        return dvert
#-------------------------------------------------------------------#
def isgaussian(datadir="./",subfol="./",index=2,start=1000):
    mc_log=np.loadtxt(datadir+"/"+subfol+"/mc_log")
    data=mc_log[start:,index]
    avg=np.mean(mc_log[start:,index])
    std=np.std(mc_log[start:,index])
    nbin=50
    hist=np.histogram((data-avg)/std,bins=nbin,density=True)
    cen = hist[1][1:] + hist[1][:-1]
    plt.subplots()
    plt.plot(cen*0.5,hist[0],'.-')
    xx = (np.linspace(0,nbin,nbin+1)-nbin/2)*8/nbin
    f = 1/(np.sqrt(2*np.pi))*np.exp(-0.5*xx**2)
    plt.semilogy(xx,f)
    plt.show()
#-------------------------------------------------------------------#
def Ks_spring(datadir="Zn1.0/",KbT=1e0,subfol="./",start=100000):
    mc_log=np.loadtxt(datadir+subfol+"/mc_log")
    z0=mc_log[start:,-2]-mc_log[start:,-1]
    print(np.mean(z0))
    return KbT/np.var(z0)