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
        running=isrunning()
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
def avg_quantity_tz(tz_start=None,tz_end=None,index=2,step=0.01,
                    datadir="./",subfol="rerun/",start=1000,
                    nch=20,error=True):
    nruns=int((tz_start-tz_end)/step)+2
    tzall=np.linspace(tz_start,tz_end,nruns)
    mc_log=np.empty(nruns,dtype=object)
    for ifol,tz in enumerate(tzall):
        folder = str(float("{0:.3f}".format(tz)))+"/"
        mc_log[ifol]=np.loadtxt(datadir+folder+subfol+"/mc_log")
    # ------------------ compute average -------------------- #
    dvert=tzall[0]-tzall
    mc_noafm = np.loadtxt(datadir+"noafm/"+subfol+"/mc_log")
    baseE = np.mean(mc_noafm[start:,index])
    # print(baseE)
    avgE=np.zeros(nruns)
    std=np.zeros(nruns)
    err=np.zeros(nruns)
    for ifol in range(nruns):
        tot_ener=mc_log[ifol][start:,index]
        Nmc=tot_ener.shape[0]
        chsz=int(Nmc/nch)
        Ech=np.zeros(chsz)
        for ich in range(nch):
            Ech[ich] = np.mean(tot_ener[ich*chsz:(ich+1)*chsz])
        err[ifol]=np.std(Ech)
        # print(err[ifol])
        avgE[ifol]=np.mean(Ech)
    if error:
        return dvert,avgE,baseE,err
    else:
        return dvert,avgE,baseE
# # ---------------------------------------------------------------- #
# def FF_tz(tz_start,tz_end,step=0.02,datadir="../"):
#     nruns=int((tz_start-tz_end)/0.02)+2
#     tzall=np.linspace(tz_start,tz_end,nruns)
#     mc_log=np.empty(nruns,dtype=object)
#     for ifol,tz in enumerate(tzall):
#         folder = str(float("{0:.3f}".format(tz)))
#         mc_log[ifol]=np.loadtxt(datadir+folder+"/rerun/mc_log")
#     # ------------------ compute mean_energy -------------------- #
#     dvert=tzall[0]-tzall
#     avgE=np.zeros(nruns)
#     std=np.zeros(nruns)
#     for ifol in range(nruns):
#         tot_ener=mc_log[ifol][:,]
#         avgE[ifol] = np.mean(tot_ener)
#         std[ifol] = np.std(tot_ener)
#     return dvert,avgE
# # ---------------------------------------------------------------- #