import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
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
def partition(energy,KbT=1,running_avg=1):
    '''Returns the running average of partition function i.e. Z=<exp(-beta*E_tot)>'''
    beta=1.0/KbT
    Nensemble=energy.shape[0]
    # O(n^2) algorithm
    if running_avg:
        ZZ=np.zeros(Nensemble)
        Eminis=np.zeros(Nensemble)
        for i in range(0,Nensemble):
            Emin=np.min(energy[0:i+1])
            ZZ[i]=np.sum(np.exp(-beta*(energy[0:i+1] - Emin)))/(i+1)
            Eminis[i]=Emin
        return Eminis,ZZ
    else:
        Emin=np.min(energy)
        ZZ=np.sum(np.exp(-beta*(energy - Emin)))/Nensemble
        return Emin,ZZ
#------------------------------------------------------------------------------------#
def free_energy(energy,KbT=1,running_avg=1):
    '''Returns the free energy i.e. Z=<exp(-beta*E_tot)>, F= -KbT*np.log10(Z)'''
    beta = 1.0/KbT
    Nens = energy.shape[0]
    Emin = -1e-16
    # FF = Emin
    Eminis,ZZ=partition(energy,KbT=KbT,running_avg=running_avg)
    # Free energy = -KbT*log(ZZ)
    FF = Eminis - KbT*np.log(ZZ)
    return Eminis, FF
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