import numpy as np
from numpy import linalg as LA
#-----------------------------------------
def voronoi_area(cotJ,cotK,jsq,ksq,area):
    '''Q. Given two cotangent angles, it returns either the area due to perpendicular bisector,
        or the barycenter.'''
    sigma=0
    if cotJ>0 and cotK>0:
        if cotJ*cotK <1:
            sigma = 0.125*(cotJ*jsq+cotK*ksq)
        else:
            sigma = 0.5*area
    else:
        sigma = 0.25*area
    return sigma;
#-----------------------------------------
def foldername(file):
    '''The function return the folder name for a given filename'''
    x=file.split("/")
    name=""
    for n in x[0:-1]:
        name=name+n+"/"
    return name
#-----------------------------------------
def partition(energy,KbT):
    '''Returns the running average of partition function i.e. Z=<exp(-beta*E_tot)>'''
    beta=1.0/KbT
    Nensemble=energy.shape[0]
    ZZ=np.zeros(Nensemble)
    Emin=np.min(energy)
    # Zsum=np.exp(-beta*energy[0])
    ZZ[0]=Zsum
    for i in range(1,Nensemble):
        Zsum=Zsum+np.exp(-beta*energy[i])
        ZZ[i]=Zsum/(i+1)
    return ZZ
#-----------------------------------------
def free_energy(energy,KbT=1):
    '''Returns the running average of free energy i.e. Z=<exp(-beta*E_tot)>, F= -KbT*np.log10(Z)'''
    beta = 1.0/KbT
    Nens = energy.shape[0]
    Emin = -1e-16
    # FF = Emin
    Zred=0
    for i in range(Nens):
        Emin_old = Emin
        if Emin<FF[i]:
            Emin=FF[i]
        # Zred=Zred*np.exp(-beta*)
        FF=Emin+np.log10(Zred)
    return -KbT*np.log10(partition(energy,KbT=KbT))
#-----------------------------------------