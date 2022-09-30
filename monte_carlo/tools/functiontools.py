import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time
import paramio as pio
import os
import psutil
from scipy.spatial import SphericalVoronoi
from multipledispatch import dispatch
import glob
import visitio as vio
import numpy.linalg as la
from scipy.special import sph_harm
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
    Emini,ZZ=partition(energy,KbT=KbT)
    FF = Emini - KbT*np.log(ZZ[-1])
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
def avg_quantity_tz(folders=None,mc_log=None,index=2,datadir="./",subfol="rerun/",start=0,
                    nch=10,error=True,nopush="noafm/",index2=None,Ks=False):
    if index2 is None:
        index2=index
    if mc_log is None:
        mc_log= read_mc_log(folders,datadir=datadir,subfol=subfol)
    # ------------------ compute average -------------------- #
    nruns=mc_log.shape[0]
    avgE=np.zeros(nruns)
    std=np.zeros(nruns)
    err=np.zeros(nruns)
    for ifol in range(nruns):
        if index==0:
            tot_ener=np.abs(mc_log[ifol][start:])
        else:
            tot_ener=np.abs(mc_log[ifol][start:,index])+np.abs(mc_log[ifol][start:,index2])
            tot_ener=0.5*tot_ener
        #
        Nmc=tot_ener.shape[0]
        chsz=int(Nmc/nch)
        Ech=np.zeros(nch)
        for ich in range(nch):
            Ech[ich] = np.mean(tot_ener[ich*chsz:(ich+1)*chsz])
        err[ifol]=np.std(Ech)
        avgE[ifol]=np.mean(Ech)
    if error:
        return avgE,err
    else:
        return avgE
#---------------------------------------------------------------- #
def read_mc_log_one(folder,start=0,fname='mc_log'):
    try:
        mc_log=np.load(folder+fname+".npy")[start:]
    except:
        dd=np.loadtxt(folder+fname)
        np.save(folder+fname,dd)
        mc_log=dd[start:]
    return mc_log
#---------------------------------------------------------------- #
def read_mc_log(folders,datadir="",subfol="./",start=0,fname='mc_log'):
    '''It reads the data and returns mc_log for all the folders'''
    nruns=len(folders)
    if len(folders)==1:
        filein=datadir+folders+subfol
        mc_log=read_mc_log_one(filein,start=start,fname=fname)
    mc_log=np.empty(nruns,dtype=object)
    for i,fol in enumerate(folders):
        filein=datadir+fol+subfol
        mc_log[i]=read_mc_log_one(filein,start=start,fname=fname)
    return mc_log
#---------------------------------------------------------------- #
def Ks_vs_Y3d(datadir="./",subsubfol="rerun/",start=0,
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
        print(m)
        mm[i]=m
        bb[i]=b
        dvert_all[i]=dvert
        # plt.legend()
        os.chdir("../")
    return Y3d, mm, bb, dvert_all
#-------------------------------------------------------------------#
def dvert(folders=None,mc_log=None,datadir="./",subfol="rerun/",start=0,
          nch=10,error=True,nopush="noafm/",index1=-3,index2=-2):
    if mc_log is None:
        mc_log=read_mc_log(folders,datadir=datadir,subfol=subfol)
    nruns=mc_log.shape[0]
    mc_nopush = np.loadtxt(datadir+nopush+subfol+"/mc_log")
    zoavg=np.mean(mc_nopush[start:,index1])-np.mean(mc_nopush[start:,index2])
    dvert=np.zeros(nruns)
    err=np.zeros(nruns)
    for ifol in range(nruns):
        zz=mc_log[ifol][start:,index1]-mc_log[ifol][start:,index2]
        Nmc=zz.shape[0]
        chsz=int(Nmc/nch)
        zch=np.zeros(nch)
        for ich in range(nch):
            zch[ich] = np.mean(zz[ich*chsz:(ich+1)*chsz])
        # print(zch)
        err[ifol]=np.std(1-zch/zoavg)
        dvert[ifol]=np.mean(1-zch/zoavg)
    if error:
        return dvert,err
    else:
        return dvert
#-------------------------------------------------------------------#
def isgaussian(folders=None,mc_log=None,datadir="./",subfol="./",index=2,start=0):
    if mc_log is None:
        mc_log= read_mc_log(folders,datadir=datadir,subfol=subfol)
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
#-------------------------------------------------------------------------------------------#
def Ks_spring(folders=None,mc_log=None,datadir="./",subfol="./",KbT=1e0,start=10000,nch=10,
              error=True,index1=-3,index2=-2):
    if mc_log is None:
        mc_log = read_mc_log(folders,datadir=datadir,subfol=subfol)
    nruns = mc_log.shape[0]
    Ks=np.zeros(nruns)
    err=np.zeros(nruns)
    for irun in range(nruns):
        z0=mc_log[irun][start:,index1]-mc_log[irun][start:,index2]
        ## Divide the data into 10 bins
        Nmc=zz.shape[0]
        chsz=int(Nmc/nch)
        Ksch=np.zeros(nch)
        for ich in range(nch):
            Ksch[ich] = KbT/np.var(zz[ich*chsz:(ich+1)*chsz])
        err[irun]=np.std(Ksch)
        Ks[irun]=np.mean(Ksch)
    return Ks,err
#-------------------------------------------------------------------------------------------#
def height_field(files,Np=None,nbin=50,radius=1):
    if Np is None:
        pfile = foldername(files[0])+"/para_file.in"
        pdict=pio.read_param(pfile)
        Np,radius=pdict['N'],pdict['radius']
    hfluc_all=[]
    for i,file in enumerate(files):
        pos=np.loadtxt(file,skiprows=5,max_rows=Np)
        hfluc=np.sqrt(pos[:,0]**2+pos[:,1]**2+pos[:,2]**2) - radius
        hfluc_all.extend(hfluc)
    hist=np.histogram(hfluc_all,bins=nbin,density=True)
    return hist
# -----------------------------------------------------------------#
def cart2sph(points):
    ''' Computes theta and phi for the given points.'''
    Np=points.shape[0]
    theta = np.zeros(Np)
    phi = np.zeros(Np)
    for ip,pos in enumerate(points):
        x,y,z=pos[0],pos[1],pos[2]
        #
        xy=np.sqrt(x**2+y**2)
        theta[ip] = np.arctan2(y,x)
        phi[ip] = np.arctan2(xy,z);
    return theta,phi
# -----------------------------------------------------------------#
def project_onto_sph(points,radius=1):
    ''' The function take the points (3 dimensional vector) and 
        projects on an unit sphere. '''
    Np = points.shape[0]
    proj_pnts=np.zeros([Np,3])
    for ip,pos in enumerate(points):
        # x,y,z=pos[0],pos[1],pos[2]
        rad_act=la.norm(pos)
        proj_pnts[ip,0]=pos[0]/rad_act  
        proj_pnts[ip,1]=pos[1]/rad_act
        proj_pnts[ip,2]=pos[2]/rad_act
    return proj_pnts
# ----------------------------------------------------------------- #
def height_field_lm(h_theta_phi,theta,phi,area,l,m,radius=1):
    # m=l
    Np = h_theta_phi.shape[0]
    h_lm=0
    # theta=list(map(math.degrees,theta))
    # phi=list(map(math.degrees,phi))
    for ip in range(Np):
        h_lm = h_lm + h_theta_phi[ip]*sph_harm(m,l,theta[ip],phi[ip]).real*area[ip]
    return h_lm
# ----------------------------------------------------------------- #
def spectra(infile,Np=5120,lmax=10):
    # if Np is None:
    #     pname=foldername(infile)+"/para_file.in"
    #     Np=pio.read_param(pname)['Np']
    points,cells = vio.vtk_to_data(infile,Np=Np)
    theta,phi = cart2sph(points)
    proj_pnts = project_onto_sph(points)
    h_theta_phi = la.norm(points,axis=1)-la.norm(proj_pnts,axis=1)
    # h_theta_phi = sph_harm(0,5,theta,phi).real
    # print(np.mean(h_theta_phi))
    sv = SphericalVoronoi(proj_pnts)
    area=sv.calculate_areas()
    h_lm=np.zeros(lmax)
    for l in range(0,lmax):
        h_lm[l]= height_field_lm(h_theta_phi,theta,phi,area,l,m=0)
    # h_lm=np.zeros(lmax*lmax+2*lmax)
    # count=0
    # for l in range(0,lmax):
    #     for m in range(-l,l+1):
    #        h_lm[count]= height_field_lm(h_theta_phi,theta,phi,area,l,m)
    #        count = count+1
        # h_lm[l]=h_lm[l]/(2*l+1)
    # print(count,lmax*lmax+2*lmax)
    return h_lm
# ----------------------------------------------------------------- #
def stat_error(data,nch=10):
    '''
    Takes a one dimensional function and returns error by binning.
    It assumes that the quantity of interest is average of the data.
    It returns 1) the average quantity 2) error in the quantity
    '''
    Nmc=data.shape[0]
    chsz=int(Nmc/nch)
    if (len(data.shape)>1):
        datach=np.zeros([nch,data.shape[1]])
    else:
        datach=np.zeros(nch)
    for ich in range(nch):
        datach[ich]=np.mean(data[ich*chsz:(ich+1)*chsz],axis=0)    
    return np.std(datach,axis=0),np.mean(datach,axis=0)