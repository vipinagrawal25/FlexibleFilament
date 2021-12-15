from pylab import *
import ipywidgets as wid
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.fftpack import dst
from scipy.interpolate import interp1d

def curvatureplot(FILE='output',omega=3,length=1,tmin=0,tmax=-1):
    dd = loadtxt(FILE+'/kappa.txt')
    dd_SS = loadtxt(FILE+'/material_point.txt')
    Y = dd_SS[tmin:tmax]
    nrowcol=dd_SS.shape
    nsnap=nrowcol[0]
    Np = nrowcol[1]
    # Now normalize the material point coordinate
    for isnap in range(0,nsnap-1,1):
        Y[isnap,:] = Y[isnap,:]*length/Y[isnap,-1]

    timeAxis = loadtxt(FILE+'/time.txt')
    timeAxis = timeAxis[tmin:tmax]
    IndexAxis = arange(Np)
    heightAxis = dd[tmin:tmax,1:Np+1]
    
    (IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
    fig = figure()
    ax = fig.add_subplot(111)
    surf = ax.contourf(X,Y,heightAxis)
    cbar = colorbar(surf)
    ax.set_aspect(20)
    ax.set_xlabel('time')
    ax.set_ylabel('material_point')
    plt.show()
    plt.savefig('curvatureplot.eps')
    close()
    
def SineCurvature(FILE,aspect_ratio=1):
    dd = loadtxt(FILE+'/curvature.txt')
    dd_SS = loadtxt(FILE+'/material_point.txt')
    Y = dd_SS
    nrowcol=dd_SS.shape
    nsnap=nrowcol[0]
    Np=nrowcol[1]
    # Now normalize the material point coordinate
    for isnap in range(0,nsnap-1,1):
        Y[isnap,:] = Y[isnap,:]/Y[isnap,-1]

    Z = dd[:,1:Np+1]        #Removed the first entry
    Zsine=zeros((nsnap,Np))
    #
    for isnap in range(0,nsnap):
        Zsine[isnap,:] = dst(Z[isnap,:],type=1)

    timeAxis = loadtxt(FILE+'/time.txt')
    IndexAxis = arange(Np)
    # heightAxis = ndimage.rotate(heightAxis,90)
    # imshow(heightAxis,cmap=plt.cm.gray)
    # colorbar()
    (IndexAxis,X) = meshgrid(IndexAxis, timeAxis)
    fig = figure()
    ax = fig.add_subplot(111)
    # surf = ax.contourf(X,Y,heightAxis,vmin=10**(-5),vmax=10**5,cmap='plasma')
    surf = ax.contourf(X,Y,Zsine)
    cbar = colorbar(surf)
    # ax.set_aspect(aspect_ratio)
    ax.set_xlabel('time')
    ax.set_ylabel('material_point')
    plt.show()
    plt.savefig('sinecurvature.eps')
    close()

def MSD_plot(FILE='',omega=1.5):
    file = loadtxt(FILE+'MSD.txt')
    tt = loadtxt(FILE+'output/time.txt')

    ttsclf = (tt[2]-tt[1])*omega/(2*pi)
    icy = int(1/ttsclf)

    cycle_index = arange(0,tt.size,icy)
    plt.plot(cycle_index,file[cycle_index],'.-') 
    # plt.plot(file)
    # plt.ylim((0, 0.1))
    plt.savefig('MSD_cycle.eps')
    plt.show()

# This function calculates Mean square displacement and substracts the translation of the rod after every cycle 
# if there is any. For this, basically I substract the co-ordinate of middle point.
def MSD_no_trans(FILE='',omega=3):
    tt = loadtxt(FILE+'output/time.txt')
    dd_ini = loadtxt(FILE+'output/var0.txt')
    nrowcol=dd_ini.shape
    Np = nrowcol[0]
    
    sclf = (tt[2]-tt[1])*omega/(2*pi)
    dcy = int(1/sclf)
    cycle_index = arange(0,tt.size,dcy)
    
    MSDnew=zeros(int(tt.size/dcy))
    ndips=int(tt.size/dcy)
    MidP=int(Np/2)

    for idip in range(0,ndips):
        isnap=cycle_index[idip]
        dd=loadtxt(FILE+'output/var'+str(isnap)+'.txt')
        dd[:,0]=dd[:,0]-dd[MidP,0]-dd_ini[:,0]
        dd[:,1]=dd[:,1]-dd[MidP,1]-dd_ini[:,1]
        dd[:,2]=dd[:,2]-dd[MidP,2]-dd_ini[:,2]
        #
        for ip in range(0,Np-1):
            MSDnew[idip]=MSDnew[idip]+dd[ip,0]*dd[ip,0]+dd[ip,1]*dd[ip,1]+dd[ip,2]*dd[ip,2]
            MSDnew[idip]=sqrt(MSDnew[idip])

    # plt.plot(dips,file[dips],'o')
    # plt.legend(['First Line'])
    plt.plot(cycle_index[0:-1],MSDnew,'.-')
    # plt.ylim([1.2, 1.4])
    plt.legend(['MSD cycle','MSD after substracting the translation'])
    plt.savefig('MSD_without_translation.eps')
    close()

def Energy(AA,HH,FILE='output'):
    ddBE=loadtxt(FILE+'/curvature.txt')
    ddSE = loadtxt(FILE+'/material_point.txt')

    nrowcol=ddBE.shape
    nrows=nrowcol[0]
    ncols=nrowcol[1]
    #
    aa = 1./(ncols-2)           # aa is the equlibrium distance
    BE= zeros(nrows)            # Bending energy
    SE=zeros(nrows)             # Stretching energy
    TE=zeros(nrows)             # Total Energy
    time = zeros(nrows)
    #
    for i in range(0,nrows-1):
        time[i] = ddBE[i,0]
        for j in range(1,ncols):
            BE[i] = BE[i]+AA/aa*ddBE[i,j]
            SE[i] = SE[i]+HH/(2*aa)*(ddSE[i,j-1]-aa*(j-1))**2
            # TE[i] = BE[i]+SE[i]
    return time,BE,SE

def LeebyL(Np,Folder='output'):
    ddMP = loadtxt(Folder+'/material_point.txt')
    time = loadtxt(Folder+'/time.txt')
    nrowcol = ddMP.shape
    nsnap=nrowcol[0]

    LeebyL = zeros(nsnap)
    for isnap in range(0,nsnap,1):
        dd = loadtxt(Folder+'/var'+str(isnap)+'.txt')
        Lee = (dd[0,0]-dd[-1,0])**2 + (dd[0,1]-dd[-1,1])**2 + (dd[0,2]-dd[-1,2])**2
        Lee = sqrt(Lee)
        #
        L = 0
        for ip in range(0,Np-1):
            L = L+(dd[ip+1,0]-dd[ip,0])**2 + (dd[ip+1,1]-dd[ip,1])**2 + (dd[ip,2]-dd[ip+1,2])**2
        LeebyL[isnap] = Lee/L

    plt.plot(time,LeebyL,'o-')
    plt.plot(time,Lee,'o-')
    # plt.ylim((0, 0.1))
    plt.savefig('LeebyL.eps')
    plt.show()

    return LeebyL

def PowerSpec(hh,Delta,deali=True):
    hh=hh-NP.mean(hh)
    NN = int(hh.size)
    NNBy2 = int(NN/2)
    # Check whether hh is one dimensional, or power of 2
    ff = (1/(NN*Delta))*arange(0,NNBy2+1,dtype=int)
    HH = fft.rfft(hh)
    Pxx = abs(HH)
    
    if(deali):
        ff = ff[0:int(0.65*ff.size)]        #De-aliasing the data
        Pxx = Pxx[0:int(0.65*Pxx.size)]
        HH = HH[0:int(0.65*HH.size)]
    
    return ff, HH, Pxx
    
def CurvatureSign(kappasqr,yy,zz,eps):
    # kappasqr[kappasqr<th*max(kappasqr)]=0
    NN=yy.size
    firstIn=2
    sign=zeros(NN)
    znew=zeros(NN)
    # First shift everything by the linear line joining the first and last line.
    # yy=yy-yy[0]-(yy[0]-yy[-1])/(zz[0]-zz[-1])*zz
    # Now interpolate things
    for iN in range(0,NN):
        znew[iN] = zz[0]+(zz[-1]-zz[0])*iN/NN
    ff = interp1d(zz,yy)
    ynew=ff(znew)
    # MidP=int(NN/2)
    # ynew=ynew-(ynew[MidP]+(znew-znew[MidP])*(ynew[MidP+1]-ynew[MidP-1])/(znew[MidP+1]-znew[MidP-1]))
    # Assign sign based on the slope in starting. Remember first two points have zero curvature so just go two 
    #more points to make sure slope calculation is right and implement it.
    ynewdiff=diff(yy)
    znewdiff=diff(zz)
    vectors=zeros([NN-1,2])
    vectors[:,0]=ynewdiff
    vectors[:,1]=znewdiff
    vectors=vectors/sqrt(vectors[:,0]**2+vectors[:,1]**2)
    cp=cross(vectors[0:-1],vectors[1:])
    cpdiff=diff(cp)
    cpdiff[cpdiff<0]=-1
    cpdiff[cpdiff>0]=1
    sign[2:NN-1]=cpdiff
    return sign

def GetCurv(Folder='output/',code='CPU',dim=3,wDataMeth=1):
    # This function will take the square root of curvature. Sign of the final thing would be decided
    # by double derivative of the position vector. If the function is convex, curvature can be negative
    # else the curvature would be positive.
    dd=loadtxt(Folder+'curvature.txt')
    time=dd[:,0]
    nrowcol=dd.shape
    nsnap = nrowcol[0]
    print(nsnap)
    NN=nrowcol[1]-1
    kappa=zeros([nsnap,NN+1])
    kappa[:,0]=time
    tangent=zeros([NN-1,2])
    position=zeros([NN,2])
    if (code=='CPU'):
        for isnap in range(1,nsnap):
            dd = loadtxt(Folder+'var'+str(isnap)+'.txt')
            if wDataMeth==1 and dim==3:
                position=dd[:,1:3]
            elif wDataMeth==1 and dim==2:
                position=dd[:,0:2]
            elif wDataMeth==2 and dim==2:
                position[:,0] = dd[2::4]
                position[:,1] = dd[4::4]
            tangent=diff(position,axis=0)
            for iN in range(NN-1):
                tangent[iN,:]=tangent[iN,:]/sqrt( tangent[iN,0]**2 + tangent[iN,1]**2 )
            kappa[isnap,1:NN-1]=cross(tangent[0:-1],tangent[1:])
    elif (code == 'GPU'):
        dd=loadtxt(Folder+"PSI")
        zz=zeros(NN)
        yy=zeros(NN)
        for isnap in range(1,nsnap):
            for iN in range(0,NN):
                yy[iN] = dd[isnap,3*iN+1]
                zz[iN] = dd[isnap,3*iN+3]
            tangent[:,0]=diff(yy)
            tangent[:,1]=diff(zz)
            for iN in range(NN-1):
                tangent[iN,:]=tangent[iN,:]/sqrt( tangent[iN,0]**2 + tangent[iN,1]**2 )
            kappa[isnap,1:NN-1]=cross(tangent[0:-1],tangent[1:])            
            # kappa[isnap,1:NN+1]=sqrt(curvsqr[isnap,:])*CurvatureSign(curvsqr[isnap,:],yy,zz,eps)
    savetxt(Folder+'kappa.txt',kappa,fmt='%.5e')

def GetCurv2(Folder='output/',code='CPU',dim=3,wDataMeth=1):
    # This function will take the square root of curvature. Sign of the final thing would be decided
    # by double derivative of the position vector. If the function is convex, curvature can be negative
    # else the curvature would be positive.
    dd=loadtxt('revolve.txt')
    # time=dd[:,0]
    nsnap = int(dd[1])
    print(nsnap)
    NN=256
    kappa=zeros([nsnap,NN+1])
    # kappa[:,0]=time
    tangent=zeros([NN-1,2])
    position=zeros([NN,2])
    if (code=='CPU'):
        for isnap in range(1,nsnap):
            kappa[isnap,0]=nsnap/50;
            dd = loadtxt(Folder+'var'+str(isnap)+'.txt')
            if wDataMeth==1 and dim==3:
                position=dd[:,1:3]
            elif wDataMeth==1 and dim==2:
                position=dd[:,0:2]
            elif wDataMeth==2 and dim==2:
                position[:,0] = dd[2::4]
                position[:,1] = dd[4::4]
            tangent=diff(position,axis=0)
            for iN in range(NN-1):
                tangent[iN,:]=tangent[iN,:]/sqrt( tangent[iN,0]**2 + tangent[iN,1]**2 )
            kappa[isnap,1:NN-1]=cross(tangent[0:-1],tangent[1:])
    elif (code == 'GPU'):
        dd=loadtxt(Folder+"PSI")
        zz=zeros(NN)
        yy=zeros(NN)
        for isnap in range(1,nsnap):
            for iN in range(0,NN):
                yy[iN] = dd[isnap,3*iN+1]
                zz[iN] = dd[isnap,3*iN+3]
            tangent[:,0]=diff(yy)
            tangent[:,1]=diff(zz)
            for iN in range(NN-1):
                tangent[iN,:]=tangent[iN,:]/sqrt( tangent[iN,0]**2 + tangent[iN,1]**2 )
            kappa[isnap,1:NN-1]=cross(tangent[0:-1],tangent[1:])            
            # kappa[isnap,1:NN+1]=sqrt(curvsqr[isnap,:])*CurvatureSign(curvsqr[isnap,:],yy,zz,eps)
    savetxt(Folder+'kappa.txt',kappa,fmt='%.5e')

def SineTransform(dd,length=1.28,dia=0.005):
    NN=dd.size
    CurvSine=1/NN*dst(dd,type=1)
    wavenrs = length/(dia*NN)*linspace(0,NN,NN)
    return wavenrs,CurvSine

def vel_tracer(dd,vel_abs,time,Xtracer,wave='sine',height=1.28,ShearRate=2,sigma=1.5,dia=0.005):
    # Instant velocity is accurate upto the order of h^2. We can have different formulae if want more
    # accurate results
    NN = dd[:,0].size
    vel_amb = zeros([NN,3])
    tracer_amb = zeros(3)
    omega=ShearRate*sigma
    if wave=='sine':
        for i in range(0,NN):
            vel_amb[i,1] = ShearRate*(height-dd[i,2])*sin(omega*time)
        tracer_amb[1] = ShearRate*(height-Xtracer[2])*sin(omega*time)
    else:
        print("What type of shear rate is given? Specify that right")
        raise ValueError
    vel_ins = vel_abs - vel_amb
    
    dij = [[1,0,0],[0,1,0],[0,0,1]]   # Unit second rank tensor
    Vtracer = zeros((1,3))
    GG = zeros((3,3))
    
    for i in range(0,NN):
        Xi = Xtracer-dd[i,:] 
        XiXj = outer(Xi,Xi)                     # XiXj second order tensor
        rr = linalg.norm(Xi)
        GG = (dij/rr + XiXj/(rr**3))            # Green's function for the current point
        del2G = (dij/(rr**3) - 3*XiXj/(rr**5))  # Double derivative of Green's function
        Vtracer =  Vtracer + 3*dia/8*tensordot(vel_ins[i,:],GG,axes=1) + (dia**3)/32*tensordot(vel_ins[i,:],del2G,axes=1)
    return Vtracer

def VtraceTimeSer(Folder='output/',Xtracer=[0.01,0,0.5],code='CPU',sigma=1.5):
    if code=='CPU':
        dd= loadtxt(Folder+'var1.txt')
        height=max(dd[2,:])
        nrowcol=dd.shape
        NN=nrowcol[0]
        time = loadtxt(Folder+'time.txt')
        ndiag = time.size
        Tmax = time[-1]
        position = zeros([NN,3])
        vel_abs = zeros([NN,3])
        Vtracer = zeros([ndiag,3])
        for fnumber in range(1,ndiag):
            print(fnumber)
            file = loadtxt(Folder+'var'+str(fnumber)+'.txt')
            position=file[:,0:3]
            vel_abs=file[:,3:6]
            Vtracer[fnumber-1,:]=vel_tracer(position,vel_abs,time[fnumber],Xtracer,height=height,sigma=sigma)
    elif code=='GPU':
        dd_pos=loadtxt(Folder+'PSI')
        dd_vel=loadtxt(Folder+'VEL')
        height=max(dd_pos[0,:])
        nrowcol=dd_pos.shape
        NN=int((nrowcol[1]-1)/3)
        time=dd_pos[:,0]
        ndiag=time.size
        Tmax=time[-1]
        position=zeros([NN,3])
        vel_abs=zeros([NN,3])
        Vtracer = zeros([ndiag,3])
        for fnumber in range(1,ndiag):
            print(fnumber)
            for iN in range(NN):
                position[iN,0]=dd_pos[fnumber,3*iN+1]
                position[iN,1]=dd_pos[fnumber,3*iN+2]
                position[iN,2]=dd_pos[fnumber,3*iN+3]
                vel_abs[iN,0]=dd_vel[fnumber,3*iN+1]
                vel_abs[iN,1]=dd_vel[fnumber,3*iN+2]
                vel_abs[iN,2]=dd_vel[fnumber,3*iN+3]
            Vtracer[fnumber-1,:]=vel_tracer(position,vel_abs,time[fnumber],Xtracer,height=height,sigma=sigma)
    print("I am done")
    return time, Vtracer

def reScale(file,factor,file2="downscale.txt",wDataMeth=1):
    dd = loadtxt(file);
    Nrows,Ncols = dd.shape
    dd2 = dd[0:factor:-1,:]
    return dd,dd2
    
def getpos(dd,dim=2,wDataMeth=1):
    NN = dd.size/(dim*2)
    position=NP.zeros([NN,2])
    if wDataMeth==1 and dim==3:
        position=dd[:,1:3]
    elif wDataMeth==1 and dim==2:
        position=dd[:,0:2]
    elif wDataMeth==2 and dim==2:
        position[:,0] = dd[2::4]
        position[:,1] = dd[4::4]
    elif wDataMeth==2 and dim==3:
        position[:,0] = dd[4::6] 
        position[:,1] = dd[6::6]
    else:
        print("I have no idea what do you want?")
    return position