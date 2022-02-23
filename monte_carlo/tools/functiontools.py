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
class MESH:
    def __init__(self,R,cells):
            self.Np=R.shape[0]
            self.R=R
            self.cells=cells
            self.cmlst=None
            self.node_nbr=None
            self.bond_nbr=None
            self.radius=1.0
            self.sp_curv=2/self.radius
            # cmlst,node_nbr,bond_nbr=neighbours(R,cells)
    #
    def neighbours(self):
        simpl=self.sort_simplices()
        Np=self.Np
        r1=simpl[:,0]
        r2=simpl[:,1]
        r3=simpl[:,2]
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
    #
    def sort_simplices(self):
        cells=self.cells
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
    #
    def assign_nbrs(self):
        cmlst,node_nbr,bond_nbr=self.neighbours()
        self.cmlst=np.array(cmlst).astype(np.int)
        self.node_nbr=np.array(node_nbr).astype(np.int)
        self.bond_nbr = np.array(bond_nbr).astype(tuple)
    #
    def curvature(self):
        '''Given a snapshot, the function computes curvature for all points.'''
        Np = self.Np
        curv = np.zeros(Np)
        R=self.R
        if self.cmlst is None or self.node_nbr is None or self.bond_nbr is None:
            print("# Seems like you are calling the function for the first time, as we have not \
                     constructed the neighbour list.\n Dont worry, I am going to do it now.")
            self.assign_nbrs()
            print("# SUCCESSS:neighbour list is succesfully constructed.")
        # start=self.cmlst[i]
        # end=self.cmlst[i+1]
        # Initialize things
        for i in range(Np):
            start=self.cmlst[i]
            end=self.cmlst[i+1]
            count=0
            cot_times_rij=np.zeros(3)
            sigma_i=0
            for j in self.node_nbr[start:end]:
                k=int(self.bond_nbr[start+count][0])
                kp=int(self.bond_nbr[start+count][1])
                # Compute all the position difference 
                rij=R[i]-R[j]
                rik=R[i]-R[k]
                rjk=R[j]-R[k]
                rikp=R[i]-R[kp]
                rjkp=R[j]-R[kp]
                # compute all the length for both triangles.
                lijsq=np.inner(rij,rij)
                liksq=np.inner(rik,rik)
                ljksq=np.inner(rjk,rjk)
                likpsq=np.inner(rikp,rikp)
                ljkpsq=np.inner(rjkp,rjkp)
                # compute area of both triangles
                area_ijk=0.5*LA.norm(np.cross(rij,rjk))
                area_ijkp=0.5*LA.norm(np.cross(rij,rjkp))
                # compute all the cot angles.
                cot_jk = 0.25*(lijsq+ljksq-liksq)/area_ijk;
                cot_k = 0.25*(ljksq+liksq-lijsq)/area_ijk;
                cot_jkp = 0.25*(lijsq+ljkpsq-likpsq)/area_ijkp;
                cot_kp =  0.25*(ljkpsq+likpsq-lijsq)/area_ijkp;
                cot_sum=0.5*(cot_k+cot_kp)
                #
                cot_times_rij=cot_times_rij+cot_sum*rij
                #
                count=count+1
                sigma_i=sigma_i+voronoi_area(cot_jk,cot_k,liksq,lijsq,area_ijk)
                sigma_i=sigma_i+voronoi_area(cot_jkp,cot_kp,likpsq,lijsq,area_ijkp);
            sigma_i=sigma_i/2
            curv[i]=(1/sigma_i)*LA.norm(cot_times_rij)
        return curv
    #
    def compute_obtuse(self):
        """For a given triangular mesh, the function returns the list of all the obtuse triangles.
        Takes the triangles and points as input and returns an array of shape number of triangles 
        #
        The output is 1 if one of the triangle is greater than
        90, and 0 if all the angle is acute (less than 90)
        """
        # angle theta for triangle ijk
        cells=self.cells
        Np=self.Np
        points=self.R
        isgt90 = np.zeros(np.shape(cells)[0],dtype=np.float64)
        piby2 = 0.5*np.pi
        a_ij_ik=np.zeros(Np)
        a_ji_jk=np.zeros(Np)
        a_ki_kj=np.zeros(Np)
        for i, tris in enumerate(cells):
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
    import numpy as np
#-----------------------------------------
def write_param(fname="../para_file.in",paramdict=None):
    if paramdict is None:
        paramdict={"N":5120,"coef_bending":2.5,"coef_stretching":25,"coef_vol_expansion":1e6,
                "radius":1.0,"pos_bot_wall":-1.05,"sigma":0.17,"epsilon":4,
                "Dfac":4,"kBT":1,"is_restart":1,"mc_total_iters":10000,"mc_dump_iter":10,
                "afm_N":0,"extent_l":-0.1,"extent_r":0.1,"extent_t":-0.1,"extent_b":0.1,
                "tip_radius":0.2,"tip_pos_z":1.05,"afm_sigma":0.17,"afm_epsilon":4}
    with open(fname, "w") as file:
        file.write("## Membrane parameters\n")
        file.write("N\tcoef_bending\tcoef_stretching\tcoef_vol_expansion\n")
        file.write("%d %3.2f %3.2f %3.2f\n" %(paramdict['N'],paramdict['coef_bending'],
                    paramdict['coef_stretching'],paramdict['coef_vol_expansion']))
        file.write("radius\tpos_bot_wall\tsigma\tepsilon\n")
        file.write("%3.2f %3.2f %3.2f %3.2f\n" %(paramdict['radius'],paramdict['pos_bot_wall'],
                    paramdict['sigma'],paramdict['epsilon']))
        file.write("## Montecarlo parameters\n")
        file.write("Dfac\tkBT is_restart mc_total_iters\tmc_dump_iter\n")
        file.write("%d %3.2f %d %d %d\n" %(paramdict['Dfac'],paramdict['kBT'],paramdict['is_restart'],
                    paramdict['mc_total_iters'],paramdict['mc_dump_iter']))
        file.write("## Afm Tip parameters\n")
        file.write("N\textent_l\textent_r\textent_t\textent_b\n")
        file.write("%d %1.1f %1.1f %1.1f %1.1f\n" %(paramdict['afm_N'],paramdict['extent_l'],
                    paramdict['extent_r'],paramdict['extent_t'],paramdict['extent_b']))
        file.write("tip_radius\ttip_pos_z\tafm_sigma\tafm_epsilon\n")
        file.write("%2.2f %2.2f %2.2f %2.2f" %(paramdict['tip_radius'],paramdict['tip_pos_z'],
                    paramdict['afm_sigma'],paramdict['afm_epsilon']))
#-----------------------------------------
def read_param(fname='../para_file.in'):
    ## Membrane parameters
    N,coef_bending,coef_stretching,coef_vol_expansion = np.loadtxt(fname,skiprows=2,max_rows=1)
    N=int(N)
    radius,pos_bot_wall,sigma,epsilon = np.loadtxt(fname,skiprows=4,max_rows=1)
    ## Montecarlo parameters
    Dfac,kBT,is_restart,mc_total_iters,mc_dump_iter = np.loadtxt(fname,skiprows=7,max_rows=1)
    Dfac=int(Dfac)
    is_restart=int(is_restart)
    mc_total_iters=int(mc_total_iters)
    mc_dump_iter=int(mc_dump_iter)
    ## Afm Tip parameters
    afm_N,extent_l,extent_r,extent_t,extent_b = np.loadtxt(fname,skiprows=10,max_rows=1)
    afm_N=int(afm_N)
    tip_radius,tip_pos_z,afm_sigma,afm_epsilon = np.loadtxt(fname,skiprows=12,max_rows=1)
    paramdict={"N":N,"coef_bending":coef_bending,"coef_stretching":coef_stretching,
                "coef_vol_expansion":coef_vol_expansion,
                "radius":radius,"pos_bot_wall":pos_bot_wall,"sigma":sigma,"epsilon":epsilon,
                "Dfac":Dfac,"kBT":kBT,"is_restart":is_restart,"mc_total_iters":mc_total_iters,
                "mc_dump_iter":mc_dump_iter,
                "afm_N":afm_N,"extent_l":extent_l,"extent_r":extent_r,"extent_t":extent_t,
                "extent_b":extent_b,
                "tip_radius":tip_radius,"tip_pos_z":tip_pos_z,"afm_sigma":afm_sigma,
                "afm_epsilon":afm_epsilon}
    return paramdict
#-------------------------------------------------------------
def change_param(fname="../para_file.in",**kwargs):
    ''' The function change one parameter from the input file and overwrites the new parameter file.'''
    paramdict=read_param(fname=fname)
    for key,value in kwargs.items():
        paramdict[key]=value
    write_param(fname=fname,paramdict=paramdict)
#-------------------------------------------------------------