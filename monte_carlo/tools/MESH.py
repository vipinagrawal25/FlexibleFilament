import numpy as np
import numpy.linalg as LA
import quaternion
#################################### CLASS ###################################
class Mesh:
    cmlst=None
    node_nbr=None
    bond_nbr=None
    radius=1.0
    sp_curv=2/radius
    #
    def __init__(self,R,cells):
            self.Np=R.shape[0]
            self.R=R
            self.cells=cells
            # cmlst,node_nbr,bond_nbr=neighbours(R,cells)
    #
    def __neighbours(self):
        simpl=self.__sort_simplices()
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
    def __sort_simplices(self):
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
        cmlst,node_nbr,bond_nbr=self.__neighbours()
        self.cmlst=np.array(cmlst).astype(np.int)
        self.node_nbr=np.array(node_nbr).astype(np.int)
        self.bond_nbr = np.array(bond_nbr).astype(tuple)
        self.__sort_nbrs()
    #
    def curvature(self):
        '''Given a snapshot, the function computes curvature for all points.'''
        Np = self.Np
        curv = np.zeros([Np,3])
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
            curv[i]=(1/sigma_i)*(cot_times_rij)
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
    #
    def __sort_nbrs(self):
        zhat = np.array([0.,0.,1.])
        for i in range(self.Np):
            nbrs=self.node_nbr[self.cmlst[i]:self.cmlst[i+1]]  # neighbours of ith node
            vector=self.R[i]
            # I will rotate the coordinate system about this vector
            vhat = np.cross(vector,zhat)       
            vnorm = LA.norm(vhat)
            # If the vector is already lying at z-axis then there is no need to rotate.
            if vnorm>1e-16:
                vhat = vhat/vnorm
                theta = polar(vector)
                # Rotate all the neighbours of a point.
                rotated=rotate(self.R[nbrs],vhat,theta)
                # Since all the voronoi cells are rotated, sort them in anticlockwise direction
                sorted_indices = sort_2Dpoints_theta(rotated[:,0],rotated[:,1])[0]
                self.node_nbr[self.cmlst[i]:self.cmlst[i+1]]=nbrs[sorted_indices]
                bond_nbrs = self.bond_nbr[self.cmlst[i]:self.cmlst[i+1]]
                self.bond_nbr[self.cmlst[i]:self.cmlst[i+1]]=bond_nbrs[sorted_indices]
################################# OTHER FUNCTIONS ###############################
def sort_2Dpoints_theta(x,y):
    len_x = len(x)
    len_y = len(y)
    if len_x!=len_y:
        raise Exception("")
    #
    xsort=np.zeros(len_x)
    ysort=np.zeros(len_y)
    #

    theta=np.arctan2(x,y)+np.pi
    indices=np.linspace(0,len_x-1,len_x)
    xyth=np.transpose(np.array([x,y,theta,indices]))
    #
    xysort = np.asarray(sorted(xyth, key=lambda x: (x[2])))
    return xysort[:,3].astype(int),np.array([xysort[:,0],xysort[:,1]])
#-----------------------------------------------------------------------------#
def polar(xyz):
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    XsqPlusYsq = x**2 + y**2
    return np.arctan2(np.sqrt(XsqPlusYsq),z)
#----------------------------------------------------------------------------#
def rotate(vector,nhat,theta):
    '''rotate a vector about nhat by angle theta'''
    cos_thby2=np.cos(theta/2)
    sin_thby2=np.sin(theta/2)
    q=np.quaternion(cos_thby2,nhat[0]*sin_thby2,nhat[1]*sin_thby2,nhat[2]*sin_thby2)
    q_inv=np.quaternion(cos_thby2,-nhat[0]*sin_thby2,-nhat[1]*sin_thby2,-nhat[2]*sin_thby2)
    nn=vector.shape[0]
    rot_vec=np.zeros([nn,3])
    for i in range(nn):
        q_vec=np.quaternion(0,vector[i][0],vector[i][1],vector[i][2])
        rot_vec[i]=quater2vec(q*q_vec*q_inv)
    return rot_vec
#----------------------------------------------------------------------------#
def quater2vec(qq,precision=1e-16):
    if qq.w>1e-8:
        print("# ERROR: Quaternion has non-zero scalar value.\n \
               # Can not convert to vector.")
        exit(1)
    return np.array([qq.x,qq.y,qq.z])
#----------------------------------------------------------------------------#