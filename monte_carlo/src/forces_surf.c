#include "../include/global.h"
#include "../include/subroutine.h"

#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

 double inner_product(POSITION s1, POSITION s2){
    return s1.x*s2.x + s1.y*s2.y + s1.z*s2.z;
}

POSITION Position_add(POSITION s1, 
        POSITION s2, double fac){
    /* Returns s1 + fac*s2*/

    POSITION add;

    add.x = s1.x + fac*s2.x;
    add.y = s1.y + fac*s2.y;
    add.z = s1.z + fac*s2.z;

    return add;

}


POSITION cross_product(POSITION s1, 
        POSITION s2){

    POSITION crosprod;
    crosprod.x = s1.y*s2.z - s1.z*s2.y;
    crosprod.y = s1.z*s2.x - s1.x*s2.z;
    crosprod.z = s1.x*s2.y - s1.y*s2.x;
    return crosprod;

}
double cotangent(POSITION si, POSITION sk, POSITION sj){
    POSITION drik, drjk, cross;
    double cot_theta, inp_1;  
    double inner_prod;
    //
    drik = Position_add(si , sk, -1e0);
    drjk = Position_add(sj , sk, -1e0);
    cross = cross_product(drik, drjk);
    inner_prod = inner_product(drik, drjk);
    cot_theta = inner_prod/sqrt(inner_product(cross,cross));
    //
    return cot_theta;
}
//
double cotangent(double a, double b, double c){
    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double cot_theta=0.25*(a*a+b*b-c*c)/area;
    return cot_theta;
}
double volume_ipart(POSITION *pos, 
        int *node_nbr, int2* bond_nbr,
        int num_nbr, int idx, MBRANE_para para){

    int i, j, k, kp, it, itn;
    double volume1, volume2;
    POSITION area1, area2, rk, ri, rj, rkp;
    POSITION rij, rijk, rik, rjkp;
    POSITION rijkp;
    POSITION rjk, pt; 
    double dir_norm, ini_vol;
    double inp_r1, inp_r2, inp_r3;
    FILE *fid;

    volume1 = 0e0;
    volume2 = 0e0;
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k  = bond_nbr[i].i1; 
        kp = bond_nbr[i].i2; 
        ri = pos[idx]; rj = pos[j]; 
        rk = pos[k]; rkp = pos[kp];

        rij = Position_add(ri, rj, 1e0);
        rijk = Position_add(rij, rk, 1e0);
        rijkp = Position_add(rij, rkp, 1e0);

        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;
        rijkp.x = rijkp.x/3e0; 
        rijkp.y = rijkp.y/3e0; rijkp.z = rijkp.z/3e0;

        rij  = Position_add(rj, ri, -1e0);
        rjk  = Position_add(rj , rk, -1e0);
        rjkp = Position_add(rj , rkp, -1e0);


        area1 = cross_product(rij, rjk);
        area2 = cross_product(rij, rjkp);
        double ip1 = inner_product(area1,rijk);
        double ip2 = inner_product(area2,rijkp);

        if(sign(ip1) == 1) volume1 = volume1 + ip1;
        if(sign(ip1) == -1) volume2 = volume2 + ip1;

        if(sign(ip2) == 1) volume1 = volume1 + ip2;
        if(sign(ip2) == -1) volume2 = volume2 + ip2;
    }
    volume1 = volume1/3e0;
    return volume1;

}

double stretch_energy_ipart(POSITION *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para){
    //
    double HH;
    double idx_ener;
    POSITION rij;
    double mod_rij;
    int i,j;
    //
    idx_ener = 0e0;
    HH = para.coef_str/(para.av_bond_len*para.av_bond_len);
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        rij = Position_add(pos[idx], pos[j], -1e0);
        mod_rij = sqrt(inner_product(rij, rij));
        //
        idx_ener = idx_ener + (mod_rij - lij_t0[i])*(mod_rij - lij_t0[i]);
    }
    //
    return 0.5*idx_ener*HH; 
}
//
double norm(POSITION rr){
    return sqrt(inner_product(rr,rr));
}
//
double bending_energy_ipart(POSITION *pos, int *node_nbr, int2 *bond_nbr, int num_nbr,
                            int idx, MBRANE_para para, string method){
    double bend_ener,curvature,sigma_i;
    POSITION cot_times_rij;
    double BB=para.coef_bend;
    double curv_t0 = 2e0/para.radius;
    cot_times_rij.x = 0e0;
    cot_times_rij.y = 0e0;
    cot_times_rij.z = 0e0;
    //
    if (method=="old"){
        double idx_ener;
        double rij_dot_rij;
        int i, k, kp;
        int j;
        double cot_sum, curv_t0;
        // double curvature;
        POSITION cot_theta_rij, rij, rik, rikp;
        //
        sigma_i = 0e0;
        //
        // for (i =0; i < num_nbr; i++){
        //     j = node_nbr[i];
        //     k  = bond_nbr[i].i1; 
        //     kp = bond_nbr[i].i2;
        //     //
        //     cot_sum=0.5*(cotangent(pos[idx],pos[k],pos[j]) +
        //                  cotangent(pos[idx],pos[kp],pos[j]));
        //     rij = Position_add(pos[idx], pos[j], -1e0);
        //     // rik = Position_add(pos[idx], pos[k], -1e0);
        //     // rikp = Position_add(pos[idx], pos[kp], -1e0);
        //     // area = area+sqrt(inner_product(cross_product(rij,rik),cross_product(rij,rik)))+
        //     //             sqrt(inner_product(cross_product(rij,rikp),cross_product(rij,rikp)));
        //     rij_dot_rij = inner_product(rij, rij);
        //     cot_times_rij  = Position_add(cot_times_rij, rij, cot_sum);
        //     // cot_times_rij  = Position_add(cot_times_rij, rij, 1/sqrt(3));
        //     sigma_i = sigma_i + rij_dot_rij*cot_sum;
        // }
        // sigma_i *= 0.25;        
        int jdx,kdx,kpdx;
        POSITION xij,xik,xjk,xikp,xjkp;
        double cot_jdx_k,cot_jdx_kp,cot_kdx,cot_kpdx;
        double lijsq,liksq,ljksq,likpsq,ljkpsq;
        for (i =0; i < num_nbr; i++){
            jdx = node_nbr[i];
            kdx  = bond_nbr[i].i1; 
            kpdx = bond_nbr[i].i2;
            //
            cot_jdx_k = cotangent(pos[idx],pos[jdx],pos[kdx]);
            cot_kdx = cotangent(pos[idx],pos[kdx],pos[jdx]);
            cot_jdx_kp = cotangent(pos[idx],pos[jdx],pos[kpdx]);
            cot_kpdx = cotangent(pos[idx],pos[kpdx],pos[jdx]);
            //
            cot_sum=0.5*(cotangent(pos[idx],pos[kdx],pos[jdx]) +
                         cotangent(pos[idx],pos[kpdx],pos[jdx]));
            //
            xij = Position_add(pos[idx], pos[jdx], -1e0);
            xik = Position_add(pos[idx], pos[kdx], -1e0);
            xjk = Position_add(pos[jdx], pos[kdx], -1e0);
            xikp = Position_add(pos[idx], pos[kpdx], -1e0);
            xjkp = Position_add(pos[jdx], pos[kpdx], -1e0);
            //
            lijsq = inner_product(xij,xij);
            liksq = inner_product(xik,xik);
            ljksq = inner_product(xjk,xjk);
            likpsq = inner_product(xikp,xikp);
            ljkpsq = inner_product(xjkp,xjkp);
            //
            // rik = Position_add(pos[idx], pos[k], -1e0);
            // rikp = Position_add(pos[idx], pos[kp], -1e0);
            // area = area+sqrt(inner_product(cross_product(rij,rik),cross_product(rij,rik)))+
            //             sqrt(inner_product(cross_product(rij,rikp),cross_product(rij,rikp)));
            rij_dot_rij = inner_product(xij, xij);
            cot_times_rij  = Position_add(cot_times_rij, xij, cot_sum);
            // cot_times_rij  = Position_add(cot_times_rij, rij, 1/sqrt(3));
            sigma_i = sigma_i + 0.125*(cot_jdx_k*liksq+
                                       cot_kdx*lijsq+
                                       cot_jdx_kp*likpsq+
                                       cot_kpdx*lijsq);
        }
        sigma_i *= 0.5;
        // is_attractive[idx]=sigma_i<0;
        // if (sigma_i<0){
        // sigma_i=area/3;
        // cout << area/3 << endl;
        // }
    }else{
        double cot_jdx_k,cot_jdx_kp,cot_kdx,cot_kpdx,cot_idx_k,cot_idx_kp;
        double area_ijk,area_ijkp;
        double lijsq,liksq,ljksq,likpsq,ljkpsq;
        // POSITION cot_times_rij;
        POSITION xij,xik,xjk,xikp,xjkp;
        int jdx,kdx,kpdx;
        double cot_sum;
        sigma_i = 0e0;
        for (int j = 0; j < num_nbr; j++){
            jdx = node_nbr[j];
            kdx  = bond_nbr[j].i1;
            kpdx = bond_nbr[j].i2;
            //
            xij = Position_add(pos[idx], pos[jdx], -1e0);
            xik = Position_add(pos[idx], pos[kdx], -1e0);
            xjk = Position_add(pos[jdx], pos[kdx], -1e0);
            xikp = Position_add(pos[idx], pos[kpdx], -1e0);
            xjkp = Position_add(pos[jdx], pos[kpdx], -1e0);
            //
            lijsq = inner_product(xij,xij);
            liksq = inner_product(xik,xik);
            ljksq = inner_product(xjk,xjk);
            likpsq = inner_product(xikp,xikp);
            ljkpsq = inner_product(xjkp,xjkp);
            //
            area_ijk = 0.5*norm(cross_product(xij,xjk)); 
            area_ijkp = 0.5*norm(cross_product(xij,xjkp));
            //
            cot_jdx_k = 0.25*(lijsq+ljksq-liksq)/area_ijk;
            cot_kdx = 0.25*(ljksq+liksq-lijsq)/area_ijk;
            cot_jdx_kp = 0.25*(lijsq+ljkpsq-likpsq)/area_ijkp;
            cot_kpdx =  0.25*(ljkpsq+likpsq-lijsq)/area_ijkp;
            cot_idx_k = 0.25*(lijsq+likpsq-ljkpsq)/area_ijk;
            cot_idx_kp = 0.25*(lijsq+likpsq-ljkpsq)/area_ijkp;
            // compute sigma_i -> first check whether all angles are acute?
            //
            // cot_jdx_k = cotangent(pos[idx],pos[jdx],pos[kdx]);
            // cot_kdx = cotangent(pos[idx],pos[kdx],pos[jdx]);
            // cot_jdx_kp = cotangent(pos[idx],pos[jdx],pos[kpdx]);
            // cot_kpdx = cotangent(pos[idx],pos[kpdx],pos[jdx]);
            // //
            cot_sum = 0.5*(cot_kdx+cot_kpdx);
            cot_times_rij = Position_add(cot_times_rij, xij, cot_sum);
            // sigma_i = sigma_i + 0.125*( cot_jdx_k*liksq+
            //                             cot_kdx*lijsq+
            //                             cot_jdx_kp*likpsq+
            //                             cot_kpdx*lijsq);
            // // sigma_i=sigma_i+0.125*();
            /* cout << cot_jdx_k << endl; */
            /* cout << cot_jdx_kp << endl; */
            if (cot_jdx_k>0 && cot_kdx>0 && cot_idx_k>0){
                // add the voronoi diagram
                sigma_i=sigma_i+0.125*(cot_jdx_k*liksq+cot_kdx*lijsq);
            }else{
                if(cot_idx_k < 0){
                    sigma_i = sigma_i+0.5*area_ijk;
                }else{
                    sigma_i = sigma_i+0.25*area_ijk;
                }
            }
            // similar thing goes for the other triangle.
            if (cot_jdx_kp>0 && cot_kpdx>0 && cot_idx_kp>0){
                // add the voronoi diagram                
                sigma_i=sigma_i+0.125*(cot_jdx_kp*likpsq+cot_kpdx*lijsq);
            }else{
                // add the voronoi diagram                
                // cout << "coming here?" << endl;
                if(cot_idx_kp < 0){
                    sigma_i = sigma_i+0.5*area_ijkp;
                }else{
                    sigma_i = sigma_i+0.25*area_ijkp;
                }
            }            
        }
        sigma_i = 0.5*sigma_i;  // as everything is counted twice.
        // double xij[num_nbr],cot_alpha[num_nbr],cot_beta[num_nbr],xij_sq[num_nbr];
        // int jdx,kdx,jpdx,kpdx,jp;
        // double sigma_i=0;
        // // bool obtuse=0;      // if any cotangent is -ve or cot(alpha+beta)<0
        // for (int j = 0; j < num_nbr; ++j){
        //     // jp = (j+1)%(num_nbr);         // things are cyclic
        //     if(cot_alpha[jp] > 0 && cot_beta[j] > 0 && cot_alpha[jp]*cot_beta[j] < 1){
        //         // this means that all the angles are acute.
        //         // add the voronoi area.
        //         sigma_i=sigma_i+0.25*(cot_alpha[jp]*xij_sq[jp] + cot_beta[j]+xij_sq[j]);
        //     }else{ 
        //         // add non-voronoi area
        //         if(cot_alpha[jp]*cot_beta[j] > 1){
        //         }
        //     }
        // }
    }
    curvature = (1e0/sigma_i)*norm(cot_times_rij);
    // cout << curvature << endl;
    bend_ener = 0.5*BB*sigma_i*(curvature - curv_t0)*(curvature - curv_t0);
    // bend_ener = 0.5*BB*1/area*inner_product(cot_times_rij, cot_times_rij);
    return bend_ener;
 }

double bending_energy_ipart_neighbour(POSITION *pos, 
        MESH mesh, int idx, MBRANE_para para){
/*    Evaluates the bending energy of all the neighbours*/ 
/*    of idx particle*/
   int j, k;
   int num_nbr_j, cm_idx;
   int nbr, cm_idx_nbr;
   double be;

   be = 0e0;

   for (j = mesh.cmlist[idx]; j < mesh.cmlist[idx+1]; j++){
       nbr = mesh.node_nbr_list[j];
       num_nbr_j = mesh.cmlist[nbr+1] -  mesh.cmlist[nbr];
       cm_idx_nbr = mesh.cmlist[nbr];

       be += bending_energy_ipart(pos, 
              (int *) mesh.node_nbr_list + cm_idx_nbr,
              (int2 *) mesh.bond_nbr_list + cm_idx_nbr, num_nbr_j, 
           nbr, para);
   }
   return be;
} 

double bending_energy_total(POSITION *pos, 
        MESH mesh, MBRANE_para para){
    int idx, j, k;
    int num_nbr, cm_idx;
    double be, curv;

    be = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];


        /* fprintf(stderr, "check neighbour of second: %d %d %d \n", num_nbr, cm_idx); */

        be += bending_energy_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx), num_nbr, 
                idx, para);
    }
    return be;
}

double stretch_energy_total(POSITION *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para){
    int idx, j, k;
    int num_nbr, cm_idx;
    double se;

    se = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];

        se += stretch_energy_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (double *) (lij_t0 + cm_idx), num_nbr, 
                idx, para);
        /* printf( "stretch: %lf \n", se); */
    }
    return se*0.5e0;
 }

double lj(double sqdr, double eps){
    double r6;
    r6 = sqdr*sqdr*sqdr;
    return eps*(r6*(r6-1));
}

double lj_bottom_surface(double zz, 
        bool is_attractive, 
        double sur_pos, double eps, double sigma){

    double inv_sqdz, ds;

    if(is_attractive){
        ds = sur_pos - zz;
        inv_sqdz = (sigma*sigma)/(ds*ds);
    }
    return  lj(inv_sqdz, eps);
}



void identify_attractive_part(POSITION *pos, 
        bool *is_attractive, int N){

    int i; 
    double theta, rr;
    for(i= 0; i<N; i++){
        theta = pi - acos(pos[i].z);
        is_attractive[i] = theta < pi/6.0;
    }
}

double rep_lj_afm(POSITION pos, AFM_para afm){
    int i;
    double ener_afm, ds;
    double ds_sig_inv, r6;
    POSITION dr;

    ener_afm = 0e0;
    if(fabs(afm.tip_pos_z - pos.z) < 4*afm.sigma) {
        for(i = 0; i<afm.N; i++) {
            dr = Position_add(pos, afm.tip_curve[i], -1); 
            ds = (inner_product(dr,dr));
            ds_sig_inv = (afm.sigma*afm.sigma)/ds;
            r6 = ds_sig_inv*ds_sig_inv*ds_sig_inv;
            ener_afm = ener_afm +  afm.epsilon*r6*r6;
        }
    }
   return ener_afm;

};
double volume_enclosed_membrane(POSITION *pos, 
        int *triangles, int num_triangles){
    int i, j, k, it;
    double volume;
    POSITION area, rk, ri, rj;
    POSITION rij, rijk, rik;
    POSITION rjk, pt; 
    double dir_norm;
    double inp_r1, inp_r2, inp_r3;
    FILE *fid;
    volume = 0e0;

    i = triangles[0];
    fid = fopen("triangles.dat", "wb");
    for (it = 0; it< 3*num_triangles; it=it+3){
        i = triangles[it];
        j = triangles[it+1];
        k = triangles[it+2];
        ri = pos[i]; rj = pos[j]; rk = pos[k];
        /* ri = Position_add(pos[i], pt, -1); */
        /* rj = Position_add(pos[j], pt, -1); */
        /* rk = Position_add(pos[k], pt, -1); */
        
        rij = Position_add(ri, rj, 1e0);
        rijk = Position_add(rij, rk, 1e0);

        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;

        rij = Position_add(rj, ri, -1e0);
        rik = Position_add(rk , ri, -1e0);

        area = cross_product(rij, rik);
        dir_norm = (rij.y*rik.z - rij.z*rik.y) + 
            (rij.z*rik.x - rij.x*rik.z) +
            (rij.x*rik.y - rij.y*rik.x);

        volume = volume + 0.5*fabs(area.x*rijk.x + area.y*rijk.y + area.z*rijk.z);
    }
    return volume/3e0;
};

void identify_obtuse(POSITION *pos, int *triangles,
       double *obtuse,  int N){
    double piby2;
    int i, j, k, it;
    double a_ij_ik, a_ji_jk, a_ki_kj;
    POSITION ri, rj, rk;
    POSITION rij, rik, rkj;
    double inp_r1, inp_r2, inp_r3;

    piby2 = 0.5*pi;
    for (it = 0; it< 3*N; it=it+3){
        i = triangles[it];
        j = triangles[it+1];
        k = triangles[it+2];

        ri = pos[i]; rj = pos[j]; rk = pos[k];

        rij = Position_add(ri, rj, -1e0);
        inp_r1 = inner_product(rij, rij);

        rik = Position_add(ri , rk, -1e0);
        inp_r2 = inner_product(rik, rik);

        rkj = Position_add(rk , rj, -1e0);
        inp_r3 = inner_product(rkj, rkj);

        a_ij_ik = (inner_product(rij,rik))/(sqrt(inp_r1*inp_r2));
        a_ij_ik = acos(a_ij_ik);

        a_ji_jk = (inner_product(rij,rkj))/(sqrt(inp_r1*inp_r3));
        a_ji_jk = acos(a_ji_jk);

        a_ki_kj = -(inner_product(rik,rkj))/(sqrt(inp_r2*inp_r3));
        a_ki_kj = acos(a_ki_kj);

       bool logic = ((a_ij_ik > piby2) ||
                (a_ji_jk > piby2) || 
                (a_ki_kj > piby2));

       /* double sum_ang = a_ij_ik + a_ji_jk + a_ki_kj; */

       /* printf( "%lf %lf %lf %lf %lf\n", a_ij_ik, */
               /* a_ji_jk, a_ki_kj, sum_ang, pi); */ 
      
        if(logic) obtuse[i/3] = 1e0;
    }
}
