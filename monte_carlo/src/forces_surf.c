
#include "../include/global.h"
#include "../include/subroutine.h"

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

double cotangent(POSITION si, POSITION sj, POSITION sk){

    POSITION drik, drjk, cross;
    double cot_theta, inp_1;  
    double inner_prod;

    drik = Position_add(si , sk, -1e0);
    drjk = Position_add(sj , sk, -1e0);
    cross = cross_product(drik, drjk);
    inner_prod = inner_product(drik, drjk);
    cot_theta = inner_prod/sqrt(inner_product(cross,cross));

    return cot_theta;
}

double stretch_energy_ipart(POSITION *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para){

    double HH;
    double idx_ener;
    POSITION rij;
    double mod_rij;
    int i,j;

    idx_ener = 0e0;
    HH = para.coef_str/(para.av_bond_len*para.av_bond_len);
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        rij = Position_add(pos[idx], pos[j], -1e0);
        mod_rij = sqrt(inner_product(rij, rij));

        idx_ener = idx_ener + (mod_rij - lij_t0[i])*(mod_rij - lij_t0[i]);
    }
    
    return 0.5*idx_ener*HH; 
 }

double bending_energy_ipart(POSITION *pos, 
        int *node_nbr, int2 *bond_nbr,
        int num_nbr, int idx, MBRANE_para para){

    double idx_ener;
    double BB, rij_dot_rij;
    int i, k, kp;
    int j;
    double cot_sum, sigma_i, bend_ener, curv_t0;
    POSITION cot_theta_rij, rij, cot_times_rij;

    BB = para.coef_bend;
    sigma_i = 0e0;
    curv_t0 = 2e0/para.radius;
    cot_times_rij.x = 0e0;
    cot_times_rij.y = 0e0;
    cot_times_rij.z = 0e0;

    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k  = bond_nbr[i].i1; 
        kp = bond_nbr[i].i2; 

        cot_sum=0.5*(cotangent(pos[idx],pos[j],pos[k]) +
                cotangent(pos[idx],pos[j],pos[kp]));
        rij = Position_add(pos[idx], pos[j], -1e0);

        rij_dot_rij = inner_product(rij, rij);

        cot_times_rij  = Position_add(cot_times_rij, rij, cot_sum); 

        sigma_i = sigma_i + rij_dot_rij*cot_sum;

    } 
    sigma_i *= 0.25;
    double curvature = (1e0/sigma_i)*sqrt(inner_product(cot_times_rij, cot_times_rij));

    bend_ener = 0.5*BB*sigma_i*(curvature - curv_t0);
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
    return se;
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
        inv_sqdz = sigma/(ds*ds);
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

double volume_enclosed_membrane(POSITION *pos, 
        int *triangles, int num_triangles){
    int i, j, k, it;
    double volume;
    POSITION area, rk, ri, rj;
    POSITION rij, rijk, rik;
    double inp_r1, inp_r2, inp_r3;
    volume = 0e0;
    for (it = 0; it< 3*num_triangles; it=it+3){
        i = triangles[it];
        j = triangles[it+1];
        k = triangles[it+2];
        ri = pos[i]; rj = pos[j]; rk = pos[k];

        /*ri.x = ri.x + 20; */
        /*ri.y = ri.y + 20; */
        /*ri.z = ri.z + 20; */

        /*rj.x = rj.x + 20; */
        /*rj.y = rj.y + 20; */
        /*rj.z = rj.z + 20; */

        /*rk.x = rk.x + 20; */
        /*rk.y = rk.y + 20; */
        /*rk.z = rk.z + 20; */

        rij = Position_add(ri, rj, 1e0);
        rijk = Position_add(rij, rk, 1e0);

        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;

        rij = Position_add(ri, rj, -1e0);
        rik = Position_add(rk , ri, -1e0);
        area = cross_product(rij, rik);

        volume = volume + 0.5*abs(area.x*rijk.x + area.y*rijk.y + area.z*rijk.z);

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
