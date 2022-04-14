#include "../include/global.h"
#include "../include/subroutine.h"
#include <random>
//
std::mt19937 rng;
void init_rng(){
    uint32_t seed_val;
    rng.seed(seed_val);
}
//
bool Metropolis(double DE, MCpara mcpara){
    bool yes;
    double rand;
    std::uniform_real_distribution<> rand_real(0, 1);
    yes = (DE <= 0.e0);
    if (!yes){
        rand = rand_real(rng);
        yes = rand < exp(-DE/mcpara.kBT);
    }
    return yes;
}
//
double rand_inc_theta(double th0, 
        double dfac){
    double dth;
    double tmp_th0;
    std::uniform_real_distribution<> rand_real(-1, 1);
/*     th1 = th0+dth */
    tmp_th0 = 10;
    while (tmp_th0 > pi || tmp_th0 < 0){
        dth = (pi/dfac)*(rand_real(rng)); 
        tmp_th0 = th0 + dth;
    }
    /*if th1 > np.pi or th1 < 0 : */
    /*    //#print(th0,dth,th1) */
    /*    th1 = get_rand_th1(th0) */
    return dth;
}
//
double energy_mc_3d(POSITION *pos, MESH mesh, 
        double *lij_t0, bool *is_attractive, int idx, bool *is_be_pos,
        MBRANE_para mbrane, 
        MCpara mcpara, AFM_para afm, SPRING_para spring){
    double E_b, E_s;
    double E_stick, E_spr;
    double E_vol, E_afm;
    int cm_idx, num_nbr;

    E_b = 0; E_s = 0; E_stick = 0; E_afm = 0;
    num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
    cm_idx = mesh.cmlist[idx];
    E_b = bending_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
            (int2 *) (mesh.bond_nbr_list + cm_idx), num_nbr, 
            idx, mbrane);
    
    E_b += bending_energy_ipart_neighbour(pos, 
            mesh, idx, mbrane);
    
    E_s = stretch_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
            (double *) (lij_t0 + cm_idx), num_nbr, 
            idx, mbrane);
    
    E_stick = lj_bottom_surface(pos[idx].z, is_attractive[idx], 
            mbrane);

    E_afm = lj_afm(pos[idx], afm);
    *is_be_pos = E_b > 0;

    E_spr = spring_energy(pos[idx], idx, mesh, spring);
    return E_b + E_s + E_stick + E_afm + E_spr;
}
//
int monte_carlo_3d(POSITION *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane,
                MCpara mcpara, AFM_para afm, SPRING_para spring){
    int i, j, move;
    int num_nbr, cm_idx;
    double x_o, y_o, z_o, x_n, y_n, z_n;
    double de, Et, Eini, Efin;
    double dxinc, dyinc, dzinc;
    double vol_i, vol_f;
    double de_vol,dvol,de_pressure;
    bool is_be_pos;
    std::uniform_int_distribution<uint32_t> rand_int(0,mbrane.N-1);
    std::uniform_real_distribution<> rand_real(-1, 1);
    // 
    move = 0;
    for(i = 0; i< mcpara.one_mc_iter; i++){
        int idx = rand_int(rng);
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];

        Eini = energy_mc_3d(pos, mesh,
                lij_t0, is_attractive, 
                idx, &is_be_pos, mbrane, mcpara, afm, spring);

        vol_i = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        x_o =   pos[idx].x;
        y_o =   pos[idx].y;
        z_o =   pos[idx].z;

        dxinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));
        dyinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));
        dzinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));

        x_n = x_o + dxinc;
        y_n = y_o + dyinc;
        z_n = z_o + dzinc;

        pos[idx].x = x_n;
        pos[idx].y = y_n;
        pos[idx].z = z_n;

        Efin = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx, &is_be_pos,  mbrane, mcpara, afm, spring);

        vol_f = volume_ipart(pos,
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        dvol=0.5*(vol_f-vol_i);
        de_vol = vol_energy_change(mbrane,dvol);
        de_pressure = PV_change(mbrane,dvol);
        /* printf("%d %g %g %d\n", idx, Efin, Eini, is_attractive[idx]); */
        de = (Efin - Eini) + de_vol + de_pressure;
        // if(Metropolis(de , mcpara) && is_be_pos){
        if (Metropolis(de,mcpara)){
            move = move + 1;
            mbrane.tot_energy[0] +=  de;
            mbrane.volume[0] += dvol;
        }
        else{
            pos[idx].x = x_o;
            pos[idx].y = y_o;
            pos[idx].z = z_o;
        }
    }
    return move;
}