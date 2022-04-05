#include "include/global.h"
#include "include/subroutine.h"
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include<iomanip>
#include <fstream>

std::mt19937 rng;
void init_rng(uint32_t seed_val){
    rng.seed(seed_val);
}

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

double energy_mc_3d(POSITION *pos, MESH mesh, 
        double *lij_t0, bool *is_attractive, int idx, bool *is_be_pos,
        MBRANE_para mbrane, 
        MCpara mcpara, AFM_para afm){
    double E_b, E_s;
    double E_stick;
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
            mbrane.pos_bot_wall, mbrane.epsilon, mbrane.sigma);

    E_afm = lj_afm(pos[idx], afm);
    *is_be_pos = E_b > 0;

    return E_b + E_s + E_stick + E_afm;
}

int monte_carlo_3d(POSITION *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane, 
                MCpara mcpara, AFM_para afm){
    int i, j, move;
    int num_nbr, cm_idx;
    double x_o, y_o, z_o, x_n, y_n, z_n;
    double de, Et, Eini, Efin;
    double dxinc, dyinc, dzinc;
    double vol_i, vol_f;
    double dvol, de_vol, ini_vol;
    double KAPPA;
    bool is_be_pos;
    std::uniform_int_distribution<uint32_t> rand_int(0,mbrane.N-1);
    std::uniform_real_distribution<> rand_real(-1, 1);
    ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    KAPPA = mbrane.coef_vol_expansion;
    move = 0;
    for(i = 0; i< mcpara.one_mc_iter; i++){
        int idx = rand_int(rng);
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];

        Eini = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx, &is_be_pos, mbrane, mcpara, afm);

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
                idx, &is_be_pos,  mbrane, mcpara, afm);

        vol_f = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        dvol =  0.5*(vol_f - vol_i);
        de_vol = (2*dvol/(ini_vol*ini_vol))*(mbrane.volume[0]  - ini_vol)
            + (dvol/ini_vol)*(dvol/ini_vol);
        de_vol = KAPPA*de_vol;
        /* printf("%d %g %g %d\n", idx, Efin, Eini, is_attractive[idx]); */
        de = (Efin - Eini) + de_vol;
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

// int monte_carlo_surf2d(POSITION *Pos, 
//         Neighbours *neib, LJpara para, 
//         MCpara mcpara){
//     int i, j, move;
//     double x_o, y_o, x_n, y_n;
//     double de, Et, Eini, Efin;
//     double dxinc, dyinc;
//     bool is_sph, is_cart;

//     is_sph = false;
//     is_cart = false;
//     if(strcmp(mcpara.metric, "sph") == 0){
//         is_sph = true;
//     }
//     if(strcmp(mcpara.metric, "cart") == 0){
//         is_cart = true;
//     }


//     move = 0;

//     for(i = 0; i< mcpara.one_mc_iter; i++){
//     /* while(move < mcpara.one_mc_iter){ */
//         int idx = randint(para.N);
//         Eini =  pairlj_ipart_energy(Pos, neib[idx].list_ss,
//                 neib[idx].cnt_ss, idx, para, mcpara.metric);
//         /* Eini =  pairlj_ipart_energy_pf(Pos, idx, para, */ 
//                 /* mcpara.metric); */
//         x_o =   Pos[idx].x;
//         y_o =   Pos[idx].y;

//         if(is_cart){
//             dxinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
//             dyinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
//             x_n = fmod((x_o + dxinc + 30*para.len), para.len);
//             y_n = fmod((y_o + dyinc + 30*para.len), para.len);
//             Pos[idx].x = x_n;
//             Pos[idx].y = y_n;
//         }
//         if(is_sph){
//             dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
//             dyinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
//             x_n = x_o + dxinc;
//             y_n = fmod((y_o + dyinc + 30*2*pi), 2*pi);
//             Pos[idx].x = x_n;
//             Pos[idx].y = y_n;
//         }
//         /* Efin =  pairlj_ipart_energy_pf(Pos, idx, para, */ 
//                 /* mcpara.metric); */
 
//         Efin =  pairlj_ipart_energy(Pos, neib[idx].list_ss,
//                 neib[idx].cnt_ss, idx, para, mcpara.metric);
//         de = Efin - Eini;
//         if(Metropolis(de, mcpara)){
//             move = move + 1;
//         }
//         else{
//             Pos[idx].x = x_o;
//             Pos[idx].y = y_o;
//         } 
//     }
//     return move;
// }
// int thermalize(){
//     int i, iterations, num_moves;
//     double Ener;
//     char *conf;
//     POSITION *Pos;
//     Neighbours *neib;
//     LJpara  para;
//     MCpara mcpara;
//     FILE *fid;


//    /* define all the paras */ 

//     para.N  = 5120;
//     para.len = 2*pi;
//     para.epsilon = 1;
//     /* para.sigma = para.len/sqrt((double)para.N); */

//     para.sigma = sqrt(8*pi/(2*para.N-4));
//     para.r_cut = 4*para.sigma;

//     // define the monte carlo parameters
//     mcpara.dfac  = 32;
//     mcpara.one_mc_iter = 10*para.N;
//     mcpara.kBT = 1;
//     mcpara.metric = "sph";
//     conf = "regular";



//     Pos = calloc(para.N, sizeof(POSITION));
//     neib = calloc(para.N, sizeof(Neighbours));

//     initialize_system(Pos, para, 
//             mcpara);


//     printf("%lf  %lf \n", para.sigma, para.r_cut);

//     make_nlist_pf(Pos, neib,  para, mcpara.metric);

//     /* I have the neighbour list now let's check them */

//      for(i=0; i<para.N; i++){ 
//         /* printf("%d\n", neib[i].cnt_ss); */
//     /* } */

//     Ener = pairlj_total_energy_pf(Pos,  para, 
//             mcpara.metric);

//     iterations = 60000;
//     system("touch output/mc_log");
//     fid = fopen("output/mc_log", "a");
//     for(i=0; i<iterations; i++){
//         if(i%1000 == 0){
//             fprintf(stderr, " iter, AcceptedMoves, Energy: %d %d %g\n",
//                     i, num_moves, Ener);
//             dump_config(Pos,  para.len, i, para.N);
//         }
//         fprintf(fid, " %d %d %g\n",
//                 i, num_moves, Ener);
//         fflush(fid);
//         num_moves = monte_carlo_surf2d(Pos, neib, 
//                 para, mcpara);
//         make_nlist_pf(Pos, neib,  para, mcpara.metric);

//         Ener = pairlj_total_energy(Pos, neib, para, 
//                 mcpara.metric);
//         /*       /1* "a comment"; *1/ */
//     }
//     fclose(fid);
//     free(Pos);
//     free(neib);
//     return 0;
// }

template<typename T>
inline string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}

int main(int argc, char *argv[]){
    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves;
    double Et[5], Ener_t;
    double vol_sph, e_t, s_t,area_sph;
    POSITION *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    MESH mesh;
    MESH mes_t;
    POSITION afm_force;
    FILE *fid;
    double *lij_t0, *obtuse;
    int *triangles;
    int *triangles_t;
    string outfolder, log_file, outfile, filename;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume  forcex, forcey forcez area nPole_z sPole_z";
    int nPole, sPole;
    int ibydumpskip, mpi_err, mpi_rank;
    uint32_t seed_v;
    double Pole_zcoord;
    //
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed_v = 7*3*11*(mpi_rank+1)*rand();
    init_rng(seed_v);
    //
    outfolder = ZeroPadNumber(mpi_rank);
    cout << outfolder << endl;
    filename = outfolder + "/para_file.in";
    initialize_read_parameters(&mbrane, &afm, &mcpara, filename.c_str());
    // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"terminal.out", ios::app);
    //
    /* define all the paras */ 
    mbrane.volume = (double *)calloc(1, sizeof(double)); 
    mbrane.volume[0] = (4./3.)*pi*pow(mbrane.radius,3);
    mbrane.tot_energy = (double *)calloc(1, sizeof(double));
    mbrane.tot_energy[0] = 0e0;
    // allocate arrays
    Pos = (POSITION *)calloc(mbrane.N, sizeof(POSITION));
    mesh.cmlist = (int *)calloc(mbrane.N+1, sizeof(int));
    mesh.node_nbr_list = (int *)calloc(mbrane.num_nbr, sizeof(int)); 
    mesh.bond_nbr_list = (int2 *)calloc(mbrane.num_nbr, sizeof(int2));
    mes_t.cmlist = (int *)calloc(mbrane.N+1, sizeof(int));
    mes_t.node_nbr_list = (int *)calloc(mbrane.num_nbr, sizeof(int)); 
    mes_t.bond_nbr_list = (int2 *)calloc(mbrane.num_nbr, sizeof(int2));
    lij_t0 = (double *)calloc(mbrane.num_nbr, sizeof(double));
    triangles = (int *)calloc(mbrane.num_nbr, sizeof(int));
    triangles_t = (int *)calloc(mbrane.num_nbr, sizeof(int));
    obtuse = (double *)calloc(mbrane.num_triangles, sizeof(double));
    is_attractive = (bool *)calloc(mbrane.N, sizeof(bool));

    //
    if(!mcpara.is_restart){
        s_t = afm.sigma; 
        afm.sigma = 0.00;
        e_t = afm.epsilon; 
        afm.epsilon = 0.0;
    }
    if(!mcpara.is_restart){
        
        filename = outfolder + "/input.h5";
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, filename);
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane);
        identify_attractive_part(Pos, is_attractive, mbrane.N,mbrane.th_cr);
        max(&nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&sPole,&Pole_zcoord,Pos,mbrane.N);
    }else{
        filename = outfolder + "/input.h5";
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, filename);
        max(&nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&sPole,&Pole_zcoord,Pos,mbrane.N);
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane);
        identify_attractive_part(Pos, is_attractive, mbrane.N,mbrane.th_cr);
        filename = outfolder + "/restart.h5";
        hdf5_io_read_config((double *) Pos, (int *) mes_t.cmlist,
                (int *) mes_t.node_nbr_list, (int2 *) mes_t.bond_nbr_list,
                triangles_t, filename);
    }
    //
    // double HH = mbrane.coef_str/(mbrane.av_bond_len*mbrane.av_bond_len);
    double YY = mbrane.YY;
    double BB = mbrane.coef_bend;
    outfile_terminal << "# Foppl von Karman (FvK): "
         << YY*mbrane.radius*mbrane.radius/BB << endl;
    filename = outfolder + "/para.out";
    fid = fopen(filename.c_str(), "w");
    fprintf(fid, " Fopl von karman %lf", YY*mbrane.radius*mbrane.radius/BB);
    fclose(fid);
    //
    Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] =  bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    //
    volume_area_enclosed_membrane(Pos, triangles, 
            mbrane.num_triangles,&vol_sph,&area_sph);
    double  ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4];
    mbrane.tot_energy[0] = Ener_t;
    mbrane.volume[0] = vol_sph;
    //
    log_file=outfolder+"/mc_log";
    fid = fopen(log_file.c_str(), "a");
    if(!mcpara.is_restart)fprintf(fid, "%s\n", log_headers);
    num_moves = 0;
    for(i=0; i < mcpara.tot_mc_iter; i++){
        Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
        Et[1] =  bending_energy_total(Pos, mesh, mbrane);
        Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
        volume_area_enclosed_membrane(Pos, triangles, mbrane.num_triangles, &vol_sph, &area_sph);
        Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
        // Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4];
        // mbrane.tot_energy[0] = Ener_t;
        if(mpi_rank == 0){
        outfile_terminal << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
                " totalener = "<< mbrane.tot_energy[0] << "; volume = " << mbrane.volume[0]<< "; area = " << area_sph << endl;
        }
        fprintf(fid, " %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
                    i, num_moves, mbrane.tot_energy[0], Et[0], Et[1], Et[2], Et[3], Et[4],
                    afm_force.x, afm_force.y, afm_force.z,area_sph,Pos[nPole].z,Pos[sPole].z);
        fflush(fid);
        if(i%mcpara.dump_skip == 0){
            outfile=outfolder+"/part_"+ ZeroPadNumber(i/mcpara.dump_skip)+".vtk";
            // sprintf(outfile,"%s/part_%05d.vtk",outfolder,(int)i/mcpara.dump_skip);
            // identify_obtuse(Pos, triangles, obtuse, mbrane.num_triangles);
            visit_vtk_io( (double *) Pos, triangles, mbrane.N, outfile);
            hdf5_io_dump_restart_config((double *) Pos, (int *) mesh.cmlist,
                    (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list,  
                    triangles, mbrane, outfolder);
        }
        if(i == 10*mcpara.dump_skip && !mcpara.is_restart ){
            afm.sigma = s_t;
            afm.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane, afm);
            mbrane.tot_energy[0] += e_t;
        }
        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive, 
                mbrane, mcpara, afm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    free(mesh.bond_nbr_list);
    mpi_err = MPI_Finalize();

    outfile_terminal.close();
    return 0;
}
