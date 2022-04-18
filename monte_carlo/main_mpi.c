#include "include/global.h"
#include "include/subroutine.h"
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
/**************************************************/
int main(int argc, char *argv[]){
    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves;
    double Et[7], Ener_t;
    double vol_sph, e_t, s_t,area_sph;
    POSITION *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    MESH mesh;
    MESH mes_t;
    POSITION afm_force,spring_force[2],tot_force;
    FILE *fid;
    double *lij_t0, *obtuse;
    int *triangles;
    int *triangles_t;
    string outfolder, syscmds, log_file, outfile, para_file,filename;
    int ibydumpskip, mpi_err, mpi_rank;
    double Pole_zcoord;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume  forcex, forcey forcez area nPole_z sPole_z";
    int nPole, sPole;
    uint32_t seed_v;
    SPRING_para spring;
    //
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed_v = (uint32_t) 7*3*11*(mpi_rank+1)*rand();
    init_rng(seed_v);
    //
    outfolder = ZeroPadNumber(mpi_rank)+"/";
    cout << "I made folder "+ outfolder << endl;
    filename = outfolder + "/para_file.in";
    initialize_read_parameters(&mbrane, &afm, &mcpara, &spring, filename.c_str());
    // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
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
    afm.tip_curve = (POSITION *)calloc(afm.N, 3*sizeof(double));
    if(!mcpara.is_restart){ 
        filename = outfolder + "/input.h5";
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, filename);
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
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
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.N,mbrane.th_cr);
        filename = outfolder + "/restart.h5";
        hdf5_io_read_config((double *) Pos, (int *) mes_t.cmlist,
                (int *) mes_t.node_nbr_list, (int2 *) mes_t.bond_nbr_list,
                triangles_t, filename);
    }
    //
    double YY = mbrane.YY;
    double BB = mbrane.coef_bend;
    outfile_terminal << "# Foppl von Karman (FvK): "
        << YY*mbrane.radius*mbrane.radius/BB << endl;
    write_param(outfolder + "/para.out",mbrane,mcpara,spring);
    //
    Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] =  bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    volume_area_enclosed_membrane(Pos, triangles,
            mbrane.num_triangles,&vol_sph,&area_sph);
    double  ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
    Et[5] = spring_tot_energy_force(Pos, spring_force, mesh, spring);
    Et[6] = -mbrane.pressure*vol_sph;
    //
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4] + Et[5] + Et[6];
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
        Et[5] = spring_tot_energy_force(Pos, spring_force, mesh, spring);
        Et[6] = -mbrane.pressure*vol_sph;
        outfile_terminal << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
                " totalener = "<< mbrane.tot_energy[0] << "; volume = " << mbrane.volume[0]<< "; area = " << area_sph << endl;
        wDiag(fid, mbrane, afm, spring, mesh, i, num_moves, Et,  &afm_force,  spring_force,  area_sph,  Pos);
        if(i%mcpara.dump_skip == 0){
            outfile=outfolder+"/part_"+ ZeroPadNumber(i/mcpara.dump_skip)+".vtk";
            visit_vtk_io( (double *) Pos, triangles, mbrane.N, outfile);
            visit_vtk_io_point_data(is_attractive, mbrane.N,
                    outfile, "isattr");
            hdf5_io_dump_restart_config((double *) Pos, (int *) mesh.cmlist,
                    (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list,  
                    triangles, mbrane, outfolder);
        }
        if(i == 1*mcpara.dump_skip && !mcpara.is_restart ){
            afm.sigma = s_t;
            afm.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane, afm);
            mbrane.tot_energy[0] += e_t;
        }
        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive, 
                mbrane, mcpara, afm, spring);
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