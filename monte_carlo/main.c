#include "include/global.h"
#include "include/subroutine.h"
#include <unistd.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <random>
/**************************************************/
template<typename T>
inline string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}
//
void wHeader(FILE *fid, MBRANE_para mbrane, AFM_para afm, SPRING_para spring){
    string log_headers = "#iter acceptedmoves total_ener stretch_ener bend_ener stick_ener ";
    if(afm.icompute!=0){log_headers+="afm_ener ";}
    if (spring.icompute!=0){log_headers+="spring_energy ";}
    if(fabs(mbrane.coef_vol_expansion)>1e-16){log_headers+="ener_volume ";}
    if (fabs(mbrane.pressure)>1e-16){log_headers+="pressure_ener ";}
    if (afm.icompute!=0){log_headers+="afm_fx, afm_fy afm_fz ";}
    if (spring.icompute!=0){log_headers+="spr_north.z spr_south.z ";}
    log_headers+=" area nPole_z sPole_z";
    fprintf(fid, "%s\n", log_headers.c_str());
    fflush(fid);
}
//
void wDiag(FILE *fid, MBRANE_para mbrane, AFM_para afm, SPRING_para spring, MESH mesh, 
            int i, int num_moves, double *Et,
            POSITION *afm_force, POSITION *spring_force, double area_sph, POSITION *Pos){
    fprintf(fid, " %d %d %g %g %g %g", i, num_moves, mbrane.tot_energy[0], Et[0], Et[1], Et[2]);
    if(afm.icompute!=0){fprintf(fid, " %g", Et[3]);}
    if (spring.icompute!=0){fprintf(fid, " %g", Et[5]);}
    if(fabs(mbrane.coef_vol_expansion)>1e-16){fprintf(fid, " %g", Et[4]);}
    if (fabs(mbrane.pressure)>1e-16){fprintf(fid, " %g", Et[6]);}
    if (afm.icompute!=0){fprintf(fid, " %g %g %g", afm_force->x,afm_force->y,afm_force->z );}
    if (spring.icompute!=0){fprintf(fid, " %g %g", spring_force[0].z,spring_force[1].z);}
    fprintf(fid, " %g %g %g\n",area_sph,Pos[mesh.nPole].z,Pos[mesh.sPole].z);
    fflush(fid);
}
//
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
    SPRING_para spring;
    MESH mesh;
    MESH mes_t;
    POSITION afm_force,spring_force[2],tot_force;
    FILE *fid;
    double *lij_t0, *obtuse;
    int *triangles;
    int *triangles_t;
    string outfolder, syscmds, log_file, outfile, para_file,filename;
    // int nPole,sPole;
    int ibydumpskip;
    double Pole_zcoord;
    if(argc!=3){
        printf("\n\n mayday.. requires an argument <parameter file> <output folder>\n\n");
        exit(0);
    }else{
        para_file=argv[1];
        outfolder=argv[2];
    }
    // read the input file
    initialize_read_parameters(&mbrane, &afm, &mcpara, &spring, para_file);
    //
    syscmds="mkdir "+outfolder;
    system(syscmds.c_str());
    filename = outfolder + "/para.out";
    write_param(filename,mbrane,mcpara);
    syscmds="cp "+para_file+" "+outfolder+"/";
    system(syscmds.c_str());
    init_rng();
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
    //
    if(!mcpara.is_restart){
        s_t = afm.sigma; 
        afm.sigma = 0.00;
        e_t = afm.epsilon;
        afm.epsilon = 0.0;
    }
    if(!mcpara.is_restart){
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, "input/input.h5");
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.N, mbrane.th_cr);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
    }else{
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, "input/input.h5");
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.N, mbrane.th_cr);
        hdf5_io_read_config((double *) Pos, (int *) mes_t.cmlist,
                (int *) mes_t.node_nbr_list, (int2 *) mes_t.bond_nbr_list,
                triangles_t, outfolder+"/restart.h5");
    }
    //
    double YY=mbrane.YY; 
    double BB = mbrane.coef_bend;
    cout << "# Foppl von Karman (FvK): " << YY*mbrane.radius*mbrane.radius/BB << endl;
    Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] =  bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    //
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
    wHeader(fid,mbrane,afm,spring);
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
        cout << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
                " totalener = "<< mbrane.tot_energy[0] << "; volume = " << mbrane.volume[0]<< "; area = " << area_sph << endl;
        wDiag(fid, mbrane, afm, spring, mesh, i, num_moves, Et,  &afm_force,  spring_force,  area_sph,  Pos);
        if(i%mcpara.dump_skip == 0){
            outfile=outfolder+"/part_"+ ZeroPadNumber(i/mcpara.dump_skip)+".vtk";
            visit_vtk_io( (double *) Pos, triangles,
                    mbrane.N, outfile);
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
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    free(mesh.bond_nbr_list);
    return 0;
}