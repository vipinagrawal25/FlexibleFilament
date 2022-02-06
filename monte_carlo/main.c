#include "include/global.h"
#include "include/subroutine.h"
bool Metropolis(double DE, MCpara mcpara){
    bool yes;
    double rand;

    yes = (DE <= 0.e0);
    if (!yes){
        rand = drand48();
        yes = rand < exp(-DE/mcpara.kBT);
    }
    return yes;
}

double rand_inc_theta(double th0, 
        double dfac){
    double dth;
    double tmp_th0;


/*     th1 = th0+dth */
    tmp_th0 = 10;
    while (tmp_th0 > pi || tmp_th0 < 0){
        dth = (pi/dfac)*(2*drand48() - 1); 
        tmp_th0 = th0 + dth;
    }
    /*if th1 > np.pi or th1 < 0 : */
    /*    //#print(th0,dth,th1) */
    /*    th1 = get_rand_th1(th0) */
    return dth;
}

double energy_mc_3d(POSITION *pos, MESH mesh, 
        double *lij_t0, bool *is_attractive, int idx, 
        MBRANE_para mbrane, 
        MCpara mcpara, AFM_para afm){
    double E_b, E_s;
    double E_stick;
    double E_vol, E_afm;
    int cm_idx, num_nbr;


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

    E_afm = rep_lj_afm(pos[idx], afm);

    if(E_b < 0e0){
        return 2e14;
    } else {
        return E_b + E_s + E_stick + E_afm;
    }
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

    ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    KAPPA = mbrane.coef_vol_expansion;
    move = 0;
    for(i = 0; i< mcpara.mc_iter; i++){
        int idx = randint(mbrane.N);
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];

        Eini = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx,  mbrane, mcpara, afm);

        vol_i = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        x_o =   pos[idx].x;
        y_o =   pos[idx].y;
        z_o =   pos[idx].z;

        dxinc = (mcpara.delta/mcpara.dfac)*(2*drand48() - 1);
        dyinc = (mcpara.delta/mcpara.dfac)*(2*drand48() - 1);
        dzinc = (mcpara.delta/mcpara.dfac)*(2*drand48() - 1);

        x_n = x_o + dxinc;
        y_n = y_o + dyinc;
        z_n = z_o + dzinc;

        pos[idx].x = x_n;
        pos[idx].y = y_n;
        pos[idx].z = z_n;

        Efin = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx,   mbrane, mcpara, afm);

        vol_f = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        dvol =  0.5*(vol_f - vol_i);
        de_vol = (2*dvol/(ini_vol*ini_vol))*(mbrane.volume[0]  - ini_vol);
        de_vol = KAPPA*de_vol;


        de = (Efin - Eini) + de_vol;
        if(Metropolis(de , mcpara)){
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

int monte_carlo_surf2d(POSITION *Pos, 
        Neighbours *neib, LJpara para, 
        MCpara mcpara){
    int i, j, move;
    double x_o, y_o, x_n, y_n;
    double de, Et, Eini, Efin;
    double dxinc, dyinc;
    bool is_sph, is_cart;

    is_sph = false;
    is_cart = false;
    if(strcmp(mcpara.metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(mcpara.metric, "cart") == 0){
        is_cart = true;
    }


    move = 0;

    for(i = 0; i< mcpara.mc_iter; i++){
    /* while(move < mcpara.mc_iter){ */
        int idx = randint(para.N);
        Eini =  pairlj_ipart_energy(Pos, neib[idx].list_ss,
                neib[idx].cnt_ss, idx, para, mcpara.metric);
        /* Eini =  pairlj_ipart_energy_pf(Pos, idx, para, */ 
                /* mcpara.metric); */
        x_o =   Pos[idx].x;
        y_o =   Pos[idx].y;

        if(is_cart){
            dxinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
            dyinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
            x_n = fmod((x_o + dxinc + 30*para.len), para.len);
            y_n = fmod((y_o + dyinc + 30*para.len), para.len);
            Pos[idx].x = x_n;
            Pos[idx].y = y_n;
        }
        if(is_sph){
            dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
            dyinc = (para.sigma/mcpara.dfac)*(2*drand48() - 1);
            x_n = x_o + dxinc;
            y_n = fmod((y_o + dyinc + 30*2*pi), 2*pi);
            Pos[idx].x = x_n;
            Pos[idx].y = y_n;
        }
        /* Efin =  pairlj_ipart_energy_pf(Pos, idx, para, */ 
                /* mcpara.metric); */
 
        Efin =  pairlj_ipart_energy(Pos, neib[idx].list_ss,
                neib[idx].cnt_ss, idx, para, mcpara.metric);
        de = Efin - Eini;
        if(Metropolis(de, mcpara)){
            move = move + 1;
        }
        else{
            Pos[idx].x = x_o;
            Pos[idx].y = y_o;
        } 
    }
    return move;
}

int dump_config(POSITION *Pos, double len, 
        int iter, int N){
    char part_file[80];
    int i;
    double x_n, y_n;
    FILE *fid;


    sprintf(part_file,"output/part_%04d.dat",iter);

    fid = fopen(part_file, "wb");

    for(i=0; i < N; i++){
        x_n = fmod((Pos[i].x + 30*len), len);
        y_n = fmod((Pos[i].y + 30*len), len);
        fprintf(fid, "%lf  %lf \n", x_n, y_n);
        fflush(fid);
    }
    fclose(fid);
    return 1;
}

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
//     mcpara.mc_iter = 10*para.N;
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

int main(int argc, char *argv[]){
    int i, iterations, num_moves;
    double Ener_s, Ener_b, Ener_t;
    double vol_sph, e_t, s_t;
    POSITION *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    MESH mesh;
    FILE *fid;
    double *lij_t0, *obtuse;
    int *triangles;
    char syscmds[128], outfile[89], outfolder[89], para_file[32];
    char log_file[64];

    if(argc!=3){
        printf("\n\n mayday.. requires an argument <parameter file> <output folder>\n\n");
        exit(0);
    }else{
        sscanf(argv[1],"%s",&para_file);
        sscanf(argv[2],"%s",&outfolder);
    }

    sprintf(syscmds, "mkdir %s",outfolder);
    system(syscmds);
    sprintf(syscmds, "cp %s %s/%s",para_file,outfolder,para_file);
    system(syscmds);

    // read the input file
    initialize_read_parameters(&mbrane, &afm, &mcpara, para_file);


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
    lij_t0 = (double *)calloc(mbrane.num_nbr, sizeof(double));
    triangles = (int *)calloc(mbrane.num_nbr, sizeof(int)); 
    obtuse = (double *)calloc(mbrane.num_triangles, sizeof(double)); 
    is_attractive = (bool *)calloc(mbrane.N, sizeof(bool)); 
    afm.tip_curve = (POSITION *)calloc(afm.N, 3*sizeof(double)); 

    s_t = afm.sigma; 
    afm.sigma = 0.00;
    e_t = afm.epsilon; 
    afm.epsilon = 0.0;

    hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
            (int *) mesh.node_nbr_list, (int *) mesh.bond_nbr_list, 
            triangles, "input/input.h5" );

    initialize_afm_tip(afm);
    sprintf(log_file, "%s/afm_tip.vtk", outfolder);
    visit_vtk_io_afm_tip((double *) afm.tip_curve, 
            afm.N, log_file);
    identify_attractive_part(Pos, is_attractive, mbrane.N);
    initialize_eval_lij_t0(Pos, mesh, lij_t0, mbrane);

    Ener_s =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Ener_b =  bending_energy_total(Pos, mesh, mbrane);
    Ener_t = Ener_s + Ener_b;
    mbrane.tot_energy[0] = Ener_t;

    vol_sph  = volume_enclosed_membrane(Pos, triangles, 
            mbrane.num_triangles);

    mbrane.volume[0] = vol_sph;

    iterations = 60000;

    sprintf(log_file, "%s/mc_log", outfolder);
    fid = fopen(log_file, "a");
    num_moves = 0;
    for(i=0; i<iterations; i++){
        if(i%100 == 0){
            fprintf(stderr, "iter, AcceptedMoves, ener, ener, volume: %d %d %g %g\n",
                    i, num_moves, mbrane.tot_energy[0], mbrane.volume[0]);
            identify_obtuse(Pos, triangles, obtuse, mbrane.num_triangles);
            sprintf(outfile,"%s/part_%05d.vtk",outfolder,i);
            visit_vtk_io( (double *) Pos, triangles, 
                    mbrane.N, outfile, "miketesting");
            visit_vtk_io_cell_data(obtuse, mbrane.num_triangles,
                    outfile, "is90");
            visit_vtk_io_point_data(is_attractive, mbrane.N,
                    outfile, "attr");
            fprintf(fid, " %d %d %g %g\n",
                    i, num_moves, mbrane.tot_energy[0], mbrane.volume[0]);
            fflush(fid);
        }
        if(i > 6000){
            afm.sigma = s_t;
            afm.epsilon = e_t;
        }
        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive, 
                mbrane, mcpara, afm);
    }
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    free(mesh.bond_nbr_list);
    return 0;
}
