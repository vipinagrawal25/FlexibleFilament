#ifndef subroutine_h
#define subroutine_h
#include <iostream>
#include <cstdlib>
//forces_lj.c
void make_nlist();
void make_nlist_pf();
// bool len_check_ss();
bool len_check_ss(POSITION s1, POSITION s2, double len, double new_rc);
bool len_check_pf(POSITION s1, POSITION s2, double len, double new_rc, char* metric);
double cal_length();
double pairlj_ipart_energy(POSITION *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric);
double pairlj_total_energy();
double pairlj_ipart_energy_pf(POSITION *Pos, 
        int i_p, LJpara para, char *metric);
double pairlj_total_energy_pf(POSITION *Pos, LJpara para, char *metric);
/* double pairlj(); */

//forces_surf.c
double bending_energy_total(POSITION *pos, MESH mesh, MBRANE_para para);
double bending_energy_ipart(POSITION *pos, 
        int *node_nbr, int2 *bond_nbr,
        int num_nbr, int idx, MBRANE_para para);
double bending_energy_ipart_neighbour(POSITION *pos, 
        MESH mesh, int idx, MBRANE_para para);
double stretch_energy_total(POSITION *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para);
void identify_obtuse(POSITION *pos, int *triangles, 
       double *obtuse,  int N);
double stretch_energy_ipart(POSITION *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para);

//initialise.c
void initialize_system();
void initialize_read_config();
void initialize_eval_lij_t0();
int randint(int n);
//visit_io.c
void visit_vtk_io(double *points, 
        int *triangles, 
        int Np, char filename[], 
        char dataname[]);
void visit_vtk_io_points_data();
void visit_vtk_io_cell_data(double *data, 
        int Np, char filename[], 
        char dataname[]);
//hdf5_io
void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int *bond_nbr, int *triangles,
        char input_file[]);
//misc
bool FileExists(const std::string &s);
#endif