#ifndef subroutine_h
#define subroutine_h
//2DMd_sphere_rod.c

int main();

//forces_lj.c
void make_nlist();
void make_nlist_pf();
bool len_check_ss();
bool len_check_pf();
double cal_length();
double pairlj_ipart_energy();
double pairlj_total_energy();
double pairlj_ipart_energy_pf();
double pairlj_total_energy_pf();
/* double pairlj(); */

//forces_surf.c

double bending_energy_total();
double bending_energy_ipart();
double bending_energy_ipart_neighbour();
double stretch_energy_total();
double stretch_energy_ipart();
void identify_obtuse();

//initialise.c
void initialize_system();
void initialize_read_config();
void initialize_eval_lij_t0();
int randint();


//visit_io.c
void visit_vtk_io(); 
void visit_vtk_io_points_data();
void visit_vtk_io_cell_data();


//hdf5_io
void hdf5_io_read_config();


#endif
