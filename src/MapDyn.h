#include "model.h"
#ifndef FILE_MAPDYN_SEEN
#define FILE_MAPDYN_SEEN
// #include "main.cu"
#define MaxIter 64
/*--------------------------------------------------*/
struct MV{
  double time;    // Current time of the simulation
  double dt;      // Current dt for the adaptive scheme
  int period;      // Number of iteration of the map.
  int iorbit;
  bool istab;		// Do you want to calculate the stability of the orbit?
  // bool irel_orb;
  int mapdim;
  string guess_space;
  // int ndiag;
  // int ldiag;
};
const int size_MV = 2*sizeof( double ) + 3*sizeof( int ) + sizeof(bool);
// add the size of guess_space.
/*--------------------------------------------------*/
void map_one_iter(double *y);
bool periodic_orbit(double ytrans_all[], double fytrans[], 
                    double yall[], double fy[], double time[]);
bool IsOrbit(double y[]);
// void Jacobian(double DerM[][ndim],double x[]);
void assign_map_param(void);
void GG(double y[]);
void map_multiple_iter(double y[]);
void coordinate_transform(double *y_trans, double *y) __attribute__((weak));
void inv_coordinate_transform(double *y, double *y_trans) __attribute__((weak));
/*--------------------------------------------------*/
extern MV MM;					// Saying that MM exists globally
/*--------------------------------------------------*/
#endif /* !MAPDYN_SEEN */