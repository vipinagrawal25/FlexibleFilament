#include "model.h"
#ifndef FILE_MAPDYN_SEEN
#define FILE_MAPDYN_SEEN
// #include "main.cu"
#define MaxIter 64
/*--------------------------------------------------*/
void periodic_orbit(double y0[], double fyvec[]);
bool IsOrbit(double y[]);
void Jacobian(double DerM[][ndim],double x[]);
void assign_map_param(void);
void GG(double y[]);
void map_multiple_iter(double y[]);
/*--------------------------------------------------*/
struct MV{
  double time;    // Current time of the simulation
  double dt;      // Current dt for the adaptive scheme
  int period;      // Number of iteration of the map.
  int iorbit;
  bool istab;		// Do you want to calculate the stability of the orbit?
  bool irel_orb;
  // int ndiag;
  // int ldiag;
};
const int size_MV = 2*sizeof( double ) + 2*sizeof( int ) + 2*sizeof(bool);
/*--------------------------------------------------*/
extern MV MM;					// Saying that MM exists globally
/*--------------------------------------------------*/
#endif /* !MAPDYN_SEEN */