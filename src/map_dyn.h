#ifndef FILE_MAP_DYN_SEEN
#define FILE_MAP_DYN_SEEN
#include "model.h"
// #include "main.cu"
/*--------------------------------------------------*/
struct MV{
  double time;    // Current time of the simulation
  double dt;      // Current dt for the adaptive scheme
  double tdiag;   // Next time to save the diagnostic
  int niter;      // Number of iteration of the map.
  int ndiag;
  int ldiag;
};
const int size_EV = 3*sizeof( double ) + 3*sizeof( int );
/* These are device pointers I need for the time-stepping code. Their size
depends on what kind of algorithm I use. For example, for euler
time stepping dev_kk is allocated to the same size as PSI. For 4th order
Runge-Kutta it is allocated to 5 times the size of PSI. */
void ode2map(int Nsize, char *algo, MV *MM)
#endif /* !EVOLVE_SEEN */