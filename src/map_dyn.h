#ifndef FILE_MAP_DYN_SEEN
#define FILE_MAP_DYN_SEEN
// #include "main.cu"
struct MV{
  double time;    // Current time of the simulation
  double dt;      // Current dt for the adaptive scheme
  int period;      // Number of iteration of the map.
  // int ndiag;
  // int ldiag;
};
const int size_MV = 2*sizeof( double ) + sizeof( int );
/*--------------------------------------------------*/
/* These are device pointers I need for the time-stepping code. Their size
depends on what kind of algorithm I use. For example, for euler
time stepping dev_kk is allocated to the same size as PSI. For 4th order
Runge-Kutta it is allocated to 5 times the size of PSI. */
#endif /* !EVOLVE_SEEN */