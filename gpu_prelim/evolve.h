#ifndef FILE_EVOLVE_SEEN
#define FILE_EVOLVE_SEEN
#include "model.h"
#include "cuda.h"
/*--------------------------------------------------*/
struct EV{
  double tmax;
  double time;
  double dt;
  double tdiag;
  int ndiag;
};
const int size_EV = 4*sizeof( double ) + sizeof( int );
/* These are device pointers I need for the time-stepping code. Their size
depends on what kind of algorithm I use. For example, for euler
time stepping dev_kk is allocated to the same size as PSI. For 4th order
Runge-Kutta it is allocated to 5 times the size of PSI. */
void pre_evolve( int Nsize, char *algo, EV *TT,  EV **dev_tt );
void wevolve( EV *TT, char *fname ) ;
void evolve( double PSI[], double dev_psi[], 
             EV TT, EV *dev_tt,
             MPARAM *dev_param ,
             double DIAG[], double dev_diag[], int size_diag,
             CRASH BUG,  CRASH *dev_bug ) ;
 #endif /* !EVOLVE_SEEN */
