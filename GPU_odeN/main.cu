#include <stdio.h>
#include <stdlib.h>
#include "odeN.h"
#include "evolve.h"
#include "model.h"
int main( void ){
  double *dev_psi, *PSI;
  alloc_odeN( &PSI, &dev_psi );
  iniconf( PSI, dev_psi );
  pre_evolve( ndim, "rnkt4" );
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi );
  printf( "#... time evolution finished \n");
  free_odeN( &PSI, &dev_psi );
}
