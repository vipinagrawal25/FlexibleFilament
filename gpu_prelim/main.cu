#include <stdio.h>
#include <stdlib.h>
#include "chain.h"
#include "evolve.h"
#include "model.h"
int main( void ){
  double *dev_psi, *PSI;
  alloc_chain( &PSI, &dev_psi );
  iniconf( PSI, dev_psi );
  pre_evolve( ndim, "euler" );
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi );
  printf( "#... time evolution finished \n");
  free_chain( &PSI, &dev_psi );
}
