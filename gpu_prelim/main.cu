#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"
#include "evolve.h"
#include "model.h"
/*==============================================*/
int main( void ){
  CRASH BUG, *dev_bug ; 
  double *dev_psi, *PSI;
  MPARAM PARAM, *dev_param;
  double *dev_diag, *DIAG;
  EV TT, *dev_tt;
  /*------------------------------------------------*/
  set_crash(  &BUG, &dev_bug ); 
  alloc_chain( &PSI, &dev_psi );
  set_param( &PARAM, &dev_param ) ;
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
  pre_evolve(  ndim, "euler" , &TT, &dev_tt  ) ;
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi, 
             TT,  dev_tt,
             dev_param ,
             DIAG, dev_diag, size_diag,
             BUG,  dev_bug ) ;
  printf( "#... time evolution finished \n");
  free_chain( &PSI, &dev_psi );
}
