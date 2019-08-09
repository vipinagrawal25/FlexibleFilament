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
  int Nblock, Nthread;
  /*------------------------------------------------*/
  if ( NN < 128 ) {
    Nthread= 1;
    Nblock = NN ;
    } else{
    // Otherwise: we launch the threads differently:
    Nthread = 128;
    Nblock = (NN+127)/128;
  }
  printf( "#-I shall launch %d blocks, each with %d threads\n", Nblock, Nthread );
  /*------------------------------------------------*/

  set_crash(  &BUG, &dev_bug ); 
  alloc_chain( &PSI, &dev_psi );
  set_param( &PARAM, &dev_param ) ;
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
   // pre_evolve(  ndim, "euler" , &TT, &dev_tt, Nblock, Nthread  ) ;
  // pre_evolve(  ndim, "rnkt4" , &TT, &dev_tt, Nblock, Nthread  ) ;
  pre_evolve(  ndim, "rnkf45" , &TT, &dev_tt, Nblock, Nthread  ) ;
  // setup initial configuration 
  initial_configuration( PSI, PARAM ) ;
  wPSI( PSI, TT.time ) ; 
  H2D( dev_psi, PSI, ndim );
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi, 
          &TT,  dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag, size_diag,
          BUG,  dev_bug, 
          Nblock, Nthread ) ;
  printf( "#... time evolution finished \n");
  H2D( dev_psi, PSI, ndim );
  wPSI( PSI, TT.time ) ; 
  free_chain( &PSI, &dev_psi );
}
