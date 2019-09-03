#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"
#include "evolve.h"
#include "model.h"
#include <time.h>
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
  bool Exit= check_param( PARAM );
  if (Exit){return 0;}
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
  // int size_psi = NN*sizeof(double);
  pre_evolve(  ndim, TimeScheme , &TT, &dev_tt, Nblock, Nthread  ) ;
  // setup initial configuration 
  initial_configuration( PSI, PARAM ) ;
  system("exec rm -rf data");
  system("exec mkdir data"); 
  // wPSI( PSI, TT.time ) ;   
  H2D( dev_psi, PSI, ndim );
  clock_t timer=clock();
  double ttimer = timer/CLOCKS_PER_SEC;
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi, 
          &TT,  dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag, size_diag,
          BUG,  dev_bug, 
          Nblock, Nthread ) ;
  printf( "#... time evolution finished \n");
  ttimer = clock()/CLOCKS_PER_SEC-ttimer;
  printf("Total time taken %fs\n", ttimer);
  D2H( PSI, dev_psi, ndim );
  // cudaMemcpy( PSI, dev_psi, size_psi , cudaMemcpyDeviceToHost);
  wPSI( PSI, TT.time ) ;
  free_chain( &PSI, &dev_psi );
}