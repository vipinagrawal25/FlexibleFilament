#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"
#include "model.h"
#include <time.h>
#include "evolve.h"
#include<iostream>
using namespace std;
/*==============================================*/
int main( void ){
  CRASH BUG, *dev_bug ; 
  double *dev_psi, *PSI, *dev_vel, *VEL;
  MPARAM PARAM, *dev_param;
  double *dev_diag, *DIAG;
  EV TT, *dev_tt;
  int Nblock, Nthread;
  /*------------------------------------------------*/
  double MaxThread=128;
  if ( NN <= MaxThread ) {
    Nthread= 1 ;
    Nblock = NN ;
    } else{
    // Otherwise: we launch the threads differently:
    Nthread = MaxThread;
    Nblock = (NN+MaxThread-1)/MaxThread;
  }
  printf( "#-I shall launch %d blocks, each with %d threads\n", Nblock, Nthread );
  /*------------------------------------------------*/
  set_crash(  &BUG, &dev_bug ); 
  alloc_chain( &PSI, &dev_psi );
  alloc_chain( &VEL, &dev_vel );
  set_param( &PARAM, &dev_param ) ;
  bool Exit= check_param( PARAM );
  if (Exit){return 0;}
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
  // int size_psi = NN*sizeof(double);
  pre_evolve(  ndim, TimeScheme , &TT, &dev_tt, Nblock, Nthread  ) ;
  // setup initial configuration 
  initial_configuration( PSI, PARAM );
  // system("exec rm -rf data");    
  system("exec mkdir data"); 
  // wPSI( PSI, TT.time ) ;   
  H2D( dev_psi, PSI, ndim );
  H2D(dev_vel, VEL, ndim);
  clock_t timer=clock();
  double ttimer = timer/CLOCKS_PER_SEC;
  printf( " #starting time evolution ...\n ");
  evolve( PSI, dev_psi, 
          VEL, dev_vel,
          &TT,  dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag, size_diag,
          BUG,  dev_bug, 
          Nblock, Nthread ) ;
  printf( "#... time evolution finished \n");
  post_evolve( TimeScheme  );
  ttimer = clock()/CLOCKS_PER_SEC-ttimer;
  cout << "Total time taken: " << ttimer << "s"<< endl;
  D2H( PSI, dev_psi, ndim );
  D2H(VEL,dev_vel,ndim );
  // cudaMemcpy( PSI, dev_psi, size_psi , cudaMemcpyDeviceToHost);
  wPSI( PSI, VEL, TT.time ) ;
  free_chain( &PSI, &dev_psi );
}