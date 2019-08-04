#include <stdio.h>
#include <stdlib.h>
//#include "chain.h"
//#include "evolve.h"
#include "model.h"
#include "cuda.h"
int main( void ){
  // CRASH BUG, *dev_bug ; 
  //double *dev_psi, *PSI;
  MPARAM PARAM, *dev_param;
  //double *dev_diag, *DIAG;
  //double  *dev_kk;
  //  EV TT, *dev_tt;
  /*------------------------------------------------*/
  //alloc_chain( &PSI, &dev_psi );
  set_param( &PARAM, &dev_param ) ;
  MPARAM PBACK;
  cudaMemcpy( &PBACK, dev_param,
                     size_MPARAM, cudaMemcpyDeviceToHost ) ;
  write_param( &PBACK, "pback.txt" );
}
