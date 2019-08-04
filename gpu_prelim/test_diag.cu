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
  double *dev_diag, *DIAG;
  //double  *dev_kk;
  //  EV TT, *dev_tt;
  /*------------------------------------------------*/
  //alloc_chain( &PSI, &dev_psi );
  set_param( &PARAM, &dev_param ) ;
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
  printf( " ok till here, %d \n", size_diag) ;
  double *Bdiag = (double *) malloc( size_diag ) ;
  cudaMemcpy( Bdiag, dev_diag, size_diag, cudaMemcpyDeviceToHost );
  int qdiag = PARAM.qdiag ;
  printf( " device -> host \n") ;
  for (int iN=0; iN<NN; iN++ ){
    for(int iq=0; iq<qdiag; iq++){
      printf( "%f\t%f\n", DIAG[iN*qdiag+iq], Bdiag[iN*qdiag+iq] )  ;
    }
  }
}
