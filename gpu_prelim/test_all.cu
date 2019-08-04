#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "model.h"
#include "cuda.h"
int main( void ){
  CRASH BUG, *dev_bug ; 
  //double *dev_psi, *PSI;
  MPARAM PARAM, *dev_param;
  double *dev_diag, *DIAG;
  EV TT, *dev_tt;
  /*------------------------------------------------*/
  //alloc_chain( &PSI, &dev_psi );
  set_crash(  &BUG, &dev_bug ); 
  set_param( &PARAM, &dev_param ) ;
  int size_diag = pre_diag( &DIAG , &dev_diag, PARAM );
  pre_evolve(  ndim, "euler" , &TT, &dev_tt  ) ;
  printf( "#  Testing if host device communication works \n" );
  printf( " Testing BUG " );
  printf( " Host: BUG.lstop= %d \n BUG.message = %s \n" ,
          BUG.lstop, BUG.message );
  CRASH BBUG ;
  cudaMemcpy( &BBUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost );
  printf( " Host: BBUG.lstop= %d \n BBUG.message = %s \n" ,
          BBUG.lstop, BBUG.message );
  printf ( " Testing TT \n " ) ;
  EV  TBACK;
  cudaMemcpy( &TBACK, dev_tt,
              size_EV, cudaMemcpyDeviceToHost ) ;
  wevolve( &TBACK, "Tback.txt" ) ;
  printf( " Testing PARAM \n" ) ;
  MPARAM PBACK;
  cudaMemcpy( &PBACK, dev_param,
              size_MPARAM, cudaMemcpyDeviceToHost ) ;
  write_param( &PBACK, "pback.txt" );
  /*  double *Bdiag = (double *) malloc( size_diag ) ;
  cudaMemcpy( Bdiag, dev_diag, size_diag, cudaMemcpyDeviceToHost );
  int qdiag = PARAM.qdiag ;
  for (int iN=0; iN<NN; iN++ ){
    for(int iq=0; iq<qdiag; iq++){
      printf( "%f\t%f\n", DIAG[iN*qdiag+iq], Bdiag[iN*qdiag+iq] )  ;
    }
    }*/
}
