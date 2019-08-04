#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "model.h"
#include "cuda.h"
int main( void ){
    EV TT, *dev_tt;
    /*------------------------------------------------*/
  //alloc_chain( &PSI, &dev_psi );
    pre_evolve(  ndim, "euler" , &TT, &dev_tt  ) ;
     EV  TBACK;
  cudaMemcpy( &TBACK, dev_tt,
                     size_EV, cudaMemcpyDeviceToHost ) ;
  wevolve( &TBACK, "Tback.txt" );
}
