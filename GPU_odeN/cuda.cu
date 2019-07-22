#include <stdio.h>
#include "cuda.h"
/*========================================*/
void  qdevice(int *count, cudaDeviceProp **prop ) {
  int ndevice;
  cudaGetDeviceCount( &ndevice ) ;
  *prop = (cudaDeviceProp *) malloc( ndevice*sizeof( cudaDeviceProp ) );
  for (int i=0; i<ndevice; i++){
    cudaGetDeviceProperties( prop[i], i);
  }
  *count = ndevice;
}
  /*----------------------------------------------*/
  void qfree( cudaDeviceProp *prop ){
    free( prop );
  }
