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
/*----------------------------------------------*/
__device__ void scpy( char to[], char from[]){
  int i=0;
  while ( (to[i] = from[i]) != '\0')
    i = i+1;
}
/*----------------------------------------------*/
__device__ void device_exception( struct CRASH *bug, char *mesg ){
        (*bug).lstop = 1;
        scpy( (*bug).message, mesg ) ;
}
/*-----------------------------------------------------------------------------*/
