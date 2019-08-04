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
void IStop( CRASH BUG ){
  printf( "#-I STOP, something went wrong \n") ;
  printf( "#-%s \n", BUG.message );
  exit(1);
}
/*-----------------------------------------------------------------------------*/
void set_crash( CRASH *BUG, CRASH **dev_bug ){
  (*BUG).lstop = 0;
  strcpy( (*BUG).message, " No bug yet" );
  CRASH *temp ;
  cudaMalloc( (void**)&temp, size_CRASH );
  *dev_bug = temp;
  cudaMemcpy( *dev_bug, BUG,
                     size_CRASH, cudaMemcpyHostToDevice ) ;
} 
/*-----------------------------------------------------------------------------*/
