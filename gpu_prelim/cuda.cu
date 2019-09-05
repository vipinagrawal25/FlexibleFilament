#include <stdio.h>
#include "cuda.h"
#include "model.h"
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
/*-----------------------------------------------------------------------*/
__global__ void thread_maxima( double array[], double redux[] ){
  // This is to calculate maxima of the operations going on in a thread.
  //After this the maxima should be compared among all the blocks as well.
  extern __shared__ double cache[];  // Just to create a shared variable
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  int i=blockDim.x/2;
  int cacheIndex=threadIdx.x;
  // cache[cacheIndex] = array[tid];
  // __syncthreads();
  double temp=0.;
  while(tid<NN){
    temp = max(temp,array[tid]);
    tid += blockDim.x*gridDim.x;
  }
  cache[cacheIndex]=temp;
  __syncthreads();
  while(i!=0){
    if (cacheIndex<i){
      if (cache[cacheIndex+i]>cache[cacheIndex]){
        cache[cacheIndex]=cache[cacheIndex+i];
      }
      // cache[cacheIndex]=max(cache[cacheIndex],cache[cacheIndex+i]);
      // cache[cacheIndex]=i;
    }
    i/=2;
    __syncthreads();
  }
  if (cacheIndex==0){
    redux[blockIdx.x] = cache[0];
  }
}
/*-----------------------------------------------------------------------*/
__global__ void thread_sum( double array[], double redux[] ){
  // This is to calculate maxima of the operations going on in a thread.
  //After this the maxima should be compared among all the blocks as well.
  extern __shared__ double cache[];  // Just to create a shared variable
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  int i=blockDim.x/2;
  int cacheIndex=threadIdx.x;
  // cache[cacheIndex] = array[tid];
  // __syncthreads();
  double temp=0.;
  while(tid<NN){
    temp+=array[tid];
    tid += blockDim.x*gridDim.x;
  }
  cache[cacheIndex]=temp;
  __syncthreads();
  while(i!=0){
    if (cacheIndex<i){
      cache[cacheIndex]+=cache[cacheIndex+i];
      // cache[cacheIndex]=max(cache[cacheIndex],cache[cacheIndex+i]);
      // cache[cacheIndex]=i;
    }
    i/=2;
    __syncthreads();
  }
  if (cacheIndex==0){
    redux[blockIdx.x] = cache[0];
  }
}