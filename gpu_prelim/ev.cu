#include <stdio.h>
#include "sub.h"
double *dev_b;
void preAdd( void ){
    cudaMalloc( (void**)&dev_b, ndim*sizeof( double ) );
     double  *temp = (double *) malloc( ndim*sizeof( double ) );
     for( int i = 0; i < ndim; i++){
       temp[i] = 0.;
     }
//     cudaMemcpy( dev_b,  temp, ndim * sizeof(double),
//                 cudaMemcpyHostToDevice);
//     return dev_b;
}
/*---------------------------------------------*/
__global__ void ADD( double a[], double b[] ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < ndim ){
    a[ tid ] = a[ tid ] + b[ tid ] ; 
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
}
/*---------------------------------------------*/
__global__ void setB( double a[], double b[] ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < ndim ){
    b[ tid ] = a[ tid ] ; 
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
}
/*---------------------------------------------*/
//void   AddB( double A[], double *dev_a, double *dev_b ) {
void   AddB( double A[], double *dev_a ) {
  int Nblock = 4 ;
  int Nthread = (ndim + 3)/4;
  setB<<< Nblock,  Nthread >>> ( dev_a, dev_b ) ;
  ADD<<< Nblock,  Nthread >>> ( dev_a, dev_b ) ;
  cudaMemcpy( A, dev_a, ndim*sizeof( double ),  cudaMemcpyDeviceToHost);
}
