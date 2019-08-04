#include <stdio.h>
#include "3vec.h"
#include "2Tens.h"
#define NN 8
#define pp 3
__global__ void step( double xx[], double dproj[] ); 
__device__ void smallstep( double xx[], int tid, double dproj[] );
/*=====================================*/
int main(){
  double *x, *dev_x ;
  double *proj, *dev_proj;
  int ndim = NN*pp ; 
  x = (double *) malloc( ndim*sizeof( double ) );
  cudaMalloc( (void**)&dev_x, ndim*sizeof( double ) );
  for (int i = 0; i<ndim; i++){
    x[i] = 1.+(double) i ;
  }
  cudaMemcpy( dev_x, x, ndim*sizeof(double), cudaMemcpyHostToDevice);
  //
    proj = (double *) malloc( NN*sizeof( double ) );
  cudaMalloc( (void**)&dev_proj, NN*sizeof( double ) );
  
  int Nblock = NN;
  int Nthread = 1;
  step<<<Nblock, Nthread>>>( dev_x , dev_proj );
  cudaMemcpy( x, dev_x, ndim*sizeof(double), cudaMemcpyDeviceToHost ) ;
  cudaMemcpy( proj, dev_proj, NN*sizeof(double), cudaMemcpyDeviceToHost ) ;
  for (int i = 0; i<NN; i++){
    for (int ip=0; ip<pp; ip++){
      printf( "x = %f \t ", x[pp*i + ip] );
    }
      printf( " proj = %f  \n ", proj[i] ) ;
  }
}
/*-------------------------------*/
__global__ void step( double xx[], double dproj[] ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x ;
  smallstep( xx, tid, dproj ) ;
 }
__device__ void smallstep( double xx[], int tid, double dproj[] ){
   vec3 y(1. , 0. , 0. );
  Tens2 dab(1.,0.,0.,0.,1.,0.,0.,0.,1.);
  while (tid < NN ){
  vec3 x = vec3( xx[tid*pp], xx[tid*pp+1], xx[tid*pp+2] );
   y = y + dot(dab, x );
   dproj[ tid ] = dot( y, x) ;
   xx[tid*pp ] = y.x ;
   xx[tid*pp + 1] = y.y ;
   xx[tid*pp + 2] = y.z ; 
  tid += blockDim.x * gridDim.x ;
  }
}
