#include <stdio.h>
#include<stdlib.h>
#include <math.h>
//#include "evolve.h"
//#include "chain.h"
#include "cuda.h"
#define   nn 10
struct DM{
  double x;
};
int size_DM = sizeof(double);
/* ----------------------------------------------------- */
__device__ void increment(double *psi, DM *a, int kth ){
  psi[kth] = psi[kth] + (*a).x ;
}
__global__ void kernel( double *psi, DM *a){
   int tid = threadIdx.x + blockIdx.x * blockDim.x ;
   if ( tid < nn ){
     increment( psi, a, tid ) ;
   }
}
/* ----------------------------------------------------- */
struct EV  TT;
void test_struct( void );
int main( void ){
  cudaDeviceProp *prop;
  double PSI[nn] ;
  double *dev_psi;
  // double X;
  //double *dev_x;
  DM A;
  DM *dev_a ;
  printf( "psi[i]= \t");
  for(int i=0; i<nn; i++) {
    PSI[i] = (double) i;
    printf( "%f\t", PSI[i] );
  }
  printf( "\n" );
  cudaMalloc( (void**)&dev_psi, nn*sizeof( double ) );
  cudaMemcpy( dev_psi, &PSI, nn*sizeof(double), cudaMemcpyHostToDevice );
  A.x = 10.;
  cudaMalloc( (void**)&dev_a, size_DM );
  cudaMemcpy( dev_a, &A, size_DM, cudaMemcpyHostToDevice );
  kernel<<<nn,1>>> (dev_psi, dev_a) ;
  cudaMemcpy( &PSI, dev_psi, nn*sizeof(double), cudaMemcpyDeviceToHost );
  printf( "psi[i]= \t");
  for(int i=0; i<nn; i++) {
    printf( "%f\t", PSI[i] );
  }
  printf( "\n" );
}
/*---------------------------------*/

