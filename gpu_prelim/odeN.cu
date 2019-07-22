#include <stdlib.h>
#include <stdio.h>
#include "odeN.h"
#include "model.h"
/*=======================================*/
void alloc_odeN( double **PSI, double **psi){
  double *temp;
  *PSI = (double *) malloc( ndim*sizeof( double ) );
  cudaMalloc( (void**)&temp, ndim*sizeof( double ) );
  *psi = temp;
}
/*-------------------------------------------------------------------*/
void  free_odeN( double **PSI, double **psi  ){
  free( *PSI );
  cudaFree( *psi );
}
/*-------------------------------------------------------------------*/
void H2D(double psi[], double PSI[], int Nsize){
  cudaMemcpy( psi, PSI, Nsize*sizeof(double), cudaMemcpyHostToDevice);
}
/*-------------------------------------------------------------------*/
void D2H(double PSI[], double psi[], int Nsize){
  cudaMemcpy( PSI, psi, Nsize*sizeof(double), cudaMemcpyDeviceToHost);
}
/*-------------------------------------------------------------------*/
void iniconf(  double PSI[], double psi[]){
  set_param( );
  initial_configuration( PSI ) ;
  H2D( psi, PSI, ndim );
}

