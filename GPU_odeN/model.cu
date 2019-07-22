#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "odeN.h"
#include "model.h"
__host__ void write_param( void );
struct MPARAM host_param;
struct MPARAM *dev_param;
/* ========================================= */
void set_param( void ){
  host_param.OM = 1. ;
  cudaMemcpy( dev_param, &host_param,
                     size_MPARAM, cudaMemcpyHostToDevice ) ;
  write_param( );
}
/* ===================================== */
void write_param( void ){
  printf( "# =========== Model Parameters ==========\n" );
  printf( " #Model : simple harmonic oscillator \n " ) ;
  printf( "#dimension of ODE:\n pp =  %d \n", pp ) ;
  printf( "#Number of copies:\n  NN = %d\n", NN ) ;
  printf( " omega = %f \n", host_param.OM );
  printf( " #============================\n" );
}
/* ===================================== */
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          struct MPARAM *dev_param){
  double xx, vv;
  // double OM = (*dev_param).OM ;
  double OM=0.;
  // We solve for 1-d harmonic oscillator.
  xx = psi[kelement];
  vv = psi[kelement+1];
  dpsi[kelement] = vv;
  dpsi[kelement+1] = - OM*OM*xx ;
}
/*-------------------------------------------------------------------*/
__host__ void initial_configuration( double PSI[] ){
  double xx, vv;
  for( int s=0; s<NN; s++){
   // A Harmonic Oscillator 
 /* initially all positions are zero */
    xx = 0.;
 /* and all velocities are unity */
    vv = 1.;
    PSI [pp*s ] = xx;
    PSI[pp*s+1] = vv;
  } 
}
/*-------------------------------------------------------------------*/
