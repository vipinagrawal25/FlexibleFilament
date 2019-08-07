#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "evolve.h"
double *dev_kk;
double *dev_psip, *dev_k1, *dev_k2, *dev_k3, *dev_k4;
void (*ALGO)( double [], double [],
              EV , EV *,
              MPARAM ,   MPARAM *,
              double [], double [], int,  
              CRASH , CRASH *,
              int , int  );
void euler( double PSI[], double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], int ldiag,  
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread );
void rnkt4( double PSI[], double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], int ldiag,
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread );
void pre_euler( int Nsize );
void pre_rnkt4( int Nsize );
__global__ void reduce_diag( double diag[] );
__global__ void eustep( double kk[], double psi[], EV *tt );
/*-------------------------------------------------------------------*/
void  post_euler( void  ){
  cudaFree( dev_kk );
}
/*-------------------------------------------------------------------*/
void  post_rnkt4( void  ){
  cudaFree( dev_psip );
  cudaFree( dev_k1 );
  cudaFree( dev_k2 ) ;
  cudaFree( dev_k3 ) ;
  cudaFree( dev_k4 );
}
/*-----------------------------------------------------------------------*/
void post_evolve( char *algo  ){
  if ( strcmp( algo , "euler") == 0 ){
    post_euler( );
  } else if ( strcmp( algo, "rnkt4" ) == 0 ){
    post_rnkt4(  );
  } else {
    printf( " algorithm\t%s\t not coded \n", algo);
    printf( "EXITING \n " );
    exit(1);
  }
}
/*----------------------------------------------------------------------------------*/
void pre_evolve( int Nsize, char *algo, EV *TT,  EV **dev_tt ){
  if ( strcmp( algo , "euler") == 0 ){
    pre_euler( Nsize );
    ALGO = &euler;
  } else if ( strcmp( algo, "rnkt4" ) == 0 ){
    pre_rnkt4( Nsize );
    ALGO = &rnkt4;
  } else {
    printf( " algorithm\t%s\t not coded \n", algo);
    printf( "EXITING \n " );
    exit(1);
  }
  /* Set up parameters for evolution */
  (*TT).time = 0.;
  (*TT).tprime = (*TT).time;
  (*TT).dt = 1.e-4;
  (*TT).ndiag = 2;
  (*TT).tmax =2.e-3;
  (*TT).tdiag = 0. ;
  (*TT).substep = 0 ;
  EV *temp ;
  cudaMalloc(  (void**)&temp, size_EV );
  *dev_tt = temp;
  cudaMemcpy( *dev_tt, TT,
                     size_EV, cudaMemcpyHostToDevice ) ;
  wevolve( TT, "initial.txt" );
}
/* ------------------------------------------------------------------------------*/
void wevolve( EV *TT, char *fname ){
  FILE *pout ;
  pout = fopen ( fname, "w" );
  fprintf( pout, "# =========== Evolve Parameters==========\n" );
  fprintf( pout, " time=%f \n ", (*TT).time ) ;
  fprintf( pout, " tprime=%f \n ", (*TT).tprime ) ;
  fprintf( pout, "dt = %f \n ", (*TT).dt );
  fprintf( pout, " ndiag = %d \n ", (*TT).ndiag );
  fprintf( pout, "tmax = %f \n ", (*TT).tmax );
  fprintf( pout, "tdiag = %f \n ", (*TT).tdiag );
}
  /*----------------------------------------*/
void pre_euler( int Nsize ){
  printf( " #---time-integration algorithm : EULER --\n " );
  cudaMalloc( (void**)&dev_kk, Nsize*sizeof( double ) );
  printf( "#--I have set up auxiliary storage in the device-- \n " ) ;
}
/*----------------------------------------*/
void pre_rnkt4( int Nsize ){
  printf( "#-- time-integration algorithm : RNKT4-- \n " );
  cudaMalloc( (void**)&dev_psip, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k1, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k2, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k3, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k4, Nsize*sizeof( double ) );
  printf( "--I have set up auxiliary storage in the device --\n " ) ;
}
/*----------------------------------------*/
void evolve( double PSI[], double dev_psi[], 
             EV TT, EV *dev_tt,
             MPARAM PARAM, MPARAM *dev_param ,
             double DIAG[], double dev_diag[], int size_diag,
             CRASH BUG,  CRASH *dev_bug ) {
  cudaDeviceProp *prop;
  int count;
  qdevice( &count, &prop ) ;
  printf( "#- We know device properties \n");
  /* We set the number threads here: */
  int Nthread, Nblock;
  /* If the number of elements in the chain, NN,  is smaller the maximum threads
allowed-per-block then we launch in one way */
  if ( NN < 128 ) {
    Nthread= 1;
    Nblock = NN ;
  } else{
    // Otherwise: we launch the threads differently:
    Nthread = 128;
    Nblock = (NN+127)/128;
  }
  printf( "#-I shall launch %d blocks, each with %d threads\n", Nblock, Nthread );
  while ( TT.time < TT.tmax){
    printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag) ;
    int ldiag = 0 ;
    if( TT.time >= TT.tdiag ) {
      ldiag = 1;
      TT.tdiag = TT.time +  TT.tmax/((double) TT.ndiag) ;
    }
    cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice);
    ALGO( PSI, dev_psi,
          TT, dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag, ldiag, 
          BUG, dev_bug,
          Nblock, Nthread );
    cudaMemcpy( &TT, dev_tt, size_EV, cudaMemcpyDeviceToHost);
  } // while time loop ends here
}
/*----------------------------------------------------------------*/
__global__ void eval_rhs(double kk[], double psi[],
                         EV *tt, MPARAM *param, double diag[],
                         CRASH *bug,  int ldiag ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    model_rhs( kk, psi, tid, (*tt).tprime, param, diag, bug, ldiag  );
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
}
/*--------------------------- --------------------------------- */
void euler( double PSI[], double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], int ldiag,  
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
  /* I do time-marching */
  // copy time parameters to GPU
  eval_rhs<<<Nblock,Nthread >>>( dev_kk, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, ldiag  );
    // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // reduce diagnostic
  if ( ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, TT.time, PARAM );
  }
  // take the Euler step
  eustep<<< Nblock, Nthread >>>( dev_psi, dev_kk, dev_tt ) ;
}
  /*----------------------------------------------------------------*/
__global__ void reduce_diag( double diag[] ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    tid += blockDim.x * gridDim.x ;
  }
}
/*----------------------------------------------------------------*/
__global__ void eustep( double kk[], double psi[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psi[ip+pp*tid] += kk[ip+pp*tid]*(*tt).dt ;
      if (tid == 0) { // only 0th thread updates time
        (*tt).time += (*tt).dt ;
        (*tt).tprime = (*tt).time ;
      }
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rk4_time_substep( EV TT, EV *dev_tt, int j ){
  double rk4a[4] ;
  rk4a[0] = 1./2. ;
  rk4a[1] = 1./2. ;
  rk4a[2] = 1. ;
  rk4a[3] = 0. ;
  TT.tprime = TT.time + (TT.dt)*rk4a[ j ] ;
  TT.substep = j + 1;
  cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice ) ;
}
/*----------------------------------------------------------------*/
__global__ void rk4_psi_substep( double psip[], double kin[], double psi[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  double rk4a[4] ;
  rk4a[0] = 1./2. ;
  rk4a[1] = 1./2. ;
  rk4a[2] = 1. ;
  rk4a[3] = 0. ;
  int j = (*tt).substep ;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psip[ip+pp*tid] = psi[ip+pp*tid] + kin[ip+pp*tid]*(*tt).dt*rk4a[ j ] ;
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rk4_time_step( EV TT, EV *dev_tt ){
  TT.time = TT.time + TT.dt ;
  TT.tprime = TT.time ;
  TT.substep = 0;
  cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice ) ;
}
/*----------------------------------------------------------------*/
__global__ void rk4_psi_step( double psi[],
                          double k1[], double k2[],
                          double k3[], double k4[],
                          EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psi[ip+pp*tid] = psi[ip+pp*tid] +
        (k1[ip+pp*tid]/6.+ k2[ip+pp*tid]/3. + k3[ip+pp*tid]/3. + k4[ip+pp*tid]/6. )*(*tt).dt ;
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rnkt4( double PSI[], double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], int ldiag, 
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
  /* I do time-marching */
  // 1st evaluation of rhs, diagnostic is calculated in this step
  eval_rhs<<<Nblock,Nthread >>>( dev_k1, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, ldiag  );
    // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
 // reduce diagnostic
  if ( ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    printf( "calculating diagnostics \n " );
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, TT.time, PARAM );
  }
  // take the first substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k1, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 0 ) ;
  // 2nd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k2, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the second substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k2, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 1 ) ;
  // 3rd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k3, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the third substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k3, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 2 ) ;
  // 4th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k4, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // final step
    rk4_psi_step<<<Nblock, Nthread >>>( dev_psi,
                                    dev_k1, dev_k2, dev_k3, dev_k4, dev_tt );
    rk4_time_step( TT, dev_tt ) ;
}

/*----------------------------------------------------------------*/
__global__ void rkf45_psi_substep( double psip[], double kin[], double psi[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  double rk4a[4] ;
  rk4a[0] = 1./2. ;
  rk4a[1] = 1./2. ;
  rk4a[2] = 1. ;
  rk4a[3] = 0. ;
  int j = (*tt).substep ;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psip[ip+pp*tid] = psi[ip+pp*tid] + kin[ip+pp*tid]*(*tt).dt*rk4a[ j ] ;
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rnkf45( double PSI[], double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], int ldiag, 
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
  /* I do time-marching */
  // 1st evaluation of rhs, diagnostic is calculated in this step
  eval_rhs<<<Nblock,Nthread >>>( dev_k1, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, ldiag  );
    // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
 // reduce diagnostic
  if ( ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    printf( "calculating diagnostics \n " );
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, TT.time, PARAM );
  }
  // take the first substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k1, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 0 ) ;
  // 2nd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k2, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the second substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k2, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 1 ) ;
  // 3rd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k3, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the third substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k3, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 2 ) ;
  // 4th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k4, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug, 0  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // final step
    rk4_psi_step<<<Nblock, Nthread >>>( dev_psi,
                                    dev_k1, dev_k2, dev_k3, dev_k4, dev_tt );
    rk4_time_step( TT, dev_tt ) ;
}


