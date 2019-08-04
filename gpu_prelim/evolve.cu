#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "evolve.h"
double *dev_kk;
void (*ALGO)( double [], double [], EV *, MPARAM *, double *, CRASH *);
__global__ void euler( double psi[], double k0[], EV *tt, MPARAM *,
                       double *, CRASH *);
__global__ void rnkt4( double psi[], double k0[], EV *tt, MPARAM *,
                       double *, CRASH *);
void pre_euler( int Nsize );
void pre_rnkt4( int Nsize );
/*-----------------------------------------------------------------------*/
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
  (*TT).dt = 1.e-4;
  (*TT).ndiag = 10;
  (*TT).tmax =1.e-3;
  (*TT).tdiag = (*TT).tmax/((double) (*TT).ndiag) ;
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
  cudaMalloc( (void**)&dev_kk, 5*Nsize*sizeof( double ) );
  printf( "--I have set up auxiliary storage in the device --\n " ) ;
}
/*----------------------------------------*/
void evolve( double PSI[], double dev_psi[], 
             EV TT, EV *dev_tt,
             MPARAM *dev_param ,
             double DIAG[], double dev_diag[], int size_diag,
             CRASH BUG,  CRASH *dev_bug ) {
  cudaDeviceProp *prop;
  int count;
  qdevice( &count, &prop ) ;
  printf( "#- We know device properties \n");
  /* We launch the thread here: */
  int Nthread, Nblock;
  /* If the number of elements in the chain, NN,  is smaller the maximum threads
allowed-per-block then we launch in one way */
  if ( NN < prop[0].maxThreadsPerBlock ) {
    Nthread= 32;
    Nblock = (NN+31)/32 ;
  } else{
    // Otherwise: we launch the threads differently:
    Nthread = 128;
    Nblock = (NN+127)/128;
  }
  printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag  );
  printf( "#-I am launching %d threads in %d blocks\n", Nthread, Nblock );
  /* I evolve till I reach the point where the first diagnostic must
     be calculated.  This evolution is done inside the time-stepper */
  while ( TT.time < TT.tmax){
     // copy time parameters to GPU
    cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice);
    ALGO<<<Nblock,Nthread>>>( dev_psi, dev_kk, dev_tt, dev_param ,
                              dev_diag, dev_bug);
    cudaMemcpy( &TT, dev_tt, size_EV, cudaMemcpyDeviceToHost);
   // set the time for next diagnostic
    TT.tdiag = TT.time+ TT.tmax/(double)TT.ndiag ;
    printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag) ;
   // diagnostic copied to host out here
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
    if ( BUG.lstop) { IStop( BUG );}
  }
  // Once evolution is done, copy the data back to host
  D2H(PSI, dev_psi, ndim);
}
/*--------------------------- --------------------------------- */
__global__ void euler( double psi[], double k0[], EV *tt,
                       MPARAM *dev_param, double diag[], CRASH *crash ){
  //double dpsi[pp];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  /* I do time-marching till time reaches the time to calculate 
     the next diagnostic */
  while (tid < NN ){
    while( (*tt).time < (*tt).tdiag ){
      eval_rhs( k0, psi, tid, (*tt).time, dev_param, diag, crash );
      /*  for (int ip=0; ip<pp; ip++){
        k0[ip+pp*tid]=dpsi[ip];
        }*/
      __syncthreads( ) ;
      /*  Everythings is part of the same chain, so we need to sync the threads.  */ 
      for ( int ip=0; ip<pp; ip++){
        psi[ip+pp*tid] += k0[ip+pp*tid]*(*tt).dt;
      }
   /* then increment time itself, by only one thread
      (This continues to work of the ODEs are not explicitly time dependent ) */   
      if (tid == 0) {
        (*tt).time += (*tt).dt; 
      }
      __syncthreads( );
    }// the while loop finished here.
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
  // then sync the threads again. 
  __syncthreads( );
}
/*-----------------------------------------------------------------------*/
__global__ void rnkt4( double psi[], double kk[], EV *tt,
                       MPARAM *dev_param, double *diag, CRASH *crash ){
  //double dpsi[pp];
  double *temp, *k1, *k2, *k3, *k4;
  temp = &kk[0];
  k1= &kk[ndim];
  k2= &kk[2*ndim];
  k3= &kk[3*ndim];
  k4= &kk[4*ndim];
  int tid = threadIdx.x + blockIdx.x * blockDim.x ;
  /* I do time-marching till time reaches the time to calculate 
     the next diagnostic */
  while (tid < NN ){
    while( (*tt).time < (*tt).tdiag ){
      // step 1 of rnkt4
      eval_rhs( k1, psi, tid, (*tt).time, dev_param, diag, crash );
      /* for (int ip=0; ip<pp; ip++){
        k1[ip+pp*tid]=dpsi[ip];
        }*/
      __syncthreads( );
      for ( int ip=0; ip<pp; ip++){
        temp[ip+pp*tid] = psi[ip+pp*tid] + k1[ip+pp*tid]*(*tt).dt/2 ;
      }
  // Make sure all the threads have done the sub-step. 
      __syncthreads( );
   // step 2 of rnkt4
      eval_rhs( k2, temp, tid, (*tt).time + (*tt).dt/2, dev_param, diag, crash );
      /* for (int ip=0; ip<pp; ip++){
        k2[ip+pp*tid] = dpsi[ip];
        }*/
      // sync the thread here before updating temp
      __syncthreads( );
      for ( int ip=0; ip<pp; ip++){
        temp[ip+pp*tid] = psi[ip+pp*tid] + k2[ip+pp*tid]*(*tt).dt/2 ;
      }
  // Make sure all the threads have done the sub-step. 
      __syncthreads( );
  // step 3 of rnkt4
      eval_rhs( k3, temp, tid, (*tt).time + (*tt).dt/2, dev_param, diag, crash );
      /* for (int ip=0; ip<pp; ip++){
        k3[ip+pp*tid] = dpsi[ip];
        }*/
      // sync the thread here before updating temp
      __syncthreads( );
      for ( int ip=0; ip<pp; ip++){
        temp[ip+pp*tid] = psi[ip+pp*tid] + k3[ip+pp*tid]*(*tt).dt/2 ;
      }
  // Make sure all the threads have done the sub-step. 
      __syncthreads( );     
   // step 4 of rnkt4
      eval_rhs( k4, temp, tid, (*tt).time + (*tt).dt, dev_param, diag, crash );
      /* for (int ip=0; ip<pp; ip++){
        k4[ip+pp*tid] = dpsi[ip];
        }*/
      __syncthreads( );
      for ( int ip=0; ip<pp; ip++){
        psi[ip+pp*tid] +=  (k1[ip+pp*tid]/6. + k2[ip+pp*tid]/3.
                                      + k3[ip+pp*tid]/3. + k4[ip+pp*tid]/6. )*(*tt).dt ;
      }
  // then increment time itself, by only one thread
      if (tid == 0) {
        (*tt).time += (*tt).dt; 
      }
    }// the while loop finished here.
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
  // then sync the threads again. 
  __syncthreads( );
}
