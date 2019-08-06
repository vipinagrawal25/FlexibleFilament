#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "evolve.h"
double *dev_kk;
double *dev_temp, *dev_k1, *dev_k2, *dev_k3, *dev_k4;
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
  (*TT).ttemp = (*TT).time;
  (*TT).dt = 1.e-4;
  (*TT).ndiag = 2;
  (*TT).tmax =2.e-3;
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
  cudaMalloc( (void**)&dev_temp, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k1, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k2, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k3, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k4, Nsize*sizeof( double ) );
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
  printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag  );
  ALGO( PSI, dev_psi,
          TT, dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag,
          BUG, dev_bug,
          Nblock, Nthread );
  // Once evolution is done, copy the data back to host
  D2H(PSI, dev_psi, ndim);
}
/*----------------------------------------------------------------*/
__global__ void eval_rhs(double kk[], double psi[],
                         EV *tt, MPARAM *param, double diag[],
                         CRASH *bug,  int ldiag ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    model_rhs( kk, psi, tid, (*tt).ttemp, param, diag, crash, ldiag  );
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
}
/*--------------------------- --------------------------------- */
void euler( double PSI, double dev_psi[],
            EV TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[],
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
  /* I do time-marching */
  // copy time parameters to GPU
  cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice);
  while ( TT.time < TT.tmax){
    // copy time parameters to GPU
    int ldiag = 1;
    eval_rhs<<<Nblock,Nthread >>>( dev_kk, dev_psi,  dev_tt ,
                                   dev_param, dev_diag, dev_crash, ldiag  );
    // check if there were any bugs from rhs evaluation
    cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
    if ( BUG.lstop) { IStop( BUG );}
    // take the Euler step
    eustep<<<Nblock,Nthread>>>( dev_kk,  dev_psi, dev_tt ) ;
    cudaMemcpy( &TT, dev_tt, size_EV, cudaMemcpyDeviceToHost);
    printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag) ;
    // reduce diagnostic
    if ( ldiag ) {
      reduce_diag<<<Nblock, Nthread >>> ( dev_diag, ldiag  ) ;  
      // diagnostic copied to host out here
      cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    }
  } // time loop ends here
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
        (*tt).ttemp += (*tt).dt ;
        (*tt).time = (*tt).ttemp ;
      }
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}

/*----------------------------------------------------------------*/
__global__ void rk4one( double temp[], double k1[], double psi[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      temp[ip+pp*tid] = psi[ip+pp*tid] + k1[ip+pp*tid]*(*tt).dt/2 ;
      if (tid == 0) { // only 0th thread updates time
        (*tt).ttemp += (*tt).dt/2 ;
      }
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rnkt4( double PSI, double dev_psi[],
                       EV TT, EV *dev_tt,
                       MPARAM PARAM,   MPARAM *dev_param,
                       double DIAG[], double dev_diag[],
                       CRASH BUG, CRASH *dev_bug,
                       int Nblock, int Nthread ){
  /* I do time-marching */
  // copy time parameters to GPU
  cudaMemcpy( dev_tt, &TT, size_EV, cudaMemcpyHostToDevice);
  while ( TT.time < TT.tmax){
    // first evaluation of rhs, diagnostic is calculated in this step
    int ldiag = 1;
    eval_rhs<<<Nblock,Nthread >>>( dev_kk, dev_psi,  dev_tt ,
                                   dev_param, dev_diag, dev_crash, ldiag  );
    // check if there were any bugs from rhs evaluation
    cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
    if ( BUG.lstop) { IStop( BUG );}
    // take the first substep
    rk4one<<<Nblock,Nthread>>>( dev_temp,  dev_k1, dev_psi, dev_tt ) ;
    // second evaluation of rhs, no diagnostic calculated
    eval_rhs<<<Nblock,Nthread >>>( dev_kk, dev_psi,  dev_tt ,
                                   dev_param, dev_diag, dev_crash, 0  );
    // check if there were any bugs from rhs evaluation
    cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
    if ( BUG.lstop) { IStop( BUG );}
    // take the second substep
    rk4two<<<Nblock,Nthread>>>( dev_temp,  dev_k1, dev_psi, dev_tt ) ;

    // calculate diagnostic and update time
    euler_reduce_decide<<<Nblock, Nthread >>> ( dev_diag, ldiag  ) ;  
    cudaMemcpy( &TT, dev_tt, size_EV, cudaMemcpyDeviceToHost);
    printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", TT.tmax, TT.time, TT.tdiag) ;
   // diagnostic copied to host out here
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);

  
  /* I do time-marching till time reaches the time to calculate 
     the next diagnostic */
  while (tid < NN ){
    while( (*tt).time < (*tt).tdiag ){
   
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
