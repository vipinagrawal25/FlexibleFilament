#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "evolve.h"
#include <math.h>
#define tiny pow(10,-15)
double *dev_kk;
double *dev_psip, *dev_k1, *dev_k2, *dev_k3, *dev_k4, *dev_k5, *dev_k6, **dev_kin;
double *dev_err;
double *dev_redux, *REDUX; 
int size_redux;

void (*ALGO)( double [], double [],
              EV* , EV *,
              MPARAM ,   MPARAM *,
              double [], double [],  
              CRASH , CRASH *,
              int , int  );
void euler( double PSI[], double dev_psi[],
            EV* TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[],  
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread );
void rnkt4( double PSI[], double dev_psi[],
            EV* TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[],
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread );
void rnkf45( double PSI[], double dev_psi[],
            EV* TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[],
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread );
void pre_euler( int Nsize);
void pre_rnkt4( int Nsize);
void pre_rnkf45( int Nsize, int Nblock, int Nthread );
void eu_time_step( EV* TT, EV *dev_tt );
__global__ void reduce_diag( double diag[] );
__global__ void eu_psi_step( double psi[], double kk[], EV *tt );
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
    post_rnkt4( );
  } else {
    printf( " algorithm\t%s\t not coded \n", algo);
    printf( "EXITING \n " );
    exit(1);
  }
}
/*----------------------------------------------------------------------------------*/
void pre_evolve( int Nsize, char *algo, EV *TT,  EV **dev_tt, int Nblock, int Nthread ){
  if ( strcmp( algo , "euler") == 0 ){
    pre_euler( Nsize );
    ALGO = &euler;
  } else if ( strcmp( algo, "rnkt4" ) == 0 ){
    pre_rnkt4( Nsize );
    ALGO = &rnkt4;
  } else if (strcmp( algo, "rnkf45" ) == 0){
    pre_rnkf45( Nsize, Nblock, Nthread );
    ALGO = &rnkf45;
  }
 else {
    printf( " algorithm\t%s\t not coded \n", algo);
    printf( "EXITING \n " );
    exit(1);
  }
  /* Set up parameters for evolution */
  (*TT).time = 0.;
  (*TT).tprime = (*TT).time;
  (*TT).dt = 1.e-4;
  (*TT).ndiag = 2;
  (*TT).tmax =2.e-4;
  (*TT).tdiag = 0. ;
  (*TT).substep = 0 ;
  EV *temp ;
  cudaMalloc(  (void**)&temp, size_EV );
  *dev_tt = temp;
  cudaMemcpy( *dev_tt, TT,
                     size_EV, cudaMemcpyHostToDevice ) ;
  wevolve( *TT, "initial.txt" );
}
/* ------------------------------------------------------------------------------*/
void wevolve( EV TT, char *fname ){
  FILE *pout ;
  pout = fopen ( fname, "w" );
  fprintf( pout, "# =========== Evolve Parameters==========\n" );
  fprintf( pout, " time=%f \n ", TT.time ) ;
  fprintf( pout, " tprime=%f \n ", TT.tprime ) ;
  fprintf( pout, "dt = %f \n ", TT.dt );
  fprintf( pout, " ndiag = %d \n ", TT.ndiag );
  fprintf( pout, "tmax = %f \n ", TT.tmax );
  fprintf( pout, "tdiag = %f \n ", TT.tdiag );
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
void pre_rnkf45( int Nsize, int Nblock, int Nthread ){
  printf( "#-- time-integration algorithm : RNKF45-- \n " );
  cudaMalloc( (void**)&dev_psip, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k1, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k2, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k3, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k4, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k5, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_k6, Nsize*sizeof( double ) );
  cudaMalloc( (void**)&dev_err, NN*sizeof(double));
  cudaMalloc( (void**)&dev_kin, 6*sizeof(double*));

  // Defining array of all k pointers.
  double *KIN[6] = {dev_k1,dev_k2,dev_k3,dev_k4,dev_k5,dev_k6};
  cudaMemcpy(dev_kin,KIN,6*sizeof(double*),cudaMemcpyHostToDevice);

  size_redux = Nblock*sizeof(double);
  cudaMalloc( (void**)&dev_redux, size_redux );
  REDUX = (double *) malloc(size_redux);
  for (int iblock = 0; iblock < Nblock; ++iblock){
    REDUX[iblock] = 0.;
  }
  cudaMemcpy(dev_redux,REDUX,size_redux,cudaMemcpyHostToDevice);
  printf( "--I have set up auxiliary storage in the device --\n " ) ;
}
/*----------------------------------------*/
void evolve( double PSI[], double dev_psi[], 
             EV *TT, EV *dev_tt,
             MPARAM PARAM, MPARAM *dev_param ,
             double DIAG[], double dev_diag[], int size_diag,
             CRASH BUG,  CRASH *dev_bug,
             int Nblock, int Nthread ) {
  cudaDeviceProp *prop;
  int count;
  qdevice( &count, &prop ) ;
  printf( "#- We know device properties \n");
// copy the time data to device. 
  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice);
  while ( (*TT).time < (*TT).tmax){
    printf( "#- tmax=%f\t time=%f\t tdiag=%f \n", (*TT).tmax, (*TT).time, (*TT).tdiag) ;
      (*TT).ldiag = 0 ;
    if( (*TT).time >= (*TT).tdiag ) {
      (*TT).ldiag = 1;
      (*TT).tdiag = (*TT).time +  (*TT).tmax/((double) (*TT).ndiag) ;
    }
    ALGO( PSI, dev_psi,
          TT, dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag,
          BUG, dev_bug,
          Nblock, Nthread );
  } // while time loop ends here
}
/*----------------------------------------------------------------*/
__global__ void eval_rhs(double kk[], double psi[],
                         EV *tt, MPARAM *param, double diag[],
                         CRASH *bug ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    model_rhs( kk, psi, tid, (*tt).tprime, param, diag, bug, (*tt).ldiag  );
    tid += blockDim.x * gridDim.x ;
  }// while loop over threads finishes here.
}
/*--------------------------- --------------------------------- */
void euler( double PSI[], double dev_psi[],
            EV *TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[],  
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
// evaluate the right hand side in a kernel 
  eval_rhs<<<Nblock,Nthread >>>( dev_kk, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
    // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // reduce diagnostic
  if ( (*TT).ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, (*TT).time, PARAM );
  }
  // take the Euler step
  eu_psi_step<<< Nblock, Nthread >>>( dev_psi, dev_kk, dev_tt ) ;
  eu_time_step(TT, dev_tt);
}
  /*----------------------------------------------------------------*/
__global__ void reduce_diag( double diag[] ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    tid += blockDim.x * gridDim.x ;
  }
}
/*----------------------------------------------------------------*/
__global__ void eu_psi_step( double psi[], double kk[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psi[ip+pp*tid] += kk[ip+pp*tid]*(*tt).dt ;
    }
    tid += blockDim.x * gridDim.x ;

  } // loop over threads ends here
}
/*----------------------------------------------------------------*/
void eu_time_step( EV *TT, EV *dev_tt ){
  (*TT).time = (*TT).time + (*TT).dt ;
  (*TT).tprime = (*TT).time ;
  (*TT).substep = 0;
  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice ) ;
}
/*-----------------------------------------------------------------------*/
void rk4_time_substep( EV *TT, EV *dev_tt, int j ){
  double rk4a[4] ;
  rk4a[0] = 1./2. ;
  rk4a[1] = 1./2. ;
  rk4a[2] = 1. ;
  rk4a[3] = 0. ;
  (*TT).tprime = (*TT).time + ((*TT).dt)*rk4a[ j ] ;
  (*TT).substep = j + 1;
  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice ) ;
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
void rk4_time_step( EV *TT, EV *dev_tt ){
  (*TT).time = (*TT).time + (*TT).dt ;
  (*TT).tprime = (*TT).time ;
  (*TT).substep = 0;
  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice ) ;
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
            EV* TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], 
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){
  /* I do time-marching */
  // 1st evaluation of rhs, diagnostic is calculated in this step
  eval_rhs<<<Nblock,Nthread >>>( dev_k1, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
    // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
 // reduce diagnostic
  if ( (*TT).ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    printf( "calculating diagnostics \n " );
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, (*TT).time, PARAM );
  }
  // take the first substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k1, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 0 ) ;
  // 2nd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k2, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug);
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the second substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k2, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 1 ) ;
  // 3rd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k3, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug);
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the third substep
  rk4_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_k3, dev_psi, dev_tt ) ;
  rk4_time_substep( TT, dev_tt , 2 ) ;
  // 4th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k4, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // final step
    rk4_psi_step<<<Nblock, Nthread >>>( dev_psi,
                                    dev_k1, dev_k2, dev_k3, dev_k4, dev_tt );
    rk4_time_step( TT, dev_tt ) ;
}
/*----------------------------------------------------------------*/
__global__ void rnkf45_psi_substep( double psip[], double* kin[], double psi[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  double rnkf45b[5][5] = {{0.2,0.,0.,0.,0.},
                        {3./40,9./40,0.,0.,0.},
                        {0.3,-0.9,1.2,0.,0.},
                        {-11./54,2.5,-70./27,35./27,0},
                        {1631./55296,175./512,575./13824,44275./110592,253./4096} };

  int j = (*tt).substep ;

  // double *kk = {&kin[0],&kin[ndim],&kin[2*ndim],&kin[3*ndim],&kin[4*ndim],&kin[5*ndim]}

  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
        psip[ip+pp*tid] = psi[ip+pp*tid];
        for (int jp = 0; jp < j+1; ++jp)
        {
          psip[ip+pp*tid] += kin[jp][ip+pp*tid]*rnkf45b[j][jp];
        }
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
void rnkf45_time_substep( EV* TT, EV *dev_tt, int j ){
  double rnkf45a[5] = {0.2,0.3,0.6,1.,7./8} ;
  (*TT).tprime = (*TT).time + ((*TT).dt)*rnkf45a[ j ] ;
  (*TT).substep = j + 1;
  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice ) ;
}
/*----------------------------------------------------------------*/
__global__ void rnkf45_psi_step( double psi[], double psip[], double error[], double *kptr[], EV *tt ){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  //
  double rnkf45c[2][6] = { {37./378,0,250./621,125./594,0,512./1771},
                          {2825./27648,0.,18575./48384,13525./55296,277./14336,0.25} };
  //
  error[tid]=0;                        
  while (tid < NN ){
    for ( int ip=0; ip<pp; ip++){
      psip[ip+pp*tid] = psi[ip+pp*tid];
      for (int jp = 0; jp < 6; ++jp){
          psi[ip+pp*tid] += (rnkf45c[0][jp]*kptr[jp][ip+pp*tid])*(*tt).dt;
          psip[ip+pp*tid] += (rnkf45c[1][jp]*kptr[jp][ip+pp*tid])*(*tt).dt;   
      }
      error[tid] = max(error[tid],(psi[ip+pp*tid]-psip[ip+pp*tid])*(psi[ip+pp*tid]-psip[ip+pp*tid]));
    }
    tid += blockDim.x * gridDim.x ;
  } // loop over threads ends here
}
/*-----------------------------------------------------------------------*/
bool rnkf45_time_step( EV* TT, EV *dev_tt, double maxErr){
  // cudaMemcpy(  &TT, dev_tt, size_EV, cudaMemcpyDeviceToHost ) ;
  double tol = 0.0001;
  double truncationmax=10;  // Maximum multiplication in time step
  double truncationmin=0.1; // Minimum multiplication in time step
  bool laccept;
  double s;

  if (maxErr<tol){
    laccept=1;
    (*TT).time = (*TT).time + (*TT).dt;
    (*TT).tprime = (*TT).time ;
    (*TT).substep = 0;
    s = (pow(tol/(maxErr),0.2));
    if (s>truncationmax){ s=truncationmax; }
    (*TT).dt=s*((*TT).dt);
  }
  else {
    laccept=0;  
    (*TT).tprime = (*TT).time ;
    (*TT).substep = 0;
    s = pow(tol/maxErr,0.2);
    if (s<truncationmin){s=truncationmin;}
    (*TT).dt = s*((*TT).dt);
    (*TT).ldiag = 0;
  }

  cudaMemcpy( dev_tt, TT, size_EV, cudaMemcpyHostToDevice );
  return laccept;
}
/*-----------------------------------------------------------------------*/
__global__ void thread_maxima( double array[], double redux[]){
 // This is to calculate maxima of the operations going on in a thread. After this the maxima should be compared among all the blocks as well.
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
      cache[cacheIndex]=max(cache[cacheIndex],cache[cacheIndex+i]);
      // cache[cacheIndex]=i;
    }
    i/=2;
    __syncthreads();
  }
  if (cacheIndex==0){
    redux[blockIdx.x] = cache[0];
  }
}
/*----------------------------------------------------------------*/
double MaxDevArray(double dev_array[], int Nblock, int Nthread){
  double maxA=0.;
  // cout << dev_array[5] << endl;
  thread_maxima <<<Nblock, Nthread, Nthread*sizeof(double) >>>(dev_array,dev_redux);
  cudaMemcpy(REDUX,dev_redux,size_redux,cudaMemcpyDeviceToHost);
  // Copied the thread maxima output back to host.
  // Compare the maximum across the blocks. This operation is done in CPU for the time being.
  // Calculate maxima using STL
  // This could be done better by launching a kernel 
  for (int iblock = 0; iblock < Nblock; ++iblock){
    maxA = max(maxA,REDUX[iblock]);
    // cout << REDUX[iblock] << "\t" ;
  }
  // cout << maxA << endl;
  return maxA;
}
/*-----------------------------------------------------------------------*/
void rnkf45( double PSI[], double dev_psi[],
            EV* TT, EV *dev_tt,
            MPARAM PARAM,   MPARAM *dev_param,
            double DIAG[], double dev_diag[], 
            CRASH BUG, CRASH *dev_bug,
            int Nblock, int Nthread ){

  /* I do time-marching */
  // 1st evaluation of rhs, diagnostic is calculated in this step
  eval_rhs<<<Nblock,Nthread >>>( dev_k1, dev_psi,  dev_tt ,
                                 dev_param, dev_diag, dev_bug  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // reduce diagnostic
  if ( (*TT).ldiag ) {
    //reduce_diag<<<Nblock, Nthread >>> ( dev_diag ) ;  
    // diagnostic copied to host out here
    printf( "calculating diagnostics \n " );
    int size_diag = NN * PARAM.qdiag * sizeof(double) ;
    cudaMemcpy( DIAG, dev_diag, size_diag, cudaMemcpyDeviceToHost);
    wDIAG( DIAG, (*TT).time, PARAM );
  }
// take the first substep
  rnkf45_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_kin, dev_psi, dev_tt ) ;
  rnkf45_time_substep( TT, dev_tt , 0 ) ;
  // 2nd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k2, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug  );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the second substep
  rnkf45_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_kin, dev_psi, dev_tt ) ;
  rnkf45_time_substep( TT, dev_tt , 1 ) ;
  // 3rd evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k3, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // take the third substep
  rnkf45_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_kin, dev_psi, dev_tt ) ;
  rnkf45_time_substep( TT, dev_tt , 2 ) ;
  // 4th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k4, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug);
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}

  // take the fourth substep
  rnkf45_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_kin, dev_psi, dev_tt ) ;
  rnkf45_time_substep( TT, dev_tt , 3 ) ;
  // 5th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k5, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}

  // take the fifth substep
  rnkf45_psi_substep<<<Nblock,Nthread>>>( dev_psip,  dev_kin, dev_psi, dev_tt ) ;
  rnkf45_time_substep( TT, dev_tt , 4 ) ;
  // 5th evaluation of rhs, no diagnostic calculated
  eval_rhs<<<Nblock,Nthread >>>( dev_k6, dev_psip,  dev_tt ,
                                 dev_param, dev_diag, dev_bug );
  // check if there were any bugs from rhs evaluation
  cudaMemcpy( &BUG, dev_bug, size_CRASH, cudaMemcpyDeviceToHost);
  if ( BUG.lstop) { IStop( BUG );}
  // final step
  rnkf45_psi_step<<<Nblock, Nthread >>>( dev_psi, dev_psip, dev_err, dev_kin, dev_tt);
  // psip will store the 4th order rnkt4 solution. 
  // rnkf45_psi_step<<<Nblock, Nthread >>>( dev_psip, dev_kin, dev_tt, 1 );
  // rnkf45_error<<Nblock,Nthread>>> (dev_err,dev_psi,dev_psip);
  double maxErr = sqrt(MaxDevArray( dev_err, Nblock, Nthread));
  // printf("%lf\n",maxErr);
  bool laccept = rnkf45_time_step(TT,dev_tt,maxErr);
  if (!laccept){
    rnkf45( PSI, dev_psi,
          TT, dev_tt,
          PARAM, dev_param ,
          DIAG, dev_diag,
          BUG, dev_bug,
          Nblock, Nthread );
  }
}
