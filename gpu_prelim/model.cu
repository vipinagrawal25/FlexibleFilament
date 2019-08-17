#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cuda.h"
#include "model.h"
#include "3vec.h"
#include "2Tens.h"
using namespace std;
__device__ void dHdR(int kp, double psi[], vec3* add_FF,
                     double* add_kappasqr,  struct MPARAM *param,
                     struct CRASH * );
__device__ vec3 drag(int ip,  double psi[], vec3 EForce[],
                     struct MPARAM *param);
__device__ vec3 ext_flow( int kelement, vec3 R, double tau,
                          struct MPARAM *param  );
__device__ vec3 ext_force( int kelement, vec3 R, double tau,
                           struct MPARAM *param  );
__device__ void GetRij(double psi[], int i, int j, double *Distance,
                       vec3 *rij);
__device__  void getub(double *bk, vec3 *uk, int kp, double psi[]);
__device__ int square_wave(double t, double Tby2) ;
__device__ vec3 Uflow ( vec3 RR, struct MPARAM *param);
__device__ void device_exception( struct CRASH *bug, char mesg[] );
__device__ vec3 psi2R(double psi[], int k);
__device__ void R2psi(double psi[], int k, vec3 R);
/* ========================================= */
void alloc_chain(  double **PSI, double **dev_psi ){
  double *temp;
  *PSI = (double *) malloc( ndim*sizeof( double ) );
  cudaMalloc( (void**)&temp, ndim*sizeof( double ) );
  *dev_psi = temp;
}
/* ------------------------------------------------------------------------------*/
int pre_diag( double **DIAG , double **dev_diag, MPARAM PARAM ){
 // allocate host data for diagnostic
  int qdiag = PARAM.qdiag ; 
  int size_diag = NN * qdiag * sizeof(double) ; 
  *DIAG = (double *) malloc( size_diag ) ;
  for (int iN=0; iN<NN; iN++ ){
    for(int iq=0; iq<qdiag; iq++){
      (*DIAG)[iN*qdiag+iq] = 1. ;
    }
  }
  double *temp;
  cudaMalloc( (void**)&temp, size_diag );
  *dev_diag = temp ;
  cudaMemcpy( *dev_diag,  *DIAG, size_diag,
              cudaMemcpyHostToDevice);
  return size_diag;
}
/*--------------------------------------------------------*/
void set_param( MPARAM *PARAM, MPARAM **dev_param ){
  MPARAM *temp;
  double height = 1.;
  (*PARAM).height = height;
  (*PARAM).aa = height/(double)(NN-1);
  // distance between two nodes.
  (*PARAM).Dbyell = 0.005 ;
  (*PARAM).dd = height*(*PARAM).Dbyell ;
  /* r/l ratio for the rod has been kept constant. It should be noted that 
     the particles would also have same diameter. */
  (*PARAM).viscosity = 10;	      // Equivalent to kinematic viscosity of glycerin
  (*PARAM).Z0=0. ;		      // If we want the bottom point of the rod to be fixed.
  (*PARAM).Famp = 0. ;	      // Different force for different configuration.
  // Sigma is a dimensionless number, which is described as frequency parameter.
  (*PARAM).sigma=1.5;					
  (*PARAM).ShearRate = 1.;
  (*PARAM).omega = (*PARAM).ShearRate*(*PARAM).sigma ;
  //
  (*PARAM).factorAA = 0.15 ;
  // (*PARAM).factorAA = 0. ;
  (*PARAM).AA = (*PARAM).factorAA*pow(10,-4) ; // AA is the bending rigidity.
  (*PARAM).KK = 64.;
  double asqr = (*PARAM).aa*(*PARAM).aa ;
  (*PARAM).HH = (*PARAM).KK*(*PARAM).AA/( asqr );
// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
// double TMAX = ShearRate*10;
// double tdiag = TMAX/2000;
  (*PARAM).qdiag = 2 ;
  int qdiag = (*PARAM).qdiag ;
  (*PARAM).bcb = 1 ;      // Boundary condition at bottom
  (*PARAM).bct = 1;       // Boundary condition at top
  (*PARAM).global_drag = 1; 
  (*PARAM).iext_force = 0; // Whether to apply the external force or not
  (*PARAM).floc = 0 ;       // External force location 
  (*PARAM).iext_flow = 1;   // External flow: Yes/No
  (*PARAM).iniconf = 1;     // Configuration of the system at t = 0.
  cudaMalloc( (void**)&temp, size_MPARAM  );
  *dev_param = temp;
  cudaMemcpy( *dev_param, PARAM,
                     size_MPARAM, cudaMemcpyHostToDevice ) ;
  write_param( PARAM, "wparam.txt" );
} 
  // allocate space for crashing gracefully.
  /* cudaMalloc( (void**)&dev_bug, size_CRASH );
  BUG.lstop = 0;
  strcpy(BUG.message, " No bug yet" );
  cudaMemcpy( dev_bug, &BUG, size_CRASH, cudaMemcpyHostToDevice); */
/* ------------------------------------------------------------------------------*/
__host__ void write_param( MPARAM *PARAM, char *fname ){
  FILE *pout ;
  pout = fopen ( fname, "w" );
  printf( "# =========== Model Parameters ==========\n" );
  printf( " #Model : Elastic String \n " ) ;
  printf( "#dimension of ODE:\n pp =  %d \n", pp ) ;
  printf( "#Number of copies:\n  NN = %d\n", NN ) ;
  printf( " height = %f \n" , (*PARAM).height) ;
  printf( " #============================\n" );
  fprintf( pout, "# =========== Model Parameters ==========\n" );
  fprintf( pout, " #Model : Elastic String \n " ) ;
  fprintf( pout, "#dimension of ODE:\n pp =  %d \n", pp ) ;
  fprintf( pout, "#Number of copies:\n  NN = %d\n", NN ) ;
  fprintf( pout, " height = %f \n" , (*PARAM).height) ;
  fprintf( pout, " aa= %f \n" , (*PARAM).aa ) ;
  fprintf( pout, " Dbyell= %f \n" , (*PARAM).Dbyell ) ;
  fprintf( pout, " dd= %f \n",  (*PARAM).dd ) ;
  fprintf( pout, " viscosity = %f \n ",  (*PARAM).viscosity) ;
  fprintf( pout, " Z0 = %f \n ", (*PARAM).Z0) ;
  fprintf( pout, " Famp = %f \n ", (*PARAM).Famp ) ;
  fprintf( pout, " sigma = %f \n ", (*PARAM).sigma ) ;
  fprintf( pout, " ShearRate = %f \n ", (*PARAM).ShearRate ) ;
  fprintf( pout, " omega = %f \n ", (*PARAM).omega ) ;
  fprintf( pout, " factorAA = %f \n ", (*PARAM).factorAA ) ;
  fprintf( pout, " AA = %f \n ", (*PARAM).AA ) ;
  fprintf( pout, " KK = %f \n ", (*PARAM).KK ) ;
  fprintf( pout, " HH = %f \n ", (*PARAM).HH ) ;
  fprintf( pout, " qdiag=%d\n", (*PARAM).qdiag );
  fprintf( pout, " bcb=%d\n", (*PARAM).bcb );
  fprintf( pout, " bct=%d\n", (*PARAM).bct );
  fprintf( pout, " global_drag=%d\n", (*PARAM).global_drag );
  fprintf( pout, " iext_force=%d\n", (*PARAM).iext_force );
  fprintf( pout, " floc=%d\n", (*PARAM).floc );
  fprintf( pout, " iext_flow=%d\n", (*PARAM).iext_flow );
  fprintf( pout, " iniconf=%d\n", (*PARAM).iniconf );
  fprintf( pout, " #============================\n" );
  fclose( pout );
}
/* -----------------------------------------------------------------------------------*/
__device__ vec3 psi2R(double psi[], int k){
  vec3 Rt;
  Rt.x = psi[3*k]; 
  Rt.y = psi[3*k+1]; 
  Rt.z = psi[3*k+2];
  return Rt;
}
/* -----------------------------------------------------------------------------------*/
__device__ void R2psi(double psi[], int k, vec3 R){
  psi[3*k] = R.x; 
  psi[3*k+1] = R.y; 
  psi[3*k+2] = R.z;
}
/* -----------------------------------------------------------------------------------*/
__device__ int square_wave( double t, double Tby2) {
  int s = t/Tby2 ;
  int sw = -2*(s % 2 ) + 1;
  return sw;
}
/* -----------------------------------------------------------------------------------*/  
__device__ vec3 Uflow ( vec3 RR, struct MPARAM *param){
  double gdot = (*param).ShearRate ;
  vec3 UU( gdot*RR.y, 0., 0. );
  return UU;
}
/* -----------------------------------------------------------------------------------*/
__device__ vec3 ext_force( int kelement, vec3 R, double tau,
                           struct MPARAM *param  ){
  double omega = (*param).omega;
  double Famp = (*param).Famp;
  vec3 FF0;
  /* iext_force : implements external force on the filament
     periodic forcing at  position floc  */
      FF0.x = 0. ;
      FF0.y = 0. ;
      FF0.z = -Famp*sin(omega*tau) ;
      return FF0;
}

/* -----------------------------------------------------------------------------------*/  
__device__ vec3 ext_flow( int kelement, vec3 R, double tau,
                          struct MPARAM *param  ){
  int iext_flow = (*param).iext_flow;
  double ShearRate = (*param).ShearRate;
  double omega = (*param).omega;
  double height = (*param).height;
  vec3 UU ;
  switch( iext_flow ){
  case 1:
    //time-dependent shear U = ( ShearRate*z, 0, 0 ) * square_wave(omega*time) 
    UU.x = (height - R.z)*ShearRate*(double)square_wave( tau, M_PI/omega ) ;
    UU.y = 0. ;
    UU.z = 0;
    break;
  case 2:
    UU.x = R.z*ShearRate ;
    break;
  }
  return UU;
}
/*--------------------------------------------------------------------------------------*/
__device__ vec3 drag(int ip,  double psi[], vec3 *EForce, struct MPARAM *param){
  Tens2 dab(1.,0.,0.,0.,1.,0.,0.,0.,1.);
  vec3 dR(0., 0., 0.) ;
  double viscosity = (*param).viscosity ;
  double dd = (*param).dd;
  double onebythree = 1./3;
  double mu0 = onebythree/(M_PI*viscosity*dd);
  double c1, dsqr1;
  if ( (*param).global_drag ){
    /* mu_ij represents the one element of mobility matrix (Size: NXN). 
       Every element of the matrix itself is a 2nd rank tensor with dimension 3x3.*/
    Tens2 mu_ij, mu_ii;
    double d_rij;
    vec3 rij;
    //
    mu_ii = dab*mu0;
    /* dab is Kroneker delta in 2d. It is defined in module/2Tens file. */
    // PTens2(mu_ii);
    // rij = R[j]-R[i] and d_rij is just the norm of this value.
    for (int jp = 0; jp < NN; ++jp){
      if (jp == ip){
        dR =  dR + dot(mu_ii, *EForce );
      }else{
        GetRij(psi, ip, jp, &d_rij, &rij);
        c1 = 1./(d_rij*8*M_PI*viscosity);
        dsqr1 = 1./(d_rij*d_rij);
        mu_ij = c1*(dab + (rij*rij)*dsqr1 +
                    dd*dd/(2*d_rij*d_rij)*(dab*onebythree - (rij*rij)*dsqr1));
        dR =  dR + dot(mu_ij, *EForce );
      }
    }
  } else{
    /* if we use local drag */
    dR = (*EForce)*mu0;
  } 
  return dR;
}
/*----------------------------------------------------------------------------------------------*/
__device__ void GetRij(double psi[], int i, int j, double *Distance,
                       vec3 *rij){

  vec3 Rj =  psi2R( psi,  j);
  vec3 Ri =  psi2R( psi,  i);

// __device__ void R2psi(double psi[], int k, vec3 R);

//   vec3 Rj (psi[pp*j], psi[pp*j+1], psi[pp*j+2] );
//   vec3 Ri (psi[pp*i], psi[pp*i+1], psi[pp*i+2] );
  /*This calculate the distance at two index i and j on the elastic string*/
  *rij = Rj - Ri;
  double Dis = norm(Rj-Ri); 
  *Distance = Dis;
}
/*--------------------------------------------------------------*/
__device__  void getub(double *bk, vec3 *uk, int kp, double psi[]){
  vec3 X = psi2R( psi, kp ) ;
  vec3 Xp1 = psi2R ( psi, kp+1) ;
  vec3 dX = Xp1-X;
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ vec3 Force_FirstPoint( double psi[],
                                  struct MPARAM *param, struct CRASH *bug ){
  double AA = (*param).AA ;
  double aa = (*param).aa ;
  double HH = (*param).HH;
  int bcb = (*param).bcb ;    
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.);
  vec3 Xzero(0.,0.,0.), dX(0.,0.,0.);
  double  bkm1, bk, bkp1;
  vec3 Xbot( psi[0], psi[1], psi[2]  );
  vec3 FF;
  getub(&bk, &uk, 0, psi);
  getub(&bkp1, &ukp1, 1, psi);
  switch( bcb ){
  case 0: 
    dX = Xbot - Xzero;
    bkm1 = norm(dX);
    ukm1=dX/bkm1;
        FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
               + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
               - (ukm1/bkm1)*( dot(ukm1,uk) )
               );
          FF = FF*AA/aa;
          // Add an extra term for inextensibility constraint
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa; 
          break;
    case 1: // free
        FF = ( (uk/bk)*( dot(uk,ukp1) )  - (ukp1)/bk );
        FF = FF*AA/aa;
        // Add an extra term for inextensibility constraint
        FF = FF + ( uk*(bk-aa))*HH/aa; 
        break;
  default: // we must crash now.
    device_exception( bug, "NN=0,  bcb not implemented " );
    break;
  }
  return FF;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ vec3 Force_SecondPoint( double psi[],
                                  struct MPARAM *param, struct CRASH *bug ){
  double AA = (*param).AA ;
  double aa = (*param).aa ;
  double HH = (*param).HH;
  int bcb = (*param).bcb ;
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.);
  vec3 Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 Xbot( psi[0], psi[1], psi[2]  );
  vec3 FF;
  getub(&bkm1, &ukm1, 0, psi);
  getub(&bk, &uk, 1, psi);
  getub(&bkp1, &ukp1, 2, psi);
  switch( bcb ){
    case 0: //clamped
      dX = Xbot-Xzero;
      bkm2 = norm(dX);
      ukm2 = dX/bkm2;
      FF = (  (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
              );
      FF = FF*(AA/aa);
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa) )*HH/aa;   // Inextensibility constraint
      break;
    case 1: // free
      FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
                 + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
                 - (ukm1/bkm1)*( dot(ukm1,uk) )
                 );
      FF = FF*(AA/aa);
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;   // Inextensibility constraint
      break;
  default: // any other boundary conditions.
    device_exception( bug, "NN=1, bcb not implemented ");
    break;
  }
  return FF;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ vec3 Force_NNm2( double psi[],
                                  struct MPARAM *param, struct CRASH *bug ){
  double AA = (*param).AA ;
  double aa = (*param).aa ;
  double HH = (*param).HH;
  int bct = (*param).bct ;
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.);
  vec3 Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk;             
  vec3 FF; 
  getub(&bkm2, &ukm2, NN-4, psi);
  getub(&bkm1, &ukm1, NN-3, psi);
  getub(&bk, &uk, NN-2, psi);
  FF = (     (uk+ukm2)/bkm1 - (ukm1)/bk
             + (uk/bk)*( dot(uk,ukm1))
             - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
             );    
  FF = FF*(AA/aa);
  FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint
  return FF;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ vec3 Force_NNm1( double psi[],
                                  struct MPARAM *param, struct CRASH *bug ){
  double AA = (*param).AA ;
  double aa = (*param).aa ;
  double HH = (*param).HH;
  int bct = (*param).bct ;
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.);
  vec3 Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1;             
  vec3 FF; 
  switch( bct ){
  case 0: // clamped
    device_exception( bug , "element NN-1, bct=0 not implemented "); 
  break;
 case 1: //free 
   getub(&bkm2, &ukm2, NN-3, psi);
   getub(&bkm1, &ukm1, NN-2, psi);
   FF = (     (ukm2)/bkm1
              - (ukm1/bkm1)*( dot(ukm1,ukm2) )
              );
   FF = FF*(AA/aa);
   FF = FF - (ukm1*(bkm1-aa))*HH/aa;
   break;
 default: // any other bct .
   device_exception( bug, "element NN-1, bct not implemented ");
   break;
  }
  return FF;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ vec3 Force_rest( double* add_kappasqr, int kp, double psi[],
                                  struct MPARAM *param, struct CRASH *bug ){
  double AA = (*param).AA ;
  double aa = (*param).aa ;
  double HH = (*param).HH;
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.);
  vec3 Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;             
  vec3 FF;
  getub(&bkm2, &ukm2, kp-2, psi);
  getub(&bkm1, &ukm1, kp-1, psi);
  getub(&bk, &uk, kp, psi);
  getub(&bkp1, &ukp1, kp+1, psi);
  FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
             + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
             - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
             );
  FF = FF*(AA/aa);
  FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint
  *add_kappasqr=2.*(1.- dot(uk,ukm1))/(aa*aa);
  return FF;
}
/*-----------------------------------------------------------------------------------------------*/
__device__ void dHdR(int kp, double psi[], vec3* add_FF,
                     double* add_kappasqr,  struct MPARAM *param, struct CRASH *bug ){
     // This function calculates the force at every node which is a function of X, time.
  switch(kp){
  case 0: // 0th element of the chain.
    *add_FF = Force_FirstPoint( psi, param, bug );
    *add_kappasqr=0.;
    break;     
  case 1: //1st element of the chain.
    *add_FF = Force_SecondPoint( psi, param, bug ); 
    *add_kappasqr=0.;
    break;
  case NN-2:
    *add_FF = Force_NNm2( psi, param, bug) ; 
    *add_kappasqr=0.;
    break;
  case NN-1:
    *add_FF = Force_NNm1( psi, param, bug );
    *add_kappasqr=0.;
    break;
  /* for all other points */
  default:
    *add_FF = Force_rest( add_kappasqr, kp, psi, param, bug ) ;
    break;
  }  
}
 /* -----------------------------------------------------------------------------------*/
__device__ void model_rhs( double dpsi[], double psi[], int kelement, double tau,
                           struct MPARAM *param, double diag[], CRASH *bug, int ldiag ){
  int iext_flow = (*param).iext_flow ;
  vec3 R, dR, EForce, FF0, Rp1;  // R is the position of the beads.
  int iext_force = (*param).iext_force ;
  int floc = (*param).floc ;
  /* we are calculating two diagnostic quantities at the moment
ds : d( material coordinate) . 
kappasqr : square of local curvature. 
This number is stored in param.qdiag */
  int qdiag = (*param).qdiag ;
  double ds, kappasqr ; 
  R = psi2R(psi, kelement );
  if ( kelement == (NN-1) ){
    ds = 0.;
  } else {
    Rp1 = psi2R(psi, kelement+1) ;
    ds = norm( Rp1-R );
  }
  dHdR( kelement, psi, &EForce, &kappasqr, param, bug );
  /* write diagnostic to corresponding array */
  if (ldiag ){
  diag[ kelement*qdiag ] = ds ;
  diag[kelement*qdiag +1]  = kappasqr;
  }
  /* add external force to the filament */
  if ( (iext_force) && (kelement == floc) ){
  EForce = EForce -  ext_force( kelement, R, tau, param ) ;}
  /* calculate the viscous (possibly non-local ) drag */
  dR = drag(kelement, psi,  &EForce, param);
  /* contribution from external flow */
  // if ( iext_flow  ){ 
  //   dR = dR+ext_flow( kelement, R, tau, param ) ; 
  // }
  /*------ put the rhs back to the dpsi array ----- */
  R2psi( dpsi, kelement, dR);
}
/*-------------------------------------------------------------------*/
void initial_configuration( double PSI[], MPARAM PARAM ){
  int iniconf = PARAM.iniconf;
  switch(iniconf){
    case -1:
      // Here the thing goes if you want to start the code from middle.
      // In this case, one more parameter is needed though, the last file.
    case 1:
      /* elastic filament is on a straight line perpendicular to the flow
         with no perturbation.*/
      for (int iN=0; iN<NN; iN++){
        PSI[iN*pp] = 0.;
        PSI[iN*pp + 1] = 0.;
        PSI[iN*pp + 2] = (PARAM.aa)*(double) iN ; 
      }     
  }  
}
/*-------------------------------------------------------------------*/
void  free_chain( double **PSI, double **psi  ){
  free( *PSI );
  cudaFree( *psi );
}
/*-------------------------------------------------------------------*/
void H2D(double dev_arr[], double ARR[], int Nsize){
  cudaMemcpy( dev_arr, ARR, Nsize*sizeof(double), cudaMemcpyHostToDevice);
}
/*-------------------------------------------------------------------*/
void D2H(double ARR[], double dev_arr[], int Nsize){
  cudaMemcpy( ARR, dev_arr, Nsize*sizeof(double), cudaMemcpyDeviceToHost);
}
/*-------------------------------------------------------------------*/
void wPSI ( double PSI[], double tau ){
  FILE *fp = fopen( "data/PSI", "a" );
  fprintf( fp, "%f\t", tau ) ;
  for ( int ichain = 0; ichain< ndim; ichain++ ){ 
    fprintf( fp, "%f\t", PSI[ichain] ) ; 
  }
  fprintf( fp, "\n " );
  fprintf(fp, "#--------------------------------------#\n");
  fclose( fp );
}
/*-------------------------------------------------------------------*/
void wDIAG( double DIAG[], double tau, MPARAM PARAM ){
  FILE *fp = fopen( "data/DIAG", "a" );
  for (int idiag = 0; idiag < PARAM.qdiag ; idiag++){
    fprintf( fp, "%f\t", tau );
    for ( int iN = 0; iN< NN; iN++ ){ 
      fprintf( fp, "%f\t", DIAG[iN + NN*idiag] ) ; 
    }
    fprintf( fp, "\n " );
  }
  fclose(fp);
}
