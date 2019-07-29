#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cuda.h"
#include "chain.h"
#include "model.h"
#include "3vec.h"
#include "2Tens.h"
using namespace std;
__host__ void write_param( void );
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
__device__ int square_wave( double t, double Tby2) ;
__device__ vec3 Uflow ( vec3 RR, struct MPARAM *param);
__device__ void device_exception( struct CRASH *bug, char mesg[] );
struct MPARAM host_param ;
struct MPARAM *dev_param;
double *DIAG ;
double *dev_diag ;
struct CRASH BUG;
struct CRASH *dev_bug;
int size_diag;
/* ========================================= */
__host__ void set_param( void ){
  host_param.height = 1.;
  double height = host_param.height ;
  host_param.aa = height/(double)(NN-1);
  // distance between two nodes.
  host_param.Dbyell = 0.005 ;
  host_param.dd = height*host_param.Dbyell ;
  /* r/l ratio for the rod has been kept constant. It should be noted that 
     the particles would also have same diameter. */
  host_param.viscosity = 10;	      // Equivalent to kinematic viscosity of glycerin
  host_param.Z0=0. ;		      // If we want the bottom point of the rod to be fixed.
  host_param.Famp = 0. ;	      // Different force for different configuration.
  // Sigma is a dimensionless number, which is described as frequency parameter.
  host_param.sigma=1.5;					
  host_param.ShearRate = 1.;
  host_param.omega = host_param.ShearRate*host_param.sigma ;
  //
  host_param.factorAA = 0.15 ; 
  host_param.AA = host_param.factorAA*pow(10,-4) ; // AA is the bending rigidity.
  host_param.KK = 64.;
  double asqr = host_param.aa*host_param.aa ;
  host_param.HH = host_param.KK*host_param.AA/( asqr );
 // Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
// double TMAX = ShearRate*10;
// double tdiag = TMAX/2000;
  host_param.qdiag = 2 ;
  int qdiag = host_param.qdiag ;
  host_param.bcb = 1 ;
  host_param.bct = 1;
  host_param.global_drag = 1;
  host_param.iext_force = 0;
  host_param.floc = 0 ;
  host_param.iext_flow = 1;
  cudaMalloc( (void**)&dev_param, size_MPARAM  );
  cudaMemcpy( dev_param, &host_param,
                     size_MPARAM, cudaMemcpyHostToDevice ) ;
  write_param( );
  // allocate host data for diagnostic
  size_diag = NN*qdiag*sizeof(double) ; 
  DIAG = (double *) malloc( size_diag ) ;
  cudaMalloc( (void**)&dev_diag, size_diag );
  // allocate space for crashing gracefully.
 cudaMalloc( (void**)&dev_bug, size_CRASH );
  BUG.lstop = 0;
  strcpy(BUG.message, " No bug yet" );
  cudaMemcpy( dev_bug, &BUG, size_CRASH, cudaMemcpyHostToDevice);
}
/* ===================================== */
__host__ void write_param( void ){
  FILE *pout ;
  pout = fopen ( "wparam.txt", "w" );
  printf( "# =========== Model Parameters ==========\n" );
  printf( " #Model : Elastic String \n " ) ;
  printf( "#dimension of ODE:\n pp =  %d \n", pp ) ;
  printf( "#Number of copies:\n  NN = %d\n", NN ) ;
  printf( " height = %f \n" , host_param.height) ;
  printf( " #============================\n" );
  fprintf( pout, "# =========== Model Parameters ==========\n" );
  fprintf( pout, " #Model : Elastic String \n " ) ;
  fprintf( pout, "#dimension of ODE:\n pp =  %d \n", pp ) ;
  fprintf( pout, "#Number of copies:\n  NN = %d\n", NN ) ;
  fprintf( pout, " height = %f \n" , host_param.height) ;
  fprintf( pout, " aa= %f \n" , host_param.aa ) ;
  fprintf( pout, " Dbyell= %f \n" , host_param.Dbyell ) ;
  fprintf( pout, " dd= %f \n",  host_param.dd ) ;
  fprintf( pout, " viscosity = %f \n ",  host_param.viscosity) ;
  fprintf( pout, " Z0 = %f \n ", host_param.Z0) ;
  fprintf( pout, " Famp = %f \n ", host_param.Famp ) ;
  fprintf( pout, " sigma = %f \n ", host_param.sigma ) ;
  fprintf( pout, " ShearRate = %f \n ", host_param.ShearRate ) ;
  fprintf( pout, " omega = %f \n ", host_param.omega ) ;
  fprintf( pout, " factorAA = %f \n ", host_param.factorAA ) ;
  fprintf( pout, " AA = %f \n ", host_param.AA ) ;
  fprintf( pout, " KK = %f \n ", host_param.KK ) ;
  fprintf( pout, " HH = %f \n ", host_param.HH ) ;
  fprintf( pout, " qdiag=%d\n", host_param.qdiag );
  fprintf( pout, " bcb=%d\n", host_param.bcb );
  fprintf( pout, " bct=%d\n", host_param.bct );
  fprintf( pout, " global_drag=%d\n", host_param.global_drag );
  fprintf( pout, " iext_force=%d\n", host_param.iext_force );
  fprintf( pout, " floc=%d\n", host_param.floc );
  fprintf( pout, " iext_flow=%d\n", host_param.iext_flow );
  fprintf( pout, " #============================\n" );
  fclose( pout );
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
    UU.x = (height - R.z)*ShearRate*square_wave( tau, M_PI/omega ) ;
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
  double onebythree = 1./3.;
  double mu0 = onebythree/(M_PI*viscosity*dd);
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
        double c1 = 1/(8*M_PI*viscosity*d_rij);
        double dsqr1 = 1./(d_rij*d_rij);
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
  vec3 Rj (psi[pp*j], psi[pp*j+1], psi[pp*j+2] );
  vec3 Ri (psi[pp*i], psi[pp*i+1], psi[pp*i+2] );
  /*This calculate the distance at two index i and j on the elastic string*/
  *rij = Rj - Ri;
  double Dis = norm( Rj-Ri ); 
  *Distance = Dis;
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
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          struct MPARAM *param, double *diag, CRASH *bug ){
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
  R.x= psi[pp*kelement];
  R.y= psi[pp*kelement + 1];
  R.z= psi[pp*kelement + 2];
  if ( kelement == (NN-1) ){
    ds = 0.;
  } else {
    Rp1.x= psi[(pp+1)*kelement];
    Rp1.y= psi[(pp+1)*kelement + 1];
    Rp1.z= psi[(pp+1)*kelement + 2];
    ds = norm( Rp1-R);
  }
  dHdR( kelement, psi, &EForce, &kappasqr, param, bug );
  /* write diagnostic to corresponding array */
  diag[kelement*qdiag ] = ds ;
  diag[kelement*qdiag +1]  = kappasqr;
  /* add external force to the filament */
  if ( (iext_force) && (kelement == floc) ){
  EForce = EForce -  ext_force( kelement, R, tau, param ) ;}
  /* calculate the viscous (possibly non-local ) drag */
  dR = drag(kelement, psi,  &EForce, param);
  /* contribution from external flow */
  if ( iext_flow  ){ 
    dR = dR + ext_flow( kelement, R, tau, param  ) ; }
  /*------ put the rhs back to the dpsi array ----- */
  dpsi[pp*kelement]       = dR.x  ;
  dpsi[pp*kelement + 1] = dR.y ;
  dpsi[pp*kelement + 2] = dR.z ;
}
/*--------------------------------------------------------------*/
__device__  void getub(double *bk, vec3 *uk, int kp, double psi[]){
  vec3 X(psi[pp*kp], psi[pp*kp+1], psi[pp*kp+2]  ) ;
  vec3 Xp1(psi[(pp+1)*kp], psi[(pp+1)*kp+1], psi[(pp+1)*kp+2]  ) ;
  vec3 dX = Xp1-X;
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
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
