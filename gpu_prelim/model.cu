#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "odeN.h"
#include "model.h"
__host__ void write_param( void );
struct MPARAM host_param ;
struct MPARAM *dev_param;
double *DIAG ;
double *dev_diag ;
/* ========================================= */
void set_param( void ){
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
  host_param.FFZ0 = 0. ;	      // Different force for different configuration.
  // Sigma is a dimensionless number, which is described as frequency parameter.
  host_param.sigma=1.5;					
  host_param.ShearRate = 1.;
  host_param.omega = ShearRate*sigma;
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
  qdiag = host_param.qdiag ;
  cudaMemcpy( dev_param, &host_param,
                     size_MPARAM, cudaMemcpyHostToDevice ) ;
  write_param( );
  // allocate host data for diagnostic
  int size_diag = NN*qdiag ; 
  *DIAG = (double *) malloc( sizeof(double)*size_diag ) ;
  cudaMalloc( (void**)&dev_diag, size_diag*sizeof( double ) );
}
/* ===================================== */
void write_param( void ){
  printf( "# =========== Model Parameters ==========\n" );
  printf( " #Model : Elastic String \n " ) ;
  printf( "#dimension of ODE:\n pp =  %d \n", pp ) ;
  printf( "#Number of copies:\n  NN = %d\n", NN ) ;
  printf( " height = %f \n" , host_param.height) ;	
  // complete writing out the rest. 
  /*   double aa; 	// distance between two nodes.
  double Dbyell // diameter/length of the filament.
  double dd ;	/* r/l ratio for the rod has been kept constant. 
                   It should be noted that the particles would also have same diameter. */
  /*double viscosity ;				
  double  Z0;	  // If we want the bottom point of the rod to be fixed.
  double FFZ0 ; // Force Value on the ends
// Sigma is a dimensionless number, which is described as frequency parameter.
  double sigma ;					
  double ShearRate ;
  double omega ;
  double  factorAA ; 
  double AA ;
  double HH ;		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
  double KK; */
  printf( " #============================\n" );
}
/* ===================================== */
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          struct MPARAM *param, double *diag ){

  vec3 R, dR, EForce, FF0;  // R is the position of the beads.
  /* we are calculating two diagnostic quantities at the moment
ds : d( material coordinate) . 
kappasqr : square of local curvature. 
This number is stored in param.qdiag */
  double qdiag = param.qdiag ;
  double ds, kappasqr ; 
  double onebythree = 1./3.;
  double Curvlength = 0; // length of the filament
  R.x= psi[pp*kelement];
  R.y= psi[pp*kelement + 1];
  R.z= psi[pp*kelement + 2];
  Rp1.x= psi[(pp+1)*kelement];
  Rp1.y= psi[(pp+1)*kelement + 1];
  Rp1.z= psi[(pp+1)*kelement + 2];
  ds = norm( Rp1-R); // this is not correct yet
  dHdR( kelement, R, &EForce, &kappasqr );

  /* write diagnostic to corresponding array */
  diag[kelement*qdiag +1] = ds ;
  diag[kelement*qdiag +1]  = kappasqr;
  drag(R, dR, EForce);

  switch(conf_number){
    case 0:
      FF0.x = 0;
      FF0.y = 0;
      FF0.z = -FFZ0*sin(omega*time);
      EForce[Np-1] = EForce[Np-1]-FF0;
      break;

    case 1:
      for (int ip = 0; ip < Np; ++ip)
      {
        if (sin(omega*time) >= 0){
          dR[ip].y = dR[ip].y + ShearRate*(height - R[ip].z)*ceil(sin(omega*time));    
        }
        else{
          dR[ip].y = dR[ip].y + ShearRate*(height - R[ip].z)*floor(sin(omega*time));
        }
      }
      break;

      case 2:
      for (int ip = 0; ip < Np; ++ip)
      {
        dR[ip].y = dR[ip].y + ShearRate*(R[ip].z);          
      }
      break; 
  }
  
  // External force applied on the end point.
  // cout << FF0.z <<endl;
  // dR[Np-1] = dR[Np-1]-FF0*; 
  //dR[Np-1].y = 0;                     // Constraint that last point should always remain on z axis. 
  
  dpsi[pp*kelement] = dR.x  ;
  dpsi[pp*kelement + 1] = dR.y ;
  dpsi[pp*kelement + 2] = dR.z ;
/**************************--------------------------------------------------------------*/
void drag(vec3 X[], vec3 dX[], vec3 EForce[]){
  double onebythree = 1./3.;
  double mu0 = onebythree/(M_PI*viscosity*dd);
  if (UseRP == 'Y')
  {
    // mu_ij represents the one element of mobility matrix (Size: NXN). 
    // Every element of the matrix itself is a 2nd rank tensor and the dimension of that should 3x3.

    Tens2 mu_ij, mu_ii;
    double d_rij;
    vec3 rij;

    mu_ii = dab*mu0;    // dab is the unit 2 dimensional tensor. It is defined in module/2Tens file.
    // PTens2(mu_ii);
    // rij = R[j]-R[i] and d_rij is just the norm of this value.

    for (int ip = 0; ip < Np; ++ip)
    {
        // The mu_ij in the next line represents the mobility tensor when j is equal to i and in 
        // response to that the "for loop" is started from ip+1 .

        dX[ip] = dX[ip] + dot(mu_ii, EForce[ip]);
        
        for (int jp = ip+1; jp < Np; ++jp)
        {
            GetRij(X, ip, jp, &d_rij, &rij);
            double c1 = 1/(8*M_PI*viscosity*d_rij);
            double dsqr1 = 1./(d_rij*d_rij);
            mu_ij = c1*(dab + (rij*rij)*dsqr1 + dd*dd/(2*d_rij*d_rij)*(dab*onebythree - (rij*rij)*dsqr1));
            dX[ip] = dX[ip] + dot(mu_ij, EForce[ip]);
            dX[jp] = dX[jp] + dot(mu_ij, EForce[ip]);
        }
    }
  }

  else
  {
    for (int ip = 0; ip < Np; ++ip)
    {
        dX[ip] = EForce[ip]*mu0;
    }
  } 
 
}

/**************************/
void getub(double *bk, vec3 *uk, int kp, vec3 X[]){
  vec3 dX = X[kp+1]-X[kp];
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
}
/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr ){
     // This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 FF = *add_FF;             
 
  
 

  switch(kp){
    case 0:
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (conf_number==0 ){
          dX = X[kp-1+1]-Xzero;
          bkm1 = norm(dX);
          // cout << bkm1 << endl;
          ukm1=dX/bkm1;
          FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
               + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
               - (ukm1/bkm1)*( dot(ukm1,uk) )
               );
          FF = FF*AA/aa;
          // Add an extra term for inextensibility constraint
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa; 
          // cout << FF.z << endl;
          *add_FF = FF;
          // *add_SS = (kp+1)*bkm1;
      }
      else{
          FF = ( (uk/bk)*( dot(uk,ukp1) )  - (ukp1)/bk );
          FF = FF*AA/aa;
          // Add an extra term for inextensibility constraint
          FF = FF + ( uk*(bk-aa))*HH/aa; 
          // cout << FF.z << endl;
          *add_FF = FF;
          // cout << "Kya ye yaha aa raha hai?" << endl;
          // *add_SS = 0;
          // cout << bk << '\t' << aa << endl;  

          break;
      }
      
      *add_kappasqr=0.;
      break;     

    case 1:

      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      
      if (conf_number==0)
      {
          dX = X[kp-2+1]-Xzero;
          bkm2 = norm(dX);
          ukm2 = dX/bkm2;
          FF = (  (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
              );
          FF = FF*(AA/aa);
          // cout << FF.z << endl;
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa) )*HH/aa;   // Inextensibility constraint
          *add_FF = FF;
      }
      else
      {
          FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,uk) )
              );
          FF = FF*(AA/aa);
          // cout << FF.z << endl;
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;   // Inextensibility constraint
          // cout << "Good for you " << endl;
          *add_FF = FF;

      }

      *add_kappasqr=0.;
      // *add_SS = (kp+1)*bkm1;      
      break;

    case Np-2:
      getub(&bkm2, &ukm2, kp-2, X);
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);

      FF = (     (uk+ukm2)/bkm1 - (ukm1)/bk
          + (uk/bk)*( dot(uk,ukm1))
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );    
      FF = FF*(AA/aa);
      // cout << bk << endl;
      // cout << FF.z << endl;
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
      // cout << FF.z << endl;
      *add_kappasqr=0.;
      *add_FF = FF;  
      // *add_SS = (kp+1)*bkm1;
      break;

    case Np-1:
      getub(&bkm2, &ukm2, kp-2, X);
      getub(&bkm1, &ukm1, kp-1, X);
  
      FF = (     (ukm2)/bkm1
        - (ukm1/bkm1)*( dot(ukm1,ukm2) )
        );
      FF = FF*(AA/aa);
      // cout << bkm1 << endl;
      // cout << FF.y << endl;
      FF = FF - (ukm1*(bkm1-aa))*HH/aa;
      // cout << FF.y << endl;
      *add_kappasqr=0.;
      *add_FF = FF;
      // *add_SS = (kp+1)*bkm1;
      break;

    default:
      getub(&bkm2, &ukm2, kp-2, X);
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
      FF = FF*(AA/aa);
      // cout << bkm1 <<endl;
      // cout << FF.y << endl;
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
      // cout << FF.x << endl;
      // cout << FF.y << endl;     

      if (flag_kappa==false)
      {
        *add_kappasqr=2.*(1.- dot(uk,ukm1))/(aa*aa);
        // cout << kappasqr << endl;
        // cout << "This is also high level shit" << endl;
      }
      *add_FF = FF;
      // *add_SS = (kp+1)*bkm1;
      break;
  }  

  //   if (flag_kappa)
  //   {
  //     *add_kappasqr=2.*(1.-dot(uk,ukm1));
  //     // cout << kappasqr << endl;
  //   }
  //   *add_FF = FF;
  // }
  // *add_SS = (kp+1)*bkm1;
 }
/**************************/

















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
