#include <iostream>
#include<fstream>
#include "ode.h"
#include "modules/3vec.h"
#include "model.h"
#include <cmath>
using namespace std;

/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, double* add_SS, bool flag_kappa);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
  vec3 R[Np],dR[Np], EForce;	// R is the position of the beads.
  // double CurvSqr[Np];
  double kappasqr, MaterialCoordinate;
  for (int ip=0;ip<Np;ip++){
    R[ip].x=y[3*ip];
    R[ip].y=y[3*ip+1];
    R[ip].z=y[3*ip+2];
  }

  for (int ip=0;ip<Np;ip++){
    kappasqr=CurvSqr[ip];
    MaterialCoordinate=SS[ip];
    dHdR(ip, R, &EForce, &kappasqr, &MaterialCoordinate, flag_kappa);
    dR[ip]=EForce*1/Gamma;
    CurvSqr[ip]=kappasqr;
    SS[ip]=MaterialCoordinate;
    // cout << CurvSqr[ip] << endl;
  }
  
  // External force applied on the end point.
  vec3 FF0(0., 0., -FFZ0*sin(omega*time));
  // cout << FF0.z <<endl;
  dR[Np-1] = dR[Np-1]-FF0*1/Gamma; 
  //dR[Np-1].y = 0;                     // Constraint that last point should always remain on z axis. 
  
  for (int ip=0;ip<Np;ip++){
    rhs[3*ip]=dR[ip].x;
    rhs[3*ip+1]=dR[ip].y;
    rhs[3*ip+2]=dR[ip].z;
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
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, double* add_SS, bool flag_kappa){
	// This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 FF = *add_FF;              // Since I am passing the address of force in add_FF and the same goes for Kapppsqr

  //vec3 FF;
  Xzero.x=0.; Xzero.y=0.;Xzero.z=Z0;

  if (kp == 0)
  {
    dX = X[kp-1+1]-Xzero;
    bkm1 = norm(dX);
    ukm1=dX/bkm1;
    getub(&bk, &uk, kp, X);
    getub(&bkp1, &ukp1, kp+1, X);
    FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
         + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
         - (ukm1/bkm1)*( dot(ukm1,uk) )
         );
    FF = FF*AA/aa;
    // Add an extra term for inextensibility constraint
    FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa; 
    // cout << FF.z << endl;
    *add_kappasqr=0.;
    *add_FF = FF;
  }

  else if (kp == 1)
  {
    dX = X[kp-2+1]-Xzero;
    bkm2 = norm(dX);
    ukm2 = dX/bkm2;
    getub(&bkm1, &ukm1, kp-1, X);
    getub(&bk, &uk, kp, X);
    getub(&bkp1, &ukp1, kp+1, X);
    FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
        + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
        - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
        );
    FF = FF*(AA/aa);
    // cout << FF.z << endl;
    FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa) )*HH/aa;   // Inextensibility constraint
    *add_kappasqr=0.;
    *add_FF = FF;
  }


  else if (kp == Np-2)
  {
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
  }

  else if (kp == Np-1)
  {
    getub(&bkm2, &ukm2, kp-2, X);
    getub(&bkm1, &ukm1, kp-1, X);
  
    FF = (     (ukm2)/bkm1
        - (ukm1/bkm1)*( dot(ukm1,ukm2) )
        );
    FF = FF*(AA/aa);
    // cout << bkm1 << endl;
    FF = FF - (ukm1*(bkm1-aa))*HH/aa;
    *add_kappasqr=0.;
    *add_FF = FF;
  }

  else
  {
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

    if (flag_kappa)
    {
      *add_kappasqr=2.*(1.-dot(uk,ukm1));
      // cout << kappasqr << endl;
    }
    *add_FF = FF;
  }
  *add_SS = (kp+1)*bkm1;
}
/**************************/
void diagnos(int p){
  // cout << time<<'\t' ;
  // for (int ip=0;ip<Np;ip++){
  //   cout << y[3*ip] <<'\t' << y[3*ip+1]  <<'\t' << y[3*ip+2] << '\t' ;
  // }
  // cout << '\n';

  // ofstream outfile;
  // string l = "output/position";
  // l.append(to_string(itn));
  // l.append(".txt");
  // outfile.open(l, ios::out);

  // for (int ip = 0; ip < Np; ++ip)
  //   {
  //     outfile << y[3*ip] << '\t' << y[3*ip+1] << '\t' << y[3*ip+2] << endl ;
  //   }

  // outfile.close(); 
  cout << p << endl;
  // cout << time << '\t' << y[0] << endl;
}
/**************************/
void iniconf(double y[]){
  vec3 R[Np];	// R is the position of the beads.
  double k = 0.5;      // determines the frequency for initial configuration
  for (int ip=0;ip<Np;ip++){
    R[ip].x=0.;
    R[ip].z=aa*double(ip+1);
    R[ip].y=aa*sin(2*M_PI*k*aa*double(ip+1));
    // cout << M_PI << endl;
    y[3*ip]=R[ip].x;
    y[3*ip+1]=R[ip].y;
    y[3*ip+2]=R[ip].z;
  }
}
