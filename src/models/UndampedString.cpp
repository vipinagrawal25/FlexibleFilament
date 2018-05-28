/* This model is aimed to solve a rod using set of springs where it is considered that surrounding fluid is not affecting the rod.
Using that assumption the problem reduces to solve Mass*Acceleration = -Gradient of energy. Whereas in Damped string Gamma factor was involved
which gives rise to only first order equation as described in Wada&Netz(DOI: 10.1209/epl/i2006-10155-0) */
#include <iostream>
#include<fstream>
#include "ode.h"
#include "modules/3vec.h"
#include "model.h"
using namespace std;
/**************************/
vec3 dHdR(int kp, vec3 X[]);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
/**************************/
void eval_rhs(double time,double y[],double rhs[]){
  vec3 R[Np],dR[Np];	// R is the position of the beads.
  for (int ip=0;ip<Np;ip++){
    R[ip].x=y[3*ip];
    R[ip].y=y[3*ip+1];
    R[ip].z=y[3*ip+2];   
  }

  for (int ip=0;ip<Np;ip++){
    vec3 EForce=dHdR(ip, R);
    dR[ip]=EForce*OneByMass;
  }
  
  // External force applied on the end point.
  vec3 FF0(0., 0., FFZ0);
  // cout << FF0.z <<endl;
  dR[Np-1] = dR[Np-1]-FF0*OneByMass; // This needs to be changed.
  // cout << dR[Np-1].z << endl;
  
  for (int ip=0;ip<Np;ip++){
    rhs[3*ip]= y[3*Np+3*ip];
    rhs[3*ip+1] = y[3*Np+3*ip+1];
    rhs[3*ip+2] = y[3*Np+3*ip+2];   
  
    rhs[3*ip+3*Np]=dR[ip].x;
    rhs[3*ip+1+3*Np]=dR[ip].y;
    rhs[3*ip+2+3*Np]=dR[ip].z;
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
vec3 dHdR(int kp, vec3 X[]){
	// This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 FF;
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
    // cout << FF.z << endl;     
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
    // cout << FF.z << endl;
    FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
  }

  return FF;

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
  for (int ip=0;ip<Np;ip++){
    R[ip].x=0.;
    R[ip].y=0.;
    R[ip].z=aa*double(ip+1);
    
    y[3*ip]=R[ip].x;
    y[3*ip+1]=R[ip].y;
    y[3*ip+2]=R[ip].z;

    // In this representation, the idea is to store position vector till 3*Np-1 indexes and after that velocity would be stored into the array.
    y[3*Np+3*ip] = 0;
    y[3*Np+3*ip+1] = 0;
    y[3*Np+3*ip+2] = 0;
  }
}
