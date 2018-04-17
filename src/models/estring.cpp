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
    dR[ip]=EForce*OneByGamma;
  }
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
vec3 dHdR(int kp, vec3 X[]){
  vec3 ukm2, ukm1, uk, ukp1, Xzero;
  double bkm2, bkm1, bk, bkp1;
  vec3 FF;
  Xzero.x=0.; Xzero.y=0.;Xzero.z=Z0;
  switch (kp){
  case 0:
  case 1:
    vec3 dX = X[kp-2+1]-Xzero;
    bkm2 = norm(dX);
    ukm2=dX/bkm2;
    getub(&bkm1, &ukm1, kp-1, X);
    getub(&bk, &uk, kp, X);
    getub(&bkp1, &ukp1, kp+1, X);
    FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
	       + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
	       - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
	       );
  case Np-2:
  case Np-1:
  default:
    getub(&bkm2, &ukm2, kp-2, X);
    getub(&bkm1, &ukm1, kp-1, X);
    getub(&bk, &uk, kp, X);
    getub(&bkp1, &ukp1, kp+1, X);
    FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
		    + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
		    - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
		    );
  }
  return FF*(AA/aa);
}
/**************************/
void iniconf(double y[]){
  vec3 R[Np];	// R is the position of the beads.
  for (int ip=0;ip<Np;ip++){
    R[ip].x=0.;
    R[ip].y=0.;
    R[ip].z=aa*(double)ip;
    y[3*ip]=R[ip].x;
    y[3*ip+1]=R[ip].y;
    y[3*ip+2]=R[ip].z;
  }
}
