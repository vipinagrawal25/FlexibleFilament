#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "quaternion.h"
#include "matrix.h"
#include "RigidBody.h"
#include "rigid.h"
using namespace std;
/**************************/
unsigned int const ndim=pdim*Nensemble;
int yini=1;
void iniconf(double y[]);
void Diagnostics(double y[], double time);
vec3 inipos(int posini);
vec3 inivel(int velini);
quaternion iniq(int qini);
vec3 iniell(int ellini, quaternion qq_zero);
vec3 fluid_velocity(vec3 );
vec3 DragForce(vec3 UU,RigidBody *rb);
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double time=0.;
  double PI=3.14159;
  iniconf(y);
//----------------------------
  while(time <= TMAX){
    Diagnostics(&y[0],time);
    //    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,y[6],y[7],y[8],y[9],y[10],y[11],y[12]);
    //rnkt4(pdim,&y[0],time,dt);
    euler(pdim,&y[0],time,dt);
    time=time+dt;
  } 
  cout<<"Done, time="<<time-dt<<"\t TMAX="<<TMAX<<"\n";
}
/* ----------------------------------------*/
void eval_rhs(double time,double y[],double rhs[]){
  RigidBody rb;
  vec3 UU,vv,Force,Torque;
  quaternion omegaq,dq_dt;
  //
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2rb(&rb, &y[irb]);
    UU=fluid_velocity(rb.xx);
//----- Calculate velocity -------
    vv = rb.pp*(1./Mass);
//---Calculate Force---------
    Force=DragForce(UU,&rb);
    //PVec3(Force);
/* ---- Calculate the evolution eqn for the quaternion 
   dq_dt = (1/2) omega * q  */
    omegaq = quaternion(0.,rb.omega);
    dq_dt = omegaq * rb.qq *(1./2.);
// -- Now for zero torque -----------
    Torque=vec3(0.,0.,0.);
// ----------------------------------
    rhs[irb+0]=vv.x;
    rhs[irb+1]=vv.y;
    rhs[irb+2]=vv.z;
    //
    rhs[irb+3]=Force.x;
    rhs[irb+4]=Force.y;
    rhs[irb+5]=Force.z;
    //
    rhs[irb+6]=dq_dt.w;
    rhs[irb+7]=dq_dt.u.x;
    rhs[irb+8]=dq_dt.u.y;
    rhs[irb+9]=dq_dt.u.z;
    //
    rhs[irb+10]=Torque.x;
    rhs[irb+11]=Torque.y;
    rhs[irb+12]=Torque.z;
    }
}
/* ----------------------------------------*/
void iniconf(double y[]){
  RigidBody rb;
  int ibody,k;
  double pi;
  srand(time(NULL));
  for (ibody=0;ibody<Nensemble;ibody++){
    rb.xx=inipos(1);
    vec3 vv=inivel(1);
    rb.pp = vv*Mass;
    /* **The quaternion must be set before the anglar momentum
     is set** */
    rb.qq=iniq(q_ini);
    rb.Ell=iniell(ell_ini,rb.qq);
    rb2array(&y[ibody],&rb);
  }
}
/* ---------------------------------- */
vec3 inipos(int posini){
  vec3 xx;
  double x1,x2,x3,pi;
  switch(posini){
  case 1:
    pi=4.*atan(1.);
    x1=rand();
    x2=rand();
    x3=rand();
    xx=vec3(x1,x2,x3);
    xx = xx*(1./norm(xx) );
    break;
  default:
    exit(EXIT_FAILURE);
  }
  return xx;
}
/* -------------------------------- */
vec3 inivel(int velini){
  vec3 vv;
  double v1,v2,v3;
  switch(velini){
  case 1:
    v1=rand();
    v2=rand();
    v3=rand();
    vv=vec3(v1,v2,v3);
    vv = vv*(1./norm(vv) );
    break;
  default:
    exit(EXIT_FAILURE);
  }
  return vv;
}
/* ---------------------------------- */
quaternion iniq(int qini){
  quaternion qq;
  switch(qini){
  case 1:
    {
    //Initially body-fixed axis and space-fixed axis coincided.
    qq=qzero;
    break; }
  case 2:
    {
    /*Initially body-fixed axis are rotated by an angle psi 
      about the the x axis of the space-fixed axis */
      cout << "Initial quaternion type="<<qini<<"\n";
      cout << "Rotated about x axis with angle psi_zero="<<psi_zero<<"\n";
      vec3 nhat= vec3(1.,0.,0.);
      double psi = psi_zero;
      qq=quaternion( cos(psi/2.), nhat*sin(psi/2.) );
      cout << "Initial quaternion "<<"\n";
      PQ(qq);
      break; }
  default:
    exit(EXIT_FAILURE);
  }
  return qq;
}
/* ---------------------------------- */
vec3 iniell(int ellini, quaternion qq_zero){
  vec3 ell;
  switch(ellini){
  case 1:
    {
      //zero angular momentum    
      cout << "Initial angular momentum="<<ellini<<"\n";
      cout << "ZERO angular momentum\n";
      ell=vec3(0.,0.,0.);
      break; }
  case 2:
    {
      double RR[3][3], RRt[3][3], Ib[3][3],temp[3][3],Ispace[3][3];
      //angular velocity, omega_zero along z axis
      cout << "Initial angular momentum="<<ellini<<"\n";
      cout << "Angular velocity along z axis with omega_zero="<<omega_zero<<"\n";
      vec3 om=vec3 (0.,0.,omega_zero);
      PQ(qq_zero);
      Rot_from_Q(RR, qq_zero);
      MatTrans(RRt,RR); // the transpose of the rotation matrix
      /* The moment-of-inertia in space fixed frame */
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          Ib[j][k]=Ibody[j][k];
        }
      }
      MatMul(temp,Ib,RRt);
      MatMul(Ispace, RR, temp);
      ell=MatVec(Ispace,om);
      break; }
  default:
    exit(EXIT_FAILURE);
  }
  return ell;
}
/* ----------------------------------------*/
vec3 fluid_velocity(vec3 xx){
  vec3 UU;
  double A_abc=1.;
  double B_abc=1.;
  double C_abc=1.;
  switch(Utype){
  case(1):
    UU.x=A_abc*sin(kk*xx.z) + C_abc*cos(kk*xx.y);
    UU.y=B_abc*sin(kk*xx.x) + A_abc*cos(kk*xx.z);
    UU.z=C_abc*sin(kk*xx.y) + B_abc*cos(kk*xx.x);
    break;
  case(2):
    break;
  default:
    exit(EXIT_FAILURE);
  }
  return UU;
}
/* ----------------------------------------*/
vec3 DragForce(vec3 UU,RigidBody *rb){
  vec3 vv,Force;
  vv = rb->pp*(1./Mass);
  Force= (UU-vv)*(1./St0);
  /*cout<<"Force--vv---UU ---\n";
  PVec3(Force);
  PVec3(vv);
  PVec3(UU);
  cout<<"---------------\n"; */
  return Force;
}
/* ----------------------------------------*/
void Diagnostics(double y[], double time){
  RigidBody rb;
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2rb(&rb, &y[irb]);
  //cout << "omega and ell\n";
  //PVec3(rb->omega);
  //PVec3(rb->Ell);
    cout <<time<<"\t"<<dot(rb.omega,rb.Ell)<<"\n";
  //PVec3(rb->omega);
  }
}
/* ----------------------------------------*/
