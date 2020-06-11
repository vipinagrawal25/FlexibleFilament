#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "quat.h"
#include "model.h"
using namespace std;
/**************************/
unsigned int const ndim=pdim*Nensemble;
int yini=1;
void iniconf(double y[]);
void Diagnostics(int itn, double y[], double time);
vec3 inipos(int posini);
vec3 inivel(int velini);
quaternion iniq(int qini);
vec3 iniell(int ellini, quaternion qq_zero);
vec3 fluid_velocity(vec3 );
void U_Omega_GradU(vec3& UU, vec3& OO, double GradU[3][3], vec3 xx);
vec3 DragForce(vec3 UU,RigidBody *rb);
vec3 JefferyTorque(vec3 OO,double Sij[3][3],RigidBody *rb);
void wparam_sim(double time);
/* Diagnostic variables */
double yzero[ndim];
double rsqr;
double tauf;
double St_trans;
double nu;
double Re;
double St_rot1;
double St_rot2;
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double time=0.;
  double PI=3.14159;
  int ldiagnos=0;
  wparam_rigid_body();
  iniconf(y);
//----------------------------
  int itn=0;
  while(time <= TMAX){
    //    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,y[6],y[7],y[8],y[9],y[10],y[11],y[12]);
    if (ldiagnos==0) {Diagnostics(itn,&y[0],time);};
    rnkt4(pdim,&y[0],time,dt);
    //euler(pdim,&y[0],time,dt);
    time=time+dt;
    itn=itn+1;
    ldiagnos=itn%idiag;
  } 
  cout<<"Done, time="<<time-dt<<"\t TMAX="<<TMAX<<"\n";
}
/* ----------------------------------------*/
void eval_rhs(double time,double y[],double rhs[]){
  RigidBody rb;
  double GradU[3][3],Sij[3][3];
  vec3 UU,OO,vv,Force,Torque;
  quaternion omegaq,dq_dt;
  //
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2rb(&rb, &y[irb]);
    U_Omega_GradU(UU, OO, GradU, rb.xx);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Sij[i][j]= 0.5*(GradU[i][j]+GradU[j][i]);
      }
    }
    /* vec3 UXO = cross(UU,OO);
    cout<<"U, O, UXO"<<"\n";g
    PVec3(UU);
    PVec3g(OO);
    PVec3(UXO);
    cout<<"-----------------------"<<"\n"; */
//----- Calculate velocity -------
    vv = rb.pp*(1./Mass);
//---Calculate Force---------
    Force=DragForce(UU,&rb);
    //PVec3(Force);
// -- Now for zero torque -----------
    Torque=JefferyTorque(OO,Sij,&rb);
/* ---- Calculate the evolution eqn for the quaternion 
   dq_dt = (1/2) omega * q  */
    omegaq = quaternion(0.,rb.omega);
    dq_dt = omegaq * rb.qq *(1./2.);

// ----------------------------------
    rhs[irb+0]=rb.pp.x/Mass;
    rhs[irb+1]=rb.pp.y/Mass;
    rhs[irb+2]=rb.pp.z/Mass;
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
void U_Omega_GradU(vec3& UU, vec3& OO, double GradU[3][3], vec3 xx){
  double A_abc=1.;
  double B_abc=1.;
  double C_abc=1.;
  switch(Utype){
  case(1):
    {
      A_abc=uzero*A_abc;
      B_abc=uzero*B_abc;
      C_abc=uzero*C_abc;
      UU.x= A_abc*sin(kk*xx.z) + C_abc*cos(kk*xx.y);
      UU.y= B_abc*sin(kk*xx.x) + A_abc*cos(kk*xx.z);
      UU.z= C_abc*sin(kk*xx.y) +  B_abc*cos(kk*xx.x);
      GradU[0][0]=0;
      GradU[0][1]= -kk*C_abc*sin(kk*xx.y);
      GradU[0][2]= kk*A_abc*cos(kk*xx.z);
      //
      GradU[1][0]= kk*B_abc*cos(kk*xx.x);
      GradU[1][1]= 0.;
      GradU[1][2]= -kk*A_abc*sin(kk*xx.z);
      //
      GradU[2][0]= -kk*B_abc*sin(kk*xx.x);
      GradU[2][1]= kk*C_abc*cos(kk*xx.y);
      GradU[2][2]= 0.;
      OO.x=GradU[1][2]-GradU[2][1];
      OO.y=GradU[2][0]-GradU[0][2];
      OO.z=GradU[0][1]-GradU[1][0];
      break;}
  case(2):
    break;
  default:
    exit(EXIT_FAILURE);
  }
};
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
  double KKspace[3][3], Kbody[3][3];
  vv = rb->pp*(1./Mass);
  // The rotation matrix that we get is transpose of the 
  // matrix A used by Mortensen et al.
  for(int j=0;j<3;j++){
    for(int k=0;k<3;k++){
      Kbody[j][k]=KK[j][k];
    }
  } 
  RMRt(KKspace,Kbody,rb->RR);
  Force= MatVec(KKspace, (UU-vv)*(1./taup0) );
  /*cout<<"Force--vv---UU ---\n";
  PVec3(Force);
  PVec3(vv);
  PVec3(UU);
  cout<<"---------------\n"; */
  return Force;
}
/* ----------------------------------------*/
vec3 JefferyTorque(vec3 OO,double Sij[3][3],RigidBody *rb){
  vec3 OOp, Torquep,Torque;
  double denomx,denomy,denomz;
  double RRt[3][3], Sijp[3][3];
  MatTrans(RRt,rb->RR);
  OOp=MatVec(RRt,OO);
  RtMR(Sijp,Sij,rb->RR);
  denomx=1.;
  denomy=1.;
  denomz=1.;
  Torquep.x=Jeff_Rtwo*( (0.5)*(OOp.x-rb->omega.x)  - Jeff_D*Sijp[1][2] );
  Torquep.y=Jeff_Rtwo*( (0.5)*(OOp.y-rb->omega.y) - Jeff_D*Sijp[0][2] );
  Torquep.z=Jeff_Rone*( (0.5)*OOp.z-rb->omega.z  );
  // Transform to torque back to space fixed coordinate 
  Torque=MatVec(rb->RR,Torquep) ;
  return Torque;
}
/* ----------------------------------------*/
void Diagnostics(int itn, double y[], double time){
  RigidBody rb, rb0;
  vec3 UU,OO,vv;
  double GradU[3][3];
  double urms=0;
  if (itn==0){
    for(int iy=0;iy<ndim;iy++){
      yzero[iy]=y[iy];
    }
  };
  rsqr=0.;
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2rb(&rb, &y[irb]);
    array2rb(&rb0, &yzero[irb]);
    U_Omega_GradU(UU, OO, GradU, rb.xx);
    urms=urms+norm(UU);
    rsqr=rsqr+norm(rb.xx-rb0.xx);
  }
  tauf=kk/urms;
  St_trans=taup0/tauf;
  nu=mu*rho;
  Re=urms/(nu*kk);
  double taup_rot1=1./(Jeff_Rone/I0);
  double taup_rot2=1./(Jeff_Rtwo/I0);
  St_rot1=taup_rot1/tauf;
  St_rot2=taup_rot2/tauf;
  wparam_sim(time);
}
/* ----------------------------------------*/
void wparam_sim(double time){
  cout<<" --- Run Parameters at time="<<time<<"----------\n";
  cout<<"tauf="<<tauf<<"\n";
  cout<<"St_trans="<<St_trans<<"\n";
  cout<<"Re="<<Re<<"\n";
  cout<<"St_rot1="<<St_rot1<<"\n";
  cout<<"St_rot2="<<St_rot2<<"\n";
  cout<<"rsqr="<<rsqr<<"\n";
  cout<<" ------------------\n";
}
/* ----------------------------------------*/
