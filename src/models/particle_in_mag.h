#ifndef FILE_PTM_SEEN
#define FILE_PTM_SEEN
/*---------------------------------------*/
#include<math.h>
#include"modules/3vec.h"
#include "input.h";
using namespace std;
/* ----------------------------------------*/
unsigned int const pdim=6;
int Btype=2;
double Bzero=1.;
/* ----------------------------------------*/
struct Particle{
  /* State variables */
  vec3 xx; // position
  vec3 vv; //momentum
};
/* ----------------------------------------*/
double const kk=32.;
double const qbym=1;
unsigned int const ndim=pdim*Nensemble;
vec3 inivel(int velini);
vec3 inipos(int posini);
void array2p(Particle *PP,  double y[]);
void p2array(double y[], Particle *PP);
void iniconf(double y[]);
void getBB(vec3& BB, double GradB[3][3], vec3 xx);
void Diagnostics(int itn, double y[], double time);
void wparagm_sim(double time);
/* Diagnostic variables */
double yzero[ndim];
double rsqr,energy;

/* ----------------------------------------*/
void eval_rhs(double time,double y[],double rhs[]){
  vec3 BB,Force;
  double GradB[3][3];
  Particle PP;
  //
  array2p(&PP, &y[0]);
  getBB(BB, GradB, PP.xx);
//---Calculate Force---------
  Force=cross(PP.vv,BB)*qbym;
// ----------------------------------
  rhs[0]=PP.vv.x;
  rhs[1]=PP.vv.y;
  rhs[2]=PP.vv.z;
//
  rhs[3]=Force.x;
  rhs[4]=Force.y;
  rhs[5]=Force.z;
  //
}
/*----------------------------------------*/
void getBB(vec3& BB, double GradB[3][3], vec3 xx){
  double A_abc=1.;
  double B_abc=1.;
  double C_abc=1.;
  double Bconst=1.;
  switch(Btype){
  case(1):
    {
      A_abc=Bzero*A_abc;
      B_abc=Bzero*B_abc;
      C_abc=Bzero*C_abc;
      BB.x= A_abc*sin(kk*xx.z) + C_abc*cos(kk*xx.y);
      BB.y= B_abc*sin(kk*xx.x) + A_abc*cos(kk*xx.z);
      BB.z= C_abc*sin(kk*xx.y) +  B_abc*cos(kk*xx.x);
      GradB[0][0]=0;
      GradB[0][1]= -kk*C_abc*sin(kk*xx.y);
      GradB[0][2]= kk*A_abc*cos(kk*xx.z);
      //
      GradB[1][0]= kk*B_abc*cos(kk*xx.x);
      GradB[1][1]= 0.;
      GradB[1][2]= -kk*A_abc*sin(kk*xx.z);
      //
      GradB[2][0]= -kk*B_abc*sin(kk*xx.x);
      GradB[2][1]= kk*C_abc*cos(kk*xx.y);
      GradB[2][2]= 0.;
      break;}
  case(2):
    {
      BB.x= 0.;
      BB.y= 0.;
      BB.z= Bzero;
      GradB[0][0]=0;
      GradB[0][1]= 0.;
      GradB[0][2]= 0.;
      //
      GradB[1][0]= 0.;
      GradB[1][1]= 0.;
      GradB[1][2]= 0.;
      //
      GradB[2][0]= 0.;
      GradB[2][1]= 0.;
      GradB[2][2]= 0.;
      break;}
  default:
    exit(EXIT_FAILURE);
  }
};
/* ----------------------------------------*/
void Diagnostics(int itn, double y[], double time){
  Particle PP, PP0;
  vec3 BB;
  double GradB[3][3];
  if (itn==0){
    for(int iy=0;iy<ndim;iy++){
      yzero[iy]=y[iy];
    }
  };
  rsqr=0.;
  energy=0.;
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2p(&PP, &y[irb]);
    array2p(&PP0, &yzero[irb]);
    rsqr=rsqr+norm(PP.xx-PP0.xx);
    energy=energy+sqnorm(PP.vv);
  }
  //wparam_sim(time);
  //cout<<time<<"\t"<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<energy<<"\n";
  printf("%lf\t%lf\t%lf\t%lf\t\n",time,y[0],y[1],energy);
}
/* ----------------------------------------*/
void wparam_sim(double time){
  cout<<" --- Run Parameters at time="<<time<<"----------\n";
  cout<<"rsqr="<<rsqr<<"\n";
  cout<<"Energy="<<energy<<"\n";
}
/* ----------------------------------------*/
void iniconf(double y[]){
  Particle PP;
  int ibody,k;
  double pi;
  srand(time(NULL));
  for (ibody=0;ibody<Nensemble;ibody++){
    PP.xx=inipos(1);
    PP.vv=inivel(1);
    p2array(&y[ibody*pdim],&PP);
  }
}

/* ----------------------------------------*/
void array2p(Particle *PP,  double y[pdim]){
  /* 0,1,2 are position */
  PP->xx.x=y[0];
  PP->xx.y=y[1];
  PP->xx.z=y[2];
  /* 3,4,5 are momentum */
  PP->vv.x=y[3];
  PP->vv.y=y[4];
  PP->vv.z=y[5];
}
/* ----------------------------------------*/
void p2array(double y[pdim], Particle *PP){
  /* 0,1,2 are position */
  y[0]=PP->xx.x;
  y[1]=PP->xx.y;
  y[2]=PP->xx.z;
  /* 3,4,5 are momentum */
  y[3]=PP->vv.x;
  y[4]=PP->vv.y;
  y[5]=PP->vv.z;
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
    //v1=1.+rand();
    //v2=1.+rand();
    //v3=1.+rand();
    v1=0.;v2=1.;v3=0.;
    vv=vec3(v1,v2,v3);
    vv = vv*(1./norm(vv) );
    break;
  default:
    exit(EXIT_FAILURE);
  }
  return vv;
}
/* ----------------------------------------*/
#endif /* !FILE_PTM_SEEN */
