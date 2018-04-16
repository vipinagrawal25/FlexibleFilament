#ifndef FILE_VBD_SEEN
#define FILE_VBD_SEEN
/*---------------------------------------*/
#include<math.h>
#include <complex.h>
#include "Random.h";
#include "input.h";
using namespace std;
typedef complex<double> dcomp;
dcomp I=dcomp(0.,1.);
/* ----------------------------------------*/
unsigned int const pdim=4;
/* ----------------------------------------*/
struct Bfield{
  /* State variables */
  dcomp x; // x component of magnetic field
  dcomp y; // y component of magnetic field
};
/* ----------------------------------------*/
double const kk=1.;
double ksqr=kk*kk;
double const alpha0=1.;
double const tauc=10.;
double nu = 1./tauc;
double const Shear=4.;
double const eta=1.;
double Biamp = 1.e-5;
unsigned int const ndim=pdim*Nensemble;
void array2b(Bfield *B,  double y[]);
void b2array(double y[], Bfield *B);
void iniconf(double y[]);
void Diagnostics(int itn, double y[], double time);
void wparagm_sim(double time);
/* Diagnostic variables */
double yzero[ndim];
int const pmom=4;
double Bmom[pmom];
/* ----------------------------------------*/
double telegraph(double tauc, double time){
  int icall;
  if (time == 0){
    icall=0;
  }else{icall=1;}
  double poisson_mean=time/tauc;
  double prand=Poisson(poisson_mean,&icall);
  int ppower=(int) fmod(prand,2.);
  double tele_ran = pow(-1.,ppower);
  return tele_ran;
}
/* ----------------------------------------*/
void eval_rhs(double time,double y[],double rhs[]){
  Bfield B,dtB;
//

  array2b(&B, &y[0]);
// ----------------------------------
  double alpha=alpha0*telegraph(tauc,time);
  //dtB.x= I*kk*alpha*B.y-eta*ksqr*B.x;
  dtB.x= I*kk*alpha0*B.y-eta*ksqr*B.x;
  dtB.y= -Shear*B.x    -eta*ksqr*B.y;
// ----------------------------------
  b2array(&rhs[0],&dtB);
//
}
/* ----------------------------------------*/
void Diagnostics(int itn, double y[], double time){
  Bfield B;
  for (int i=0;i<pmom;i++){Bmom[i]=0.;}
  for (int iensemble=0;iensemble<Nensemble;iensemble++){
  array2b(&B, &y[iensemble*pdim]);
    for (int i=0;i<pmom;i++){
      Bmom[i]=Bmom[i]+pow(abs(real(B.x)),i+1);
    } 
  }  
  printf("meanB=%f\n",Bmom[0]);
   
}
/* ----------------------------------------*/
void wparam_sim(double time){

}
/* ----------------------------------------*/
void iniconf(double y[]){
  Bfield B;
  double gammad=eta*ksqr;
  double gamma_aS=sqrt(alpha0*kk*abs(Shear));
  // write out some parameters on first call 
  printf("nu=\t%f\n",nu);
  printf("gamma_d=\t%f\n",gammad);
  printf("gamma_aS=\t%f\n",gamma_aS);
  int icall=0;
  B.x=Biamp*dcomp(UniformRandom(&icall),UniformRandom(&icall));
  B.y=Biamp*dcomp(UniformRandom(&icall),UniformRandom(&icall));
  b2array(&y[0],&B);
}

/* ----------------------------------------*/
void array2b(Bfield *B,  double y[pdim]){
  /* real and imaginary part of Bx */
  B->x=dcomp(y[0],y[1]);
  /* real and imaginary part of By */
  B->y=dcomp(y[2],y[3]);
}
/* ----------------------------------------*/
void b2array(double y[pdim], Bfield *B){
  y[0]=real(B->x);
  y[1]=imag(B->x);
  y[2]=real(B->y);
  y[3]=imag(B->y);
}
/* ----------------------------------------*/
#endif /* !FILE_VBD_SEEN */
