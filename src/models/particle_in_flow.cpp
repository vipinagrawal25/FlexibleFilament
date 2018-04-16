#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "pif.h"
using namespace std;
/**************************/
unsigned int const ndim=pdim*Nensemble;
int yini=1;
void iniconf(double y[]);
void normalize(double AA[]);
void getUU(double xx[],double BB[]);
void cross(double A[],double B[], double C[]);
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double time=0.;
  iniconf(y);
//----------------------------
  while(time <= TMAX){
    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,y[0],y[1],y[2],y[3],y[4],y[5]);
    for(int k=0;k<=ndim-1;k+=pdim){
      rnkt4(pdim,&y[k],time,dt);
    //euler(ndim,&y[0],time,dt);
    }
    time=time+dt;
    //    cout << y[0] << "\n";
  }
}
/* ----------------------------------------*/
void eval_rhs(double time,double y[],double rhs[]){
  double UU[pdim/2];
  int i;
  double *xx,*vv;
  xx=&y[0];
  vv=&y[pdim/2];
//-----------------
  getUU(xx,UU);
  /*printf("%s\t%f\t%f\t%f\n","xx",y[0],y[1],y[2]);
    printf("%s\t%f\t%f\t%f\n","vv",vv[0],vv[1],vv[2]); */
  //  printf("%s\t%f\t%f\t%f\t%f\n","vXB",time,y[3],y[4],y[5]); 
  for(i=0;i<=pdim/2-1;i++){
    //printf("%d\t%d\n",i,i+3);
    rhs[i]=vv[i] + beta*UU[i];
    rhs[i+3]=(1/St)*( (1-beta)*UU[i]-vv[i] );
  }
}
/* ----------------------------------------*/
void iniconf(double y[]){
  int k;
  double *xx,*vv;
  double pi;
  srand(time(NULL));
  switch(yini){
  case 1:
    pi=4.*atan(1.);
    for(k=0;k<=ndim-1;k+= pdim){
      xx=&y[k];
      vv=&y[k+3];
      xx[0]=rand();
      xx[1]=rand();
      xx[2]=rand();
      vv[0]=rand();
      vv[1]=rand();
      vv[2]=rand();
      normalize(xx);
      xx[0]=2*pi*xx[0];
      xx[1]=2*pi*xx[1];
      xx[2]=2*pi*xx[2];
      normalize(vv);
      printf("%d\t%f\t%f\n",k,xx[0],xx[1]);
    }
    break;
  case 2:
    for(k=0;k<=ndim-1;k+= pdim){
      xx=&y[k];
      vv=&y[k+3];
      xx[0]=0;
      xx[1]=1.;
      xx[2]=0;
      vv[0]=1.;
      vv[1]=0.;
      vv[2]=0.;
      normalize(xx);
      normalize(vv);
      printf("%d\t%f\t%f\t%f\n",k,vv[0],vv[1],vv[2]);
    }
    break;
  default:
    exit(EXIT_FAILURE);
  }
}
/* ----------------------------------------*/
void normalize(double AA[]){
  double norm;
  norm = sqrt(AA[0]*AA[0]+AA[1]*AA[1]+AA[2]*AA[2]);
  AA[0]=AA[0]/norm;
  AA[1]=AA[1]/norm;
  AA[2]=AA[2]/norm;
}
/* ----------------------------------------*/
void getUU(double xx[],double UU[]){
  double A_abc=1.;
  double B_abc=1.;
  double C_abc=1.;
  switch(Utype){
  case(1):
    UU[0]=A_abc*sin(kk*xx[2]) + C_abc*cos(kk*xx[1]);
    UU[1]=B_abc*sin(kk*xx[0]) + A_abc*cos(kk*xx[2]);
    UU[2]=C_abc*sin(kk*xx[1]) + B_abc*cos(kk*xx[0]);
    break;
  case(2):
    break;
  default:
    exit(EXIT_FAILURE);
  }
}
/* ----------------------------------------*/
void cross(double AA[],double BB[], double CC[]){
  CC[0] = AA[1]*BB[2]-AA[2]*BB[1];
  CC[1] = -AA[0]*BB[2]+AA[2]*BB[0];
  CC[2] = AA[0]*BB[1]-AA[1]*BB[0];
}
