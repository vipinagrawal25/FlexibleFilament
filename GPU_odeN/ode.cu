#include <iostream>
#include<fstream>
#include "ode.h"
#include "param.h"
using namespace std;
/*********************************/
__global__ void pre_euler( void ){
  cudaMalloc( (void**)&temp, ndim*sizeof(double) );
}
/*********************************/
__global__ void euler(double *psi, double time,double dt){
  
  eval_rhs(time,&y[0],temp);
  for(idim=0;idim<=ndim-1;idim++){
    y[idim]=y[idim]+temp[idim]*dt;
  }
}
/*********************************
void rnkt4(unsigned int ndim, double *y, double time,double dt){
  double  temp[ndim],k1[ndim],k2[ndim],k3[ndim],k4[ndim];
  int idim;
    //k1=f(y,t)
  eval_rhs(time,y,k1);
  for(idim=0;idim<=ndim-1;idim++){
    temp[idim]=y[idim]+k1[idim]*dt/2.;
  }
  eval_rhs(time+(dt/2.),temp,k2);
  for(idim=0;idim<=ndim-1;idim++){
    temp[idim]=y[idim]+k2[idim]*(dt/2.);
  }
  eval_rhs(time+(dt/2.),temp,k3);
  for(idim=0;idim<=ndim-1;idim++){
    temp[idim]=y[idim]+k3[idim]*(dt/2.);
  }
  eval_rhs(time+dt,temp,k4);
  for(idim=0;idim<=ndim-1;idim++){
    y[idim]=y[idim]+dt*((k1[idim]/6.)+(k2[idim]/3.)+k3[idim]/3.+k4[idim]/6.);
  }
}*/
/* ----------------- */

