#include <iostream>
#include<fstream>
#include<cmath>
#include "ode.h"
#include "model.h"
#include "constant.h"
using namespace std;
/*********************************/
// void euler(double *y,double* vel, double time,double dt){
//   double  k1[ndim];
//   int idim;
//   eval_rhs(time,y,k1);
//   for(idim=0;idim<ndim;idim++){
//     y[idim]=y[idim]+k1[idim]*dt;  
//   }
// }
// // /*********************************/
// void rnkt2(double *y,double* vel,double time,double dt){
//   double  temp[ndim],k1[ndim],k2[ndim];
//   int idim;
// //
//   eval_rhs(time,y,k1);
//   for(idim=0;idim<ndim;idim++){
//     temp[idim]=y[idim]+k1[idim]*dt/2.;
//   }
//   eval_rhs(time+(dt/2.),temp,k2);
//   for(idim=0;idim<ndim;idim++){
//     y[idim]=y[idim]+k2[idim]*dt;
//   }
// }
// /*********************************/
// Wrapper function for rnkt4.
void rnkt4(double *y,double *add_time, double* add_dt){
  double CurvSqr[ndim],SS[ndim],vel[ndim];
  double ldiagnos=0;
  // CurvSqr[0]=NULL;
  // SS[0]=NULL;
  rnkt4(ndim,&y[0],&vel[0],add_time,add_dt,&CurvSqr[0],&SS[0],ldiagnos);
}
/*********************************/
void rnkt4(unsigned int ndim, double *y, double *vel, double *add_time, double *add_dt, double* CurvSqr, double* SS, 
          double ldiagnos){
  double temp[ndim],k1[ndim],k2[ndim],k3[ndim],k4[ndim];
  int idim;
  double dt = *add_dt;
  double time = *add_time;
  bool flag_kappa;

  if (ldiagnos){
      flag_kappa = false;
  }
  else{
      flag_kappa = true;
  }

  eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k1[idim]*dt/2.;
  }
  flag_kappa = false;
  eval_rhs(time+(dt/2.),temp,k2,flag_kappa,CurvSqr,SS);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k2[idim]*dt/2.;
  }
  eval_rhs(time+(dt/2.),temp,k3,flag_kappa,CurvSqr,SS);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k3[idim]*dt;
  }
  eval_rhs(time+dt,temp,k4,flag_kappa,CurvSqr,SS);
  for(idim=0;idim<ndim;idim++){
    y[idim]=y[idim]+dt*(  (k1[idim]/6.) + (k2[idim]/3.) + (k3[idim]/3.) + (k4[idim]/6.) );
    vel[idim]=k1[idim];
  }
  *add_time = time + dt;
}
/*********************************/
// Wrapper function for rnkf45.
void rnkf45(double *y, double *add_time, double* add_dt){
  double CurvSqr[ndim],SS[ndim],vel[ndim];
  double ldiagnos=0;
  // CurvSqr[0]=NULL;
  // SS[0]=NULL;
  rnkf45(ndim,&y[0],&vel[0],add_time,add_dt,&CurvSqr[0],&SS[0],ldiagnos);
}
/*********************************/
void rnkf45(unsigned int ndim, double *y, double *vel, double *add_time, double* add_dt, double* CurvSqr, 
            double* SS, double ldiagnos){
  // Details of method: http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
  // add_time is the address of time and the same goes for dt as well.
  double temp[ndim], k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], s, ynew[ndim];
  int idim ;
  double error,temp_error;
  double dt = *add_dt;
  // double tol_dt = pow(10,-9)*dt;
  double tol_dt = 1.e-6;
  bool flag_kappa;
  double time = *add_time;
  double epsilon = 0.84;
  double truncationmax=2.0;
  double truncationmin=0.5;
  // original fehlberg parameters
  // double ci[6] = {0,0.25,3./8,12./13,1.,1./2};
  // double aij[6][5] = {
  //   {0,0,0,0,0},
  //   {0.25,0,0,0,0},
  //   {3./32.,9./32.,0,0,0},
  //   {1932./2197.,-7200./2197.,7296./2197.,0,0},
  //   {439./216.,-8.,3680./513.,-845./4104.,0},
  //   {-8./27.,2.,-3544./2565.,1859./4104.,-11./40.}
  // };
  // double bistar[6] = {16./135.,0,6656./12825,28561./56430,-9./50,2./55};
  // double bi[6] = {25./216.,0,1408./2565.,2197./4104.,-1./5.,0};

  //Cash-Karp Parameters for evolution
  double ci[6] = {0,0.2,0.3,0.6,1.,7./8} ;
  double aij[6][5] = {
    {0,0,0,0,0},
    {0.2,0,0,0,0},
    {3./40.,9./40.,0,0,0},
    {3./10.,-9./10.,6./5.,0,0},
    {-11./54.,5/2.,-70./27,35./27.,0},
    {1631./55296,175./512.,575./13824.,44275./110592,253./4096}
  };
  double bistar[6] = {37./378,0,250./621,125./594,0,512./1771};
  double bi[6] = {2825./27648.,0,18575./48384,13525./55296.,277./14336.,0.25};
  if (ldiagnos){
      flag_kappa = true;
  }
  else{
      flag_kappa = false;
  }
  eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  // CurvSqr_Store = CurvSqr;
  for(idim=0;idim<ndim;idim++){
      temp[idim]=y[idim]+k1[idim]*dt*aij[1][0];
  }

  flag_kappa = false;
  eval_rhs(time+dt*ci[1],temp,k2, flag_kappa, CurvSqr, SS);
  
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+(aij[2][0]*k1[idim]+aij[2][1]*k2[idim])*dt;
  }
  eval_rhs(time+ci[2]*dt,temp,k3, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim]+ (aij[3][0]*k1[idim]+aij[3][1]*k2[idim]+aij[3][2]*k3[idim])*dt ;
  }
  eval_rhs(time+ci[3]*dt, temp, k4, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + aij[4][0]*dt*k1[idim] + aij[4][1]*k2[idim]*dt + aij[4][2]*k3[idim]*dt + 
                  aij[4][3]*k4[idim]*dt ;
  }
  eval_rhs(time+ci[4]*dt, temp, k5, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + aij[5][0]*k1[idim]*dt + aij[5][1]*k2[idim]*dt + aij[5][2]*k3[idim]*dt
                         + aij[5][3]*k4[idim]*dt + aij[5][4]*k5[idim]*dt ;
  }
  eval_rhs(time+dt*ci[5], temp, k6, flag_kappa, CurvSqr, SS);
  error=0;
  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + dt*(bi[0]*k1[idim] + bi[1]*k2[idim] + bi[2]*k3[idim] + bi[3]*k4[idim] + bi[4]*k5[idim] 
                         +  bi[5]*k6[idim]);
    // cout << temp[idim] << endl;
    ynew[idim] = y[idim]+ dt*(bistar[0]*k1[idim] + bistar[1]*k2[idim] + bistar[2]*k3[idim] + bistar[3]*k4[idim] 
                        + bistar[4]*k5[idim] +  bistar[5]*k6[idim]);
    temp_error = abs(temp[idim]-ynew[idim]);
    error=max(temp_error,error);
  }
  error=error+tiny;
  
  // cout << error << endl;

  if (error<tol_dt){
    *add_time=time+dt;
    for (int idim = 0; idim < ndim; ++idim){
      y[idim]=ynew[idim];   // Accept the step
      vel[idim]=k1[idim];
    }
    s = epsilon*pow((tol_dt/error),0.20);
    if (s>truncationmax){s=truncationmax;}
    *add_dt = s*dt;
    // cout << *add_dt << endl;
  }else{
    s = epsilon*pow((tol_dt/error),0.25);
    if (s<truncationmin){s=truncationmin;}
    *add_dt = s*dt;
    rnkf45(ndim, &y[0], &vel[0],add_time, add_dt, &CurvSqr[0], &SS[0], ldiagnos);
  }
} 
/*********************************/
void DP54(unsigned int ndim, double *y, double *vel, double *add_time, double* add_dt, double* CurvSqr, double* SS, 
          double ldiagnos)
{
  // In this function I have implemented Dormand-Prince Method which is more suitable than rkf45 for high order integration.
  // Details could be found in Numerical recipes book and a short description on the link: 
  // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
  double temp[ndim], s,k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], k7[ndim], ynew[ndim];
  double error,temp_error;
  int idim ;
  double dt = *add_dt;
  double tol_dt = pow(10,-5);
  bool flag_kappa;
  double time = *add_time;
  double epsilon=0.87;
  double truncationmax=2;
  double truncationmin=0.5;

  double ci[7] = {0,0.2,0.3,0.8,8./9,1.,1.} ;
  double aij[7][6] = {
    {0,0,0,0,0,0},
    {0.2,0,0,0,0,0},
    {3./40,9./40,0,0,0,0},
    {44./45,-56./15,32./9,0,0,0},
    {19372/6561,-25360./2187,64448./6561,-212./729,0,0},
    {9017./3168,-355./33,46732./5247,49./176,-5103./18656,0},
    {35./384,0,500./1113,125./192,-2187./6784,11./84}
  };
  double bi[7] = {5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
  // eval_rhs(time,y,k1,flag_kappa);
  if (ldiagnos){
      flag_kappa = true;
  }
  else{
      flag_kappa = false;
  }
  eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  // CurvSqr_Store = CurvSqr;

  for(idim=0;idim<ndim;idim++){
      temp[idim]=y[idim]+k1[idim]*dt*aij[1][0];
  }
  flag_kappa = false; 
  eval_rhs(time+dt*ci[1],temp,k2, flag_kappa, CurvSqr, SS);

  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+(aij[2][0]*k1[idim]+aij[2][1]*k2[idim])*dt;
  }
  eval_rhs(time+ci[2]*dt,temp,k3,flag_kappa,CurvSqr,SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim]+ (aij[3][0]*k1[idim]+aij[3][1]*k2[idim]+aij[3][2]*k3[idim])*dt ;
  }
  eval_rhs(time+ci[3]*dt, temp, k4, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + aij[4][0]*dt*k1[idim] + aij[4][1]*k2[idim]*dt + aij[4][2]*k3[idim]*dt + aij[4][3]*k4[idim]*dt ;
  }
  eval_rhs(time+ci[4]*dt, temp, k5, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + aij[5][0]*k1[idim]*dt + aij[5][1]*k2[idim]*dt + aij[5][2]*k3[idim]*dt + aij[5][3]*k4[idim]*dt 
                         + aij[5][4]*k5[idim]*dt ;
  }
  eval_rhs(time+dt*ci[5], temp, k6, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim){
    ynew[idim] = y[idim] + aij[6][0]*k1[idim]*dt + aij[6][1]*k2[idim]*dt + aij[6][2]*k3[idim]*dt + aij[6][3]*k4[idim]*dt 
                          + aij[6][4]*k5[idim]*dt + aij[6][5]*k6[idim]*dt  ;
  }
  eval_rhs(time+dt*ci[6], ynew, k7, flag_kappa, CurvSqr, SS);

  error= 0;
  temp_error=0;
  for (int idim = 0; idim < ndim; ++idim){
    temp[idim] = y[idim] + bi[0]*k1[idim]*dt + bi[2]*k3[idim]*dt + bi[3]*k4[idim]*dt + bi[4]*k5[idim]*dt 
                            + bi[5]*k6[idim]*dt + bi[6]*k7[idim]*dt;
    temp_error =  abs(temp[idim]-y[idim]);
    error = max(temp_error,error);
  }
  error=error+tiny;
  // cout << error << endl;
  if (error<tol_dt){
    *add_time=time+dt;
    for (int idim = 0; idim < ndim; ++idim){
      y[idim]=ynew[idim];   // Accept the step
      vel[idim]=k1[idim];
    }
    s = epsilon*pow((tol_dt/error),0.25);
    if (s>truncationmax){s=truncationmax;}
    *add_dt = s*dt;
    // cout << *add_dt << endl;
  }else{
    s = epsilon*pow((tol_dt/error),0.2);
    if (s<truncationmin){s=truncationmin;}
    *add_dt = s*dt;
    // cout << "So you mean to say that the segmentation fault is here?" << endl;
    DP54(ndim, &y[0], &vel[0],add_time, add_dt, &CurvSqr[0], &SS[0], ldiagnos);
  }
}
/* call eval_rhs(y,t,k1)
call eval_rhs(y+k1*(dt/2.),t+(dt/2.),k2)
call eval_rhs(y+k2*(dt/2.),t+(dt/2.),k3)
call eval_rhs(y+k3*dt,t+dt,k4)
y=y+dt*((k1/6.)+(k2/3.)+k3/3.+k4/6.) */
/* ----------------- */
