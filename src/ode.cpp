#include <iostream>
#include<fstream>
#include<cmath>
#include "ode.h"
#include "model.h"
using namespace std;
/*********************************/
// void euler(unsigned int ndim, double *y, double time,double dt){
//   double  k1[ndim];
//   int idim;
//   eval_rhs(time,y,k1);
//   for(idim=0;idim<ndim;idim++){
//     y[idim]=y[idim]+k1[idim]*dt;  
//   }
// }
// /*********************************/
// void rnkt2(unsigned int ndim, double *y, double time,double dt){
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
// void rnkt4(unsigned int ndim, double *y, double time, double dt){
//   double  temp[ndim],k1[ndim],k2[ndim],k3[ndim],k4[ndim];
//   int idim;
//   eval_rhs(time,y,k1);
//   for(idim=0;idim<ndim;idim++){
//     temp[idim]=y[idim]+k1[idim]*dt/2.;
//   }
//   eval_rhs(time+(dt/2.),temp,k2);
//   for(idim=0;idim<ndim;idim++){
//     temp[idim]=y[idim]+k2[idim]*dt/2.;
//   }
//   eval_rhs(time+(dt/2.),temp,k3);
//   for(idim=0;idim<ndim;idim++){
//     temp[idim]=y[idim]+k3[idim]*dt;
//   }
//   eval_rhs(time+dt,temp,k4);
//   for(idim=0;idim<ndim;idim++){
//     y[idim]=y[idim]+dt*(  (k1[idim]/6.) + (k2[idim]/3.) + (k3[idim]/3.) + (k4[idim]/6.) );
//   }
// }

void rnkf45(unsigned int ndim, double *y, double time, double* add_dt, double* CurvSqr, double* SS, double ldiagnos)
{
	// Details of method: http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf

	double temp[ndim], k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], s;
	int idim ;
	double error = 0;
	double dt = *add_dt;
	double tol_dt = pow(10,-6)*dt;
  bool flag_kappa;

  if (ldiagnos)
  {
      flag_kappa = false;
  }
  else
  {
      flag_kappa = true;
  }

	eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  // CurvSqr_Store = CurvSqr;

	for(idim=0;idim<ndim;idim++)
	{
      temp[idim]=y[idim]+k1[idim]*dt/4.;
  }
  flag_kappa = false; 
  eval_rhs(time+(dt/4.),temp,k2, flag_kappa, CurvSqr, SS);

  for(idim=0;idim<ndim;idim++)
  {
    temp[idim]=y[idim]+(3.*k1[idim]+9.*k2[idim])*dt/32.;
 	}
 	eval_rhs(time+(3/8.)*dt,temp,k3, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim]+ (1932.*k1[idim]-7200.*k2[idim]+7296.*k3[idim])*dt/2197. ;
 	}
 	eval_rhs(time+(12/13.)*dt, temp, k4, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] + 439.*dt*k1[idim]/216. - 8.*k2[idim]*dt + 3680.*k3[idim]*dt/513. - 845.*k4[idim]*dt/4104. ;
 	}
 	eval_rhs(time+dt, temp, k5, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] - 8.*k1[idim]*dt/27. + 2.*k2[idim]*dt - 3544.*k3[idim]*dt/2565. + 1859.*k4[idim]*dt/4104. - 11.*k5[idim]*dt/40. ;
 	}
 	eval_rhs(time+dt/2., temp, k6, flag_kappa, CurvSqr, SS);


 	for (int idim = 0; idim < ndim; ++idim)
 	{
    temp[idim] = y[idim] + 25.*k1[idim]*dt/216. + 1408.*k3[idim]*dt/2565. + 2197.*k4[idim]*dt/4104 - k5[idim]*dt/5.;
 		// cout << temp[idim] << endl;
    y[idim] = y[idim]+ 16.*k1[idim]*dt/135. + 6656.*k3[idim]*dt/12825. + 28561.*k4[idim]*dt/56430. - 9.*k5[idim]*dt/50. + 2.*k6[idim]*dt/55.;

 		error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
 	}
    error = sqrt(error);
    // cout << error << endl;
 	s = 0.84*pow(tol_dt/error,0.25);

  // cout << error << endl;
    // cout << s << endl;
 	*add_dt = s*dt;
  // CurvSqr = CurvSqr_Store;
}	

void DP54(unsigned int ndim, double *y, double time, double* add_dt, double* CurvSqr, double* SS, double ldiagnos)
{
  // In this function I have implemented Dormand-Prince Method which is more suitable than rkf45 for high order integration.
  // Details could be found in Numerical recipes book and a short description on the link: 
  // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
  double temp[ndim], s;
  int idim ;
  double error = 0;
  double dt = *add_dt;
  double tol_dt = pow(10,-2);
  bool flag_kappa;

  double ci[6] = {0,0.2,0.3,0.8,8./9,1} ;
  double aij[6][6] = {
    {0,0,0,0,0,0},
    {0.2,0,0,0,0,0},
    {3./40,9./40,0,0,0,0},
    {44./45,-56./15,32./9,0,0,0},
    {19372/6561,-25360./2187,64448./6561,-212./729,0,0},
    {9017./3168,-355./33,46732./5247,49./176,-5103./18656,0}
  };
  double bi[6] = {35./384,0,500./1113,125./192,-2187./6784,11./84};
  double bistar[7] = {5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};

  double k[7][ndim], k_temp[ndim];

  // eval_rhs(time,y,k1,flag_kappa);

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < ndim; ++j)
    {
      k[i][ndim] = 0;
    }
  }

  for (int j = 0; j < 6; ++j)
  {
    for (int idim = 0; idim < ndim; ++idim)
    {
      temp[idim] = y[idim];
      for (int i = 0; i < j; ++i)
      {
        temp[idim] = temp[idim] + aij[j][i]*dt*k[i][idim];       
      }
    }

    eval_rhs(time+dt*ci[j],temp,k_temp,flag_kappa,CurvSqr,SS);

    for (int idim = 0; idim < ndim; ++idim)
    {
        k[j][idim] = k_temp[idim];
    }

    if (j == 0)
    {
        flag_kappa = false;
    }

  }

  for (int idim = 0; idim < ndim; ++idim)
  {   
      for (int j = 0; j < 6; ++j)
      {
          y[idim] = y[idim] + bi[j]*k[j][idim]*dt;
      }
  }
  
  eval_rhs(time+dt,temp,k_temp,flag_kappa,CurvSqr,SS);

  for (int idim = 0; idim < ndim; ++idim)
  {
      k[6][idim] = k_temp[idim];
  }

  for (int idim = 0; idim < ndim; ++idim)
  {
      temp[idim] = y[idim];
      for (int i = 0; i < 7; ++i)
      {
          temp[idim] = temp[idim] + bistar[i]*k[i][idim]*dt;
      }

      error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
  }

  error = sqrt(error/ndim);

  s = 0.87*pow(tol_dt/error,0.20);
  // // cout << error << endl;
  // if (s > 10)
  // {
  //     s = 10;
  // }
  // else if (s<0.2)
  // {
  //     s = 0.2;
  // }

  *add_dt = s*dt;

  // cout << *add_dt << endl;

}


// void DP54(unsigned int ndim, double *y, double time, double* add_dt, double* CurvSqr, double ldiagnos)
// {
//   // In this function I have implemented Dormand-Prince Method which is more suitable than rkf45 for high order integration.
//   // Details could be found in Numerical recipes book and a short description on the link: 
//   // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
//   double temp[ndim], s;
//   int idim ;
//   double error = 0;
//   double dt = *add_dt;
//   double tol_dt = pow(10,-2)*dt;
//   bool flag_kappa;

//   double ci[6] = {0,0.2,0.3,0.8,8./9,1} ;
//   double aij[6][6] = {
//     {0,0,0,0,0,0},
//     {0.2,0,0,0,0,0},
//     {3./40,9./40,0,0,0,0},
//     {44./45,-56./15,32./9,0,0,0},
//     {19372/6561,-25360./2187,64448./6561,-212./729,0,0},
//     {9017./3168,-355./33,46732./5247,49./176,-5103./18656,0}
//   };
//   double bi[6] = {35./384,0,500./1113,125./192,-2187./6784,11./84};
//   double bistar[7] = {5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};

//   double k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], k7[ndim], ynew[ndim];

//   // eval_rhs(time,y,k1,flag_kappa);

//   for (int i = 0; i < 8; ++i)
//   {
//     for (int j = 0; j < ndim; ++j)
//     {
//       k[i][ndim] = 0;
//     }
//   }

//   for (int j = 0; j < 6; ++j)
//   {
//     for (int idim = 0; idim < ndim; ++idim)
//     {
//       temp[idim] = y[idim];
//       for (int i = 0; i < j; ++i)
//       {
//         temp[idim] = temp[idim] + aij[j][i]*dt*k[i][idim];       
//       }
//     }

//     eval_rhs(time+dt*ci[j],temp,k_temp,flag_kappa,CurvSqr);

//     for (int idim = 0; idim < ndim; ++idim)
//     {
//         k[j][idim] = k_temp[idim];
//     }

//     if (j == 0)
//     {
//         flag_kappa = false;
//     }

//   }

//   for (int idim = 0; idim < ndim; ++idim)
//   {   
//       y[idim] = 0;
//       for (int j = 0; j < 6; ++j)
//       {
//           y[idim] = y[idim] + bi[j]*k[j][idim];
//       }
//   }
  
//   eval_rhs(time+dt,temp,k_temp,flag_kappa,CurvSqr);

//   for (int idim = 0; idim < ndim; ++idim)
//   {
//       k[6][idim] = k_temp[idim];
//   }

//   for (int idim = 0; idim < ndim; ++idim)
//   {
//       temp[idim] = y[idim];
//       for (int i = 0; i < 7; ++i)
//       {
//           temp[idim] = temp[idim] + bistar[i]*k[i][idim];
//       }

//       error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
//   }

//   error = sqrt(error/ndim);

//   s = 0.87*pow(tol_dt/error,0.20);
//   // cout << error << endl;
//   if (s > 10)
//   {
//       s = 10;
//   }
//   else if (s<0.2)
//   {
//       s = 0.2;
//   }

//   *add_dt = s*dt;

//   // cout << *add_dt << endl;

// }


/* call eval_rhs(y,t,k1)
call eval_rhs(y+k1*(dt/2.),t+(dt/2.),k2)
call eval_rhs(y+k2*(dt/2.),t+(dt/2.),k3)
call eval_rhs(y+k3*dt,t+dt,k4)
y=y+dt*((k1/6.)+(k2/3.)+k3/3.+k4/6.) */
/* ----------------- */

