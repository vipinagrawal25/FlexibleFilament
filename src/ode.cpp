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
void rnkt4(unsigned int ndim, double *y, double *add_time, double *add_dt, double* CurvSqr, double* SS, double ldiagnos)
{
  double temp[ndim],k1[ndim],k2[ndim],k3[ndim],k4[ndim];
  int idim;
  double dt = *add_dt;
  double time = *add_time;
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
  }

  *add_time = time + dt;
}

void rnkf45(unsigned int ndim, double *y, double *add_time, double* add_dt, double* CurvSqr, double* SS, double ldiagnos)
{
	// Details of method: http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
  // add_time is the address of time and the same goes for dt as well.

	double temp[ndim], k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], s, yold[ndim];
	int idim ;
	double error = 0;
	double dt = *add_dt;
	// double tol_dt = pow(10,-9)*dt;
  double tol_dt = pow(10,-8);
  bool flag_kappa;
  double time = *add_time;

  double ci[6] = {0,0.25,3./8,12./13,1.,1./2} ;
  double aij[6][5] = {
    {0,0,0,0,0},
    {0.25,0,0,0,0},
    {3./32.,9./32.,0,0,0},
    {1932./2197.,-7200./2197.,7296./2197.,0,0},
    {439./216.,-8.,3680./513.,-845./4104.,0},
    {-8./27.,2.,-3544./2565.,1859./4104.,-11./40.}
  };
  double bistar[6] = {16./135.,0,6656./12825.,28561./56430.,-9./50.,2./55.};
  double bi[6] = {25./216.,0,1408./2565.,2197./4104.,-1./5.,0};

  if (ldiagnos)
  {
      flag_kappa = true;
      // cout<< "This is a high level shit" << endl;
  }
  else
  {
      flag_kappa = false;
  }

	eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  // CurvSqr_Store = CurvSqr;

	for(idim=0;idim<ndim;idim++)
	{
      temp[idim]=y[idim]+k1[idim]*dt*aij[1][0];
  }
  flag_kappa = false; 
  eval_rhs(time+dt*ci[1],temp,k2, flag_kappa, CurvSqr, SS);

  for(idim=0;idim<ndim;idim++)
  {
    temp[idim]=y[idim]+(aij[2][0]*k1[idim]+aij[2][1]*k2[idim])*dt;
 	}
 	eval_rhs(time+ci[2]*dt,temp,k3, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim]+ (aij[3][0]*k1[idim]+aij[3][1]*k2[idim]+aij[3][2]*k3[idim])*dt ;
 	}
 	eval_rhs(time+ci[3]*dt, temp, k4, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] + aij[4][0]*dt*k1[idim] + aij[4][1]*k2[idim]*dt + aij[4][2]*k3[idim]*dt + aij[4][3]*k4[idim]*dt ;
 	}
 	eval_rhs(time+ci[4]*dt, temp, k5, flag_kappa, CurvSqr, SS);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] + aij[5][0]*k1[idim]*dt + aij[5][1]*k2[idim]*dt + aij[5][2]*k3[idim]*dt + aij[5][3]*k4[idim]*dt + aij[5][4]*k5[idim]*dt ;
 	}
 	eval_rhs(time+dt*ci[5], temp, k6, flag_kappa, CurvSqr, SS);


 	for (int idim = 0; idim < ndim; ++idim)
 	{
    temp[idim] = y[idim] + bi[0]*k1[idim]*dt + bi[2]*k3[idim]*dt + bi[3]*k4[idim]*dt + bi[4]*k5[idim]*dt;
 		// cout << temp[idim] << endl;
    
    yold[idim] = y[idim];
    y[idim] = y[idim]+ bistar[0]*k1[idim]*dt + bistar[2]*k3[idim]*dt + bistar[3]*k4[idim]*dt + bistar[4]*k5[idim]*dt + bistar[5]*k6[idim]*dt;

 		error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
 	}
    error = sqrt(error);
    // error = error/ndim;
    // cout << error << endl;
 	s = 0.84*pow(tol_dt/error,0.25);

  // cout << s << endl;

  if (s>10)
  {
    // cout << s << endl;
    if (isinf(s))
    {
        s=1;
    }
    else
    {
        s=10;
    }    
    *add_time = time + dt;
    *add_dt = s*dt;
  }
  else if (s < 0.5)
  {
    *add_dt = s*dt;
    if (s<0.2)
    {
       s = 0.2; 
       rnkf45(pdim, &yold[0], add_time, add_dt, &CurvSqr[0], &SS[0], ldiagnos);
   }
    // cout << error << endl;
  }
  else
  {
    *add_time = time + dt;
    *add_dt = s*dt;
  }

  // *add_time = time + dt;
  // *add_dt = s*dt;

  // cout << error << endl;
    // cout << s << endl;
  // CurvSqr = CurvSqr_Store;

  // cout << error << endl;
}	

// void rnkf451(unsigned int ndim, double *y, double *add_time, double* add_dt, double* CurvSqr, double* SS, double ldiagnos)
// {
//   // In this function I have implemented Runge Kutta fehlberg Method which is more suitable than rkf45 for high order integration.
//   // Details could be found in Numerical recipes book and a short description on the link: 
//   // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
//   double temp[ndim], s, yold[ndim];
//   int idim ;  
//   double error = 0;
//   double dt = *add_dt;
//   double tol_dt = pow(10,-10)*dt;
//   bool flag_kappa;
//   double time = *add_time;

//   double ci[6] = {0,0.25,3./8,12./13,1.,1./2} ;
//   double aij[6][5] = {
//     {0,0,0,0,0},
//     {0.25,0,0,0,0},
//     {3./32.,9./32.,0,0,0},
//     {1932./2197.,-7200./2197.,7296./2197.,0,0},
//     {439./216.,-8.,3680./513.,-845./4104.,0},
//     {-8./27.,2.,-3544./2565.,1859./4104.,-11./40.}
//   };
//   double bistar[6] = {16./135.,0,6656./12825.,28561./56430.,-9./50.,2./55.};
//   double bi[6] = {25./216.,0,1408./2565.,2197./4104.,-1./5.,0};
    
//   double k[6][ndim], k_temp[ndim];
//   // eval_rhs(time,y,k1,flag_kappa);

//   if (ldiagnos)
//   {
//       flag_kappa = true;
//   }
//   else
//   {
//       flag_kappa = false;
//   }

//   for (int i = 0; i < 7; ++i)
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
//       for (int i = 0; i < 5; ++i)
//       {
//         temp[idim] = temp[idim] + aij[j][i]*dt*k[i][idim];       
//       }
//     }

//     eval_rhs(time+dt*ci[j],temp,k_temp,flag_kappa,CurvSqr,SS);

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
//       // temp[idim] = y[idim];
//       // for (int j = 0; j < 6; ++j)
//       // {
//       //     temp[idim] = temp[idim] + bi[j]*k[j][idim]*dt;
//       // }
//        temp[idim] = y[idim] + 25.*k[0][idim]*dt/216. + 1408.*k[2][idim]*dt/2565. + 2197.*k[3][idim]*dt/4104 - k[4][idim]*dt/5.;      
//   }
  
//   // eval_rhs(time+dt,temp,k_temp,flag_kappa,CurvSqr,SS);

//   // for (int idim = 0; idim < ndim; ++idim)
//   // {
//   //     k[6][idim] = k_temp[idim];
//   // }

//   for (int idim = 0; idim < ndim; ++idim)
//   {
//       for (int i = 0; i < 6; ++i)
//       {
//           y[idim] = y[idim] + bistar[i]*k[i][idim]*dt;     
//       }

//       // y[idim] = y[idim]+ 16*k[0][idim]*dt/135. + 6656.*k[2][idim]*dt/12825. + 28561.*k[3][idim]*dt/56430. - 9.*k[4][idim]*dt/50. + 2.*k[5][idim]*dt/55.;
//       error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);       

//   }

//   error = sqrt(error/ndim);

//   s = 0.84*pow(tol_dt/error,0.2);
//   // cout << error << endl;
//   // if (s > 10)
//   // {
//   //     s = 10;
//   // }
//   // else if (s<0.2)
//   // {
//   //     s = 0.2;
//   // }
//   // cout << s << endl;
//   *add_time = time + dt;
//   *add_dt = s*dt;

//   // If the time stepping in next iteration is lesser than half of the current time step then the current step should be 
//   // repeated
//   // if (s<1)
//   // {
//   //   /* code */
//   // }

//   // cout << *add_dt << endl;

// }


void DP54(unsigned int ndim, double *y, double *add_time, double* add_dt, double* CurvSqr, double* SS, double ldiagnos)
{
  // In this function I have implemented Dormand-Prince Method which is more suitable than rkf45 for high order integration.
  // Details could be found in Numerical recipes book and a short description on the link: 
  // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
  double temp[ndim], s,k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], k7[ndim], yold[ndim];
  double err;
  int idim ;
  double dt = *add_dt;
  double tol_dt = pow(10,-4)*dt;
  bool flag_kappa;
  double time = *add_time;

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

  // eval_rhs(time,y,k1,flag_kappa);

  if (ldiagnos)
  {
      flag_kappa = true;
  }
  else
  {
      flag_kappa = false;
  }

  eval_rhs(time,y,k1,flag_kappa,CurvSqr,SS);
  // CurvSqr_Store = CurvSqr;

  for(idim=0;idim<ndim;idim++)
  {
      temp[idim]=y[idim]+k1[idim]*dt*aij[1][0];
  }
  flag_kappa = false; 
  eval_rhs(time+dt*ci[1],temp,k2, flag_kappa, CurvSqr, SS);

  for(idim=0;idim<ndim;idim++)
  {
    temp[idim]=y[idim]+(aij[2][0]*k1[idim]+aij[2][1]*k2[idim])*dt;
  }
  eval_rhs(time+ci[2]*dt,temp,k3, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim)
  {
    temp[idim] = y[idim]+ (aij[3][0]*k1[idim]+aij[3][1]*k2[idim]+aij[3][2]*k3[idim])*dt ;
  }
  eval_rhs(time+ci[3]*dt, temp, k4, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim)
  {
    temp[idim] = y[idim] + aij[4][0]*dt*k1[idim] + aij[4][1]*k2[idim]*dt + aij[4][2]*k3[idim]*dt + aij[4][3]*k4[idim]*dt ;
  }
  eval_rhs(time+ci[4]*dt, temp, k5, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim)
  {
    temp[idim] = y[idim] + aij[5][0]*k1[idim]*dt + aij[5][1]*k2[idim]*dt + aij[5][2]*k3[idim]*dt + aij[5][3]*k4[idim]*dt + aij[5][4]*k5[idim]*dt ;
  }
  eval_rhs(time+dt*ci[5], temp, k6, flag_kappa, CurvSqr, SS);

  for (int idim = 0; idim < ndim; ++idim)
  {
    yold[idim] = y[idim];
  
    y[idim] = yold[idim] + bi[0]*k1[idim]*dt + bi[2]*k3[idim]*dt + bi[3]*k4[idim]*dt + bi[4]*k5[idim]*dt + bi[5]*k6[idim]*dt;
    // cout << temp[idim] << endl;
  }
  eval_rhs(time+dt,y,k7,flag_kappa,CurvSqr,SS);

  err = 0;
  for (int idim = 0; idim < ndim; ++idim)
  {
    temp[idim] = yold[idim] +  bistar[0]*k1[idim]*dt + bistar[2]*k3[idim]*dt + bistar[3]*k4[idim]*dt + bistar[4]*k5[idim]*dt + bistar[5]*k6[idim]*dt + bistar[6]*k7[idim]*dt;
    err = err + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
  }

  err = sqrt(err)/ndim;

  s = 0.87*pow(tol_dt/err,0.20);
  // cout << error << endl;
  // if (s > 10)
  // {
  //     s = 10;
  // }
  // else if (s<0.2)
  // {
  //     s = 0.2;
  // }
  // cout << s << endl;
  *add_time = time + dt;
  *add_dt = s*dt;

  // If the time stepping in next iteration is lesser than half of the current time step then the current step should be 
  // repeated
  // if (s<1)
  // {
  //   /* code */
  // }

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

//   double k[7][ndim], k_temp[ndim], yold[ndim];
//   // double k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], k7[ndim], ynew[ndim];

//   // eval_rhs(time,y,k1,flag_kappa);

//   for (int i = 0; i < 7; ++i)
//   {
//     for (int j = 0; j < ndim; ++j)
//     {
//       k[i][j] = 0;
//     }
//   }

//   for (int j = 0; j < 6; ++j)
//   {
//     for (int idim = 0; idim < ndim; ++idim)
//     {
//       temp[idim] = y[idim];
//       for (int i = 0; i < 6; ++i)
//       {
//         temp[idim] = temp[idim] + aij[j][i]*dt*k[i][idim];       
//       }
//     }

//     eval_rhs(time+dt*ci[j],temp,k_temp,flag_kappa,CurvSqr,);

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
//       yold[idim] = y[idim];
//       for (int j = 0; j < 6; ++j)
//       {
//           y[idim] = y[idim] + bi[j]*k[j][idim];
//       }
//   }
  
//   eval_rhs(time+dt,y,k_temp,flag_kappa,CurvSqr);

//   for (int idim = 0; idim < ndim; ++idim)
//   {
//       k[6][idim] = k_temp[idim];
//   }

//   for (int idim = 0; idim < ndim; ++idim)
//   {
//       temp[idim] = yold[idim];
//       for (int i = 0; i < 7; ++i)
//       {
//           temp[idim] = temp[idim] + bistar[i]*k[i][idim];
//       }

//       error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
//   }

//   error = sqrt(error/ndim);

//   s = 0.87*pow(tol_dt/error,0.20);
//   // cout << error << endl;
//   // if (s > 10)
//   // {
//   //     s = 10;
//   // }
//   // else if (s<0.2)
//   // {
//   //     s = 0.2;
//   // }

//   *add_time = time + dt;
//   *add_dt = s*dt;

//   // cout << *add_dt << endl;

// }


/* call eval_rhs(y,t,k1)
call eval_rhs(y+k1*(dt/2.),t+(dt/2.),k2)
call eval_rhs(y+k2*(dt/2.),t+(dt/2.),k3)
call eval_rhs(y+k3*dt,t+dt,k4)
y=y+dt*((k1/6.)+(k2/3.)+k3/3.+k4/6.) */
/* ----------------- */

