#include <iostream>
#include<fstream>
#include<cmath>
#include "ode.h"
using namespace std;
/*********************************/
void euler(unsigned int ndim, double *y, double time,double dt){
  double  k1[ndim];
  int idim;
  eval_rhs(time,y,k1);
  for(idim=0;idim<ndim;idim++){
    y[idim]=y[idim]+k1[idim]*dt;  
  }
}
/*********************************/
void rnkt2(unsigned int ndim, double *y, double time,double dt){
  double  temp[ndim],k1[ndim],k2[ndim];
  int idim;
//
  eval_rhs(time,y,k1);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k1[idim]*dt/2.;
  }
  eval_rhs(time+(dt/2.),temp,k2);
  for(idim=0;idim<ndim;idim++){
    y[idim]=y[idim]+k2[idim]*dt;
  }
}
/*********************************/
void rnkt4(unsigned int ndim, double *y, double time, double dt){
  double  temp[ndim],k1[ndim],k2[ndim],k3[ndim],k4[ndim];
  int idim;
  eval_rhs(time,y,k1);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k1[idim]*dt/2.;
  }
  eval_rhs(time+(dt/2.),temp,k2);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k2[idim]*dt/2.;
  }
  eval_rhs(time+(dt/2.),temp,k3);
  for(idim=0;idim<ndim;idim++){
    temp[idim]=y[idim]+k3[idim]*dt;
  }
  eval_rhs(time+dt,temp,k4);
  for(idim=0;idim<ndim;idim++){
    y[idim]=y[idim]+dt*(  (k1[idim]/6.) + (k2[idim]/3.) + (k3[idim]/3.) + (k4[idim]/6.) );
  }
}

void rnkf45(unsigned int ndim, double *y, double time, double* add_dt)
{
	// Details of method: http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf

	double temp[ndim], k1[ndim], k2[ndim], k3[ndim], k4[ndim], k5[ndim], k6[ndim], s;
	int idim ;
	double error = 0;
	double dt = *add_dt;
	double tol_dt = pow(10,-6)*dt;

	eval_rhs(time,y,k1);

	for(idim=0;idim<ndim;idim++)
	{
    temp[idim]=y[idim]+k1[idim]*dt/4.;
  }	
  eval_rhs(time+(dt/4.),temp,k2);
  	
  for(idim=0;idim<ndim;idim++)
  {
    temp[idim]=y[idim]+(3*k1[idim]+9*k2[idim])*dt/32.;
 	}
 	eval_rhs(time+(3/8.)*dt,temp,k3);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim]+ (1932*k1[idim]-7200*k2[idim]+7296*k3[idim])*dt/2197. ;
 	}
 	eval_rhs(time+(12/13.)*dt, temp, k4);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] + 439*dt*k1[idim]/216. - 8*k2[idim]*dt + 3680*k3[idim]*dt/513. - 845*k4[idim]*dt/4104. ;
 	}
 	eval_rhs(time+dt, temp, k5);

 	for (int idim = 0; idim < ndim; ++idim)
 	{
 		temp[idim] = y[idim] - 8*k1[idim]*dt/27. + 2*k2[idim]*dt - 3544*k3[idim]*dt/2565. + 1859*k4[idim]*dt/4104. - 11*k5[idim]*dt/40. ;
 	}
 	eval_rhs(time+dt/2., temp, k6);


 	for (int idim = 0; idim < ndim; ++idim)
 	{
    temp[idim] = y[idim] + 25*k1[idim]*dt/216. + 1408*k3[idim]*dt/2565. + 2197*k4[idim]*dt/4104 - k5[idim]*dt/5.;
 		// cout << temp[idim] << endl;
    y[idim] = y[idim]+ 16*k1[idim]*dt/135. + 6656*k3[idim]*dt/12825. + 28561*k4[idim]*dt/56430 - 9*k5[idim]*dt/50. + 2*k6[idim]*dt/55.;

 		error = error + (temp[idim]-y[idim])*(temp[idim]-y[idim]);
 	}
  error = sqrt(error);
  // cout << error << endl;
 	s = 0.84*pow(tol_dt/error,0.25);
  // cout << s << endl;
 	*add_dt = s*dt;
}	


/* call eval_rhs(y,t,k1)
call eval_rhs(y+k1*(dt/2.),t+(dt/2.),k2)
call eval_rhs(y+k2*(dt/2.),t+(dt/2.),k3)
call eval_rhs(y+k3*dt,t+dt,k4)
y=y+dt*((k1/6.)+(k2/3.)+k3/3.+k4/6.) */
/* ----------------- */

