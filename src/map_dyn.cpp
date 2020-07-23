#include <iostream>
#include<fstream>
#include "ode.h"
#include "model.h"
#include "map_dyn.h"
#include "math.h"
#define time_tol 1.e-9
#define TimeScheme "rnkf45"
/* -----------------------------------------------*/
void map_one_iter(double *y, double *vel, MV* MM);
void pre_ode2map(MV* MM);
void periodic_orbit(double* y0,double* vel,MV* MM);
template<typename T>  // This type of function can take any variable type. Quite interesting, isn't it ?
T abs(T value);
/* -----------------------------------------------*/
int main(){
	// Here I will define whether I am calling the code for fixed point or periodic orbits.
	// Idea is to use the same code for periodic orbit and fixed point both.
	// Fixed point is a periodic orbit with time-period 1 for a map. 
	MV MM;
  double y0[ndim],vel[ndim];
  //Initialize them to zero.
  for (int idim = 0; idim < ndim; ++idim){
    y0[idim]=0;
    vel[idim]=0;
  }
  // First define all the parameters.
  pre_ode2map(&MM);          
  // Now get the initial configuration (t=0) of the system.
  iniconf(y0);
  // ode2map(double *y, double *vel, MM);
  // map_one_iter(&y0[0],&vel[0],&MM);
  periodic_orbit(&y0[0],&vel[0],&MM);
  // periodic_orbit(y0,&MM);  // Function to determine whether given some conditions are periodic orbit 
                          //of any map? or Fixed point?
}
/*-----------------------------------------------*/
void pre_ode2map(MV *MM){
  // In principle, this function converts set of ODEs to a poincare map by taking poincare section.
  // For a non-automous flow such as the particular problem of Elastic string in periodically driven
  // Stokes flow, we take Poincare section sin(omega*t)=0 which is 3*N-1 dimensional hyper-plane, 
  // where 'N' is number of ODE to solve. 
  // This is also the other way of saying that for "n" iteration of the map, we integrate the ODEs for
  // "n" time period.
  /* Set up parameters for map iteration */
  (*MM).time = 0.;
  (*MM).dt = 1.e-5;
  (*MM).niter = 1;
  // (*MM).iter_method=1;   // 1 -> Newton-Raphson

  // I am commenting things for diagnostics for the time being.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve crosses 
  // Poincare section.
  // (*MM).tdiag = 0.;
  // (*MM).ldiag = 0.;
  // (*MM).ndiag = 1;
}
/* -----------------------------------------------*/
void periodic_orbit(double* y0,double* vel,MV* MM){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method to go to nearby guess.
  double* y=y0;
  int niter=(*MM).niter;
  for (int iter = 0; iter < niter; ++iter){
    map_one_iter(&y[0],&vel[0],MM);
  }
}
/* -----------------------------------------------*/
void map_one_iter(double *y, double *vel, MV* MM){
  // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
  // intermediate points as welll.
  double Tmax = 1*2*M_PI/omega;  // We typically should get this by defining Poincare section. 
                                 // which can depend on the initial condition, but not in this case.
  // The function for Poincare section should be defined in model.cpp file.
  double time = (*MM).time;
  double dt = (*MM).dt;
  double ldiag = 0;               //Theoretically it should be bool but both works.
  while(abs(time-Tmax)>time_tol){
    if( TimeScheme == "rnkt4"){
      rnkt4(&y[0], &vel[0], &time, &dt);
    }else if ( TimeScheme == "rnkf45" ){
      rnkf45(&y[0], &vel[0], &time, &dt);
    }else{
      printf( " algorithm\t%s\t not coded \n", TimeScheme);
      printf( "EXITING \n " );
      exit(1);
    }
    cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<Tmax<<"\n";
  }
  (*MM).dt=dt;
  (*MM).time=time;
}
/*-----------------------------------------------*/
template<typename T>
T abs(T value){
  if (value>=0){
    return value;
  }else{
    return value*(-1);
  }
}