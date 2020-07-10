#include <iostream>
#include<fstream>
#include "ode.h"
#include "model.h"
#include "input.h"
#include<sys/stat.h>
unsigned int const ndim=Nensemble*pdim
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
  (*MM).tdiag = 0.;
  // I am commenting things for diagnostics for the time being.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve crosses 
  // Poincare section.
  (*MM).ldiag = 0.;
  (*MM).ndiag = 1;
}

void periodic_orbit(MV *MM,y0){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method to go to nearby guess.
}

int map_one_iter(){
  // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
  // intermediate points as welll.
  
}

int main(){
	// Here I will define whether I am calling the code for fixed point or periodic orbits.
	// Idea is to use the same code for periodic orbit and fixed point both.
	// Fixed point is a periodic orbit with time-period 1 for a map. 
	MV MM;
  double y0[ndim], vel[ndim];
  //Initialize them to zero.
  for (int idim = 0; idim < ndim; ++idim){
    y0[idim]=0;
    vel[idim]=0;
  }
  // First define all the parameters.
  pre_ode2map(&MM);          
  // Now get the initial configuration (t=0) of the system.
  iniconf(double *y0, double *vel);
  // ode2map(double *y, double *vel, MM);
  periodic_orbit(y0,MM);  // Function to determine whether given some conditions are periodic orbit 
                          //of any map? or Fixed point?
}