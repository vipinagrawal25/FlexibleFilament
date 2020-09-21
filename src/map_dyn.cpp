#include <iostream>
#include <fstream>
#include "ode.h"
#include "model.h"
#include "map_dyn.h"
#include "math.h"
#include <Eigen/Eigenvalues>
#include <sys/stat.h>
#include "constant.h"
#include <memory.h>
#include <cstdlib>
#include <algorithm>
/* -----------------------------------------------*/
using namespace std;
using namespace Eigen;
/* -----------------------------------------------*/
void assign_map_param(MV* MM);
void map_multiple_iter(double y[][ndim],MV *MM);
double* map_multiple_iter(double y[], MV *MM);
void map_one_iter(double *y, double *vel, MV* MM);
void write_map_param(MV *MM, string fname);
/*-----------------------------------------------*/
/*1) Prefix a means address to that particular varibale (avar -> address to var).
  2) Idea is not to use Eigen/Spectra namespace anywhere in these files.
  all of the converting function needs to be defined utils.
  3) wData,rData needs to be defined in IO.cpp file. Let's do it later.
  4) Define NewtonRaphson and NewtonKrylov method.*/
/*-----------------------------------------------*/
void assign_map_param(MV *MM){
  // In principle, this function converts set of ODEs to a poincare map by taking poincare section.
  // For a non-automous flow such as the particular problem of Elastic string in periodically driven
  // Stokes flow, we take Poincare section sin(omega*t)=0 which is 3*N-1 dimensional hyper-plane, 
  // where 'N' is number of ODE to solve. 
  // This is also the other way of saying that for "n" iteration of the map, we integrate the ODEs for
  // "n" time period.
  /* Set up parameters for map iteration */
  (*MM).time = 0.;    // Ignore for discrete map
  (*MM).dt = 1.e-5;   // Ignore for discrete map
  (*MM).period = 1.;
  (*MM).iorbit = 0.;  // 0 if you already have the orbit, 
                      // 1 for calculating the orbit using Newton-Raphson
                      // 2 for letting the simulation evolve to a stable orbit.
  // 0 for no stability analysis, 1 for yes.
  // It computes the eigenvalues and save them in the folder.
  (*MM).istab = 1.;
  (*MM).irel_orb = 0.;  // Do you want to search for relative periodic orbits?
                        // 0 -> no, 1-> yes. Symmetry needs to be defined in model.cpp file.
  // (*MM).iter_method=1;   // 1 -> Newton-Raphson
  // I am commenting things for diagnostics for the time being.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve crosses
  // Poincare section.
  // (*MM).tdiag = 0.;
  // (*MM).ldiag = 0.;
  // (*MM).ndiag = 1;
  write_map_param(MM,"map_intials.txt");
}
/* -----------------------------------------------*/
// Think of merging this function with write_param. 
// An array/std::vector containing the parameter names can be passed.
void write_map_param(MV *MM, string fname){s
  ofstream pout(fname, ofstream::out);
  pout << "# =========== Map Parameters ==========\n";
  pout << "Period = " << MM -> period << endl;
  pout << "iorbit = " << MM -> iorbit << endl;
  pout << "istab = " << MM -> istab << endl;
  pout << "irel_orb = "<< MM -> irel_orb << endl; 
}
/* -----------------------------------------------*/
void periodic_orbit(double y[],MV* aMM){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method 
  // to go to nearby guess.
  int iorbit = (*aMM).iorbit;
  int period = (*aMM).period;
  //
  double fy[period+1][ndim], fy_old[period+1][ndim];
  double fyvec[ndim],step[ndim];
  double GradF[ndim,ndim];
  int MaxTry = (int) MaxIter/period;
  //
  switch(iorbit){
    case 0:
      break ;    // the code should never come here.
    case 1:
      // newton_raphson(A function pointer which returns f(x), function to write gradient/gradient matrix, initial guess, MaxNewtonStep);
      // for (int itry = 0; itry < MaxTry; ++itry){
      //   memcpy(fy[0],y,ndim*sizeof(double));      // Filling first row of fy[][] by initial guess.
      //   // The multiple iter function will take care of one iteration running again and again.
      //   // and saving the intermediate data as well. It is like overdo but needed,
      //   // so that we do not commit silly mistakes of not initializing time to zero etc.
      //   map_multiple_iter(fy,vel_all,aMM);
      //   // write data iff we found an orbit
      //   if(SqEr(fy[period+1],fy[0],ndim)<err_tol){
      //     cout << "Voila! you got the periodic orbit with period " << period << endl;
      //     wData(fy,vel_all,period+1);                  // Code it in your model.cpp file
      //     return;
      //   }else{
      //     //get a new guess now and call periodic_orbit function again.
      //     MatrixXd GradF(ndim,ndim);
      //     cout << "\n---- Jacobian Calculation for newton-raphson iteration: " << itry 
      //          << "----" << endl;
      //     GradF = Jacobian(y,vel,aMM);
      //     cout << "------------- Jacobian Calculation Completed ------------- " 
      //          << endl << endl;
      //     // Convert fy[iter+1][] to vecXd to use Eigen.
      //     fyvec = arrtoVecXd(fy[period],ndim);
      //     // Take a new guess.
      //     step = GradF.inverse() * fyvec;
      //     // Now take step: y = y - step;
      //     transform(y, y + ndim, VecXdtoArr(&step,ndim), y, minus<double>());
      //   }
      // }
      break;
    case 2:
      // Implement NewtonKrylov method
      break;
    case 3:
      memcpy(fy[0],y,ndim*sizeof(double));      // filling first row of fy[][] by initial guess.
      map_multiple_iter(fy,aMM);
            
      if (SqEr(fy[0],fy[period],ndim) < err_tol){
        cout << "Voila!!! You got the periodic orbit with period " << period << endl;
        wData(fy,period+1);                  // Code it in your model.cpp file
        return;
      }
      //
      for (int itry = 0; itry < MaxTry; ++itry){
        cout << "\n---------- Starting next try: " << itry << "----------" << endl;
        memcpy(fy_old,fy,ndim*(period+1)*sizeof(double));
        // Preparing fy for next iteration. 0th row is the starting point.
        memcpy(fy[0],fy[period],ndim*sizeof(double));
        map_multiple_iter(fy,aMM);
        // Now take the difference of every row of fy_old and fy. Stop the simulation
        // if you find any periodic orbit.
        for (int iter = 1; iter < period+1; ++iter){
          cout << "Mean square distance: " << SqEr(fy_old[iter],fy[iter],ndim) << endl;
          if (SqEr(fy_old[iter],fy[iter],ndim)<err_tol){
            cout << "Voila!!! You got the periodic orbit with period " << period << endl;
            wData(fy,period+1);                  // Code it in your model.cpp file 
            return;
          }
        }
      }
      break;
  }
}
/* ----------------------------------------------- */
void Jacobian(double DerM[][ndim],double x[], double vel[], MV *aMM){
  // Take a small step in every direction and calculate the differences.
  // dy/dx = (f(x+dx) - f(x-dx))/2dx
  int period = (*aMM).period;
  double yp[ndim],yn[ndim],fyp[ndim],fyn[ndim];
  // MatrixXd DerM(ndim,ndim);
  for (int idim = 0; idim < ndim; ++idim){
    cout << "Starting calculation for row " << idim+1 << endl;
    // Calculation for f(x+dx)
    memcpy(yp,x,ndim*sizeof(double));
    yp[idim] = x[idim]+delta;
    fyp = map_multiple_iter(yp,aMM);
    // Calculation for f(x-dx)
    memcpy(yn,x,ndim*sizeof(double));
    yn[idim]=x[idim]-delta;
    fyn = map_multiple_iter(yn,aMM);
    //Calculation of the derivative
    for (int jdim = 0; jdim < ndim; ++jdim){
      DerM[jdim,idim] = (fyp[jdim]-fyn[jdim])/(2*delta);
    }
  }
}
/* ----------------------------------------------- */
// Wrapper function for map_multiple_iter.
// it takes only one dimensional array as an input.
double *map_multiple_iter(double y[],MV *aMM){
  double y[ndim];
  aMM->time=0;
  int period = (*aMM).period;
  cout << "Starting map iteration" << endl;
  for (int iter = 0; iter < period; ++iter){
    map_one_iter(&y[0],aMM); 
  }
  return y;
}
/* ----------------------------------------------- */
void map_multiple_iter(double fy[][ndim],MV *aMM){
  double y[ndim],vel[ndim];
  aMM->time=0;                                  // Set initial time to zero.
  int period = aMM->period;
  memcpy(y,fy[0],ndim*sizeof(double));
  for (int iter = 0; iter < period; ++iter){
    cout << "# -------- iteration: " << iter+1 << " started" << endl;
    map_one_iter(&y[0],aMM);
    memcpy(fy[iter+1],y,ndim*sizeof(double));
  }
  cout << "# ------- Completed iterations. "<< endl;
}
/* ----------------------------------------------- */
void map_one_iter(double *y, MV* MM){
  // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
  // intermediate points as well.
  double Tmax = 1*2*M_PI/omega;  // We typically get this by defining Poincare section. 
                                 // which can depend on the initial condition, but not in the case
                                // Elastic string.
  // The function for Poincare section should be defined in model.cpp file.
  double time = MM->time;
  double dt = MM->dt;
  double ldiag = 0;               //Theoretically it should be bool but both works.
  while(abs(time-Tmax)>time_tol){
    // Time-stepper should exactly stop at Tmax.
    if (Tmax<time+dt){
      dt=Tmax-time;             
    }
    if( TimeScheme == "rnkt4"){
      rnkt4(&y[0], &time, &dt);
    }else if ( TimeScheme == "rnkf45" ){
      rnkf45(&y[0], &time, &dt);
    }else{
      printf( "Algorithm\t%s\t not coded \n", TimeScheme);
      printf( "EXITING \n " );
      exit(1);
    }
    // cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<Tmax<<"\n";
  }
  MM->dt=dt;
  MM->time=time;
}
/* ----------------------------------------------- */
// function to find out whether given initial condition is a periodic orbit or not?
bool IsOrbit(double y[],MV *aMM){
  double fy[ndim];
  fy = map_multiple_iter(y,vel,aMM);
  if (SqEr(fy,y,ndim)<err_tol){
    return 1;
  }
  else{
    return 0;
  }
}