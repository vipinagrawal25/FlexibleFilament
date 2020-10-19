#include <iostream>
#include "model.h"
#include "MapDyn.h"
#include <cmath>
#include <sys/stat.h>
#include "constant.h"
#include <memory.h>
#include <cstdlib>
#include <algorithm>
#include "NewtonKrylov.h"
#include <string>
#include <fstream>
#include "misc.h"
#include "ode.h"
#include <time.h> 
/* -----------------------------------------------*/
using namespace std;
MV MM;
/* -----------------------------------------------*/
// void map_multiple_iter(double y[][ndim],double tAll);
void map_one_iter(double *y);
void write_map_param(string fname);
void rnkt4(double *y,double *add_time, double* add_dt) __attribute__((weak));
void rnkf45(double *y,double *add_time, double* add_dt) __attribute__((weak));
void eval_rhs(double *y) __attribute__((weak));
// Define coordinate transform in model.cpp file, if you wish to.
/*-----------------------------------------------*/
/*1) Prefix a means address to that particular varibale (avar -> address to var).
  2) Idea is not to use Eigen/Spectra namespace anywhere in these files.
  all of the converting function needs to be defined utils.
  3) wData,rData needs to be defined in IO.cpp file. Let's do it later.
  4) Define NewtonRaphson and NewtonKrylov method.*/
/*-----------------------------------------------*/
void assign_map_param(){
  // In principle, this function converts set of ODEs to a poincare map by taking poincare section.
  // For a non-automous flow such as the particular problem of Elastic string in periodically driven
  // Stokes flow, we take Poincare section sin(omega*t)=0 which is 3*N-1 dimensional hyper-plane, 
  // where 'N' is number of ODE to solve. 
  // This is also the other way of saying that for "n" iteration of the map, we integrate the ODEs for
  // "n" time period.
  /* Set up parameters for map iteration */
  MM.time = 0.;    // Ignore for discrete map
  MM.dt = 1.e-5;   // Ignore for discrete map
  MM.period = 1.;
  MM.iorbit = 1.;     // 0 if you already have the orbit, 
                      // 1 for calculating the orbit using Newton-Krylov
                      // 2 for letting the simulation evolve to a stable orbit.
  // 0 for no stability analysis, 1 for yes.
  // It computes the eigenvalues and save them in the folder.
  MM.istab = 0.;
  MM.irel_orb = 0.;     // Do you want to search for relative periodic orbits?
                        // 0 -> no, 1-> yes. Symmetry needs to be defined in model.cpp file.
  MM.mapdim = Np;
  // (*MM).iter_method=1;   // 1 -> Newton-Raphson
  // I am commenting things for diagnostics for the time being.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve crosses
  // Poincare section.
  // (*MM).tdiag = 0.;
  // (*MM).ldiag = 0.;
  // (*MM).ndiag = 1;
  write_map_param("map_intials.txt");
}
/* -----------------------------------------------*/
// Think of merging this function with write_param. 
// An array/std::vector containing the parameter names can be passed.
void write_map_param(string fname){
  ofstream pout(fname, ofstream::out);
  pout << "# =========== Map Parameters ==========\n";
  pout << "Period = " << MM.period << endl;
  pout << "iorbit = " << MM.iorbit << endl;
  pout << "istab = " << MM.istab << endl;
  pout << "irel_orb = "<< MM.irel_orb << endl; 
}
/* -----------------------------------------------*/
bool periodic_orbit(double y[], double fy[]){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method 
  // to go to nearby guess.
  int iorbit = MM.iorbit;
  int period = MM.period;
  //
  double step[ndim];
  int MaxTry = (int) MaxIter/period;
  int mapdim = MM.mapdim;
  //
  switch(iorbit){
    case 0:
      break ;    // the code should never come here.
    case 1:
    {
      memcpy(fy,y,mapdim*sizeof(double));
      map_multiple_iter(fy);
      if(SqEr(fy,y,mapdim)/norm(y,mapdim) < err_tol){
          cout << "Voila! you got the periodic orbit with period " << period << endl;
          cout << "I would go ahead and save it -:)" << endl;
      }
      else{
        bool success = newton_krylov(GG,y,fy,mapdim);
        add(fy,y,fy,mapdim);
        return success;
      }
      break;
    }
    case 2:
      memcpy(fy,y,mapdim*sizeof(double));
      for (int itry = 0; itry < MaxTry; ++itry){
        cout << "\n---------- Starting next try: " << itry << "----------" << endl;
        memcpy(y,fy,mapdim*sizeof(double));
        // Preparing fy for next iteration. 0th row is the starting point.
        map_multiple_iter(fy);
      }
      break;
  }
}
/* ----------------------------------------------- */
// void Jacobian(double DerM[][ndim], double x[]){
//   // Take a small step in every direction and calculate the differences.
//   // dy/dx = (f(x+dx) - f(x-dx))/2dx
//   int period = MM.period;
//   double ypos[ndim],yneg[ndim];
//   // MatrixXd DerM(ndim,ndim);
//   for (int idim = 0; idim < ndim; ++idim){
//     cout << "Starting calculation for row " << idim+1 << endl;
//     // Calculation for f(x+dx)
//     memcpy(ypos,x,ndim*sizeof(double));
//     ypos[idim] = x[idim]+delta;
//     map_multiple_iter(ypos);
//     // Calculation for f(x-dx)
//     memcpy(yneg,x,ndim*sizeof(double));
//     yneg[idim]=x[idim]-delta;
//     map_multiple_iter(yneg);
//     //Calculation of the derivative
//     for (int jdim = 0; jdim < ndim; ++jdim){
//       DerM[idim][jdim] = (ypos[jdim]-yneg[jdim])/(2*delta);
//     }
//   }
// }
/* ----------------------------------------------- */
// void map_multiple_iter(double fy[][ndim],double tAll[]){
//   double y[ndim];
//   MM.time=0;                                  // Set initial time to zero.
//   tAll[0] = MM.time;
//   int period = MM.period;
//   memcpy(y,fy[0],ndim*sizeof(double));
//   for (int iter = 0; iter < period; ++iter){
//     cout << "# -------- iteration: " << iter+1 << " started" << endl;
//     map_one_iter(&y[0]);
//     memcpy(fy[iter+1],y,ndim*sizeof(double));
//     tAll[iter] = MM.time;
//   }
//   cout << "# ------- Completed iterations. "<< endl;
// }
/* ----------------------------------------------- */
// G(x) = f(x) - x;
void GG(double y[]){
  int mapdim = MM.mapdim;
  double fy[mapdim];
  memcpy(fy,y,mapdim*sizeof(double));
  map_multiple_iter(fy);
  for (int idim = 0; idim < mapdim; ++idim){
    y[idim] = fy[idim]-y[idim];
  }
}
/* ----------------------------------------------- */
// it takes only one dimensional array as an input.
void map_multiple_iter(double y_trans[]){
  MM.time=0;
  int period = MM.period;
  double y[ndim];
  cout << "# Starting map iteration" << endl;
  // clock_t timer=clock();
  inv_coordinate_transform(y,y_trans);
  cout << "Where is the issue? "<< endl;
  for (int iter = 0; iter < period; ++iter){
    map_one_iter(&y[0]);
  }
  coordinate_transform(y_trans,y);
  // timer = clock() - timer;
  // double timeT = timer/CLOCKS_PER_SEC;
  // cout << "Time taken by function: " << timeT << "seconds" << endl;
}
/* ----------------------------------------------- */
void map_one_iter(double *y){
  if (SysType == "continuous"){
    // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
    // intermediate points as well.
    double Tmax = 1*2*M_PI/omega;   // We typically get this by defining Poincare section. 
                                    // which can depend on the initial condition, but not in the case
                                    // Elastic string.
    // The function for Poincare section should be defined in model.cpp file.
    double time = MM.time;
    double dt = MM.dt;
    double ldiag = 0;               //Theoretically it should be bool but both works.
    while(abs(time-Tmax)>time_tol){
      // Time-stepper should exactly stop at Tmax.
      if (Tmax<time+dt){
        dt=Tmax-time;             
      }
      if(TimeScheme == "rnkt4"){
        rnkt4(y, &time, &dt);
      }else if ( TimeScheme == "rnkf45" ){
        rnkf45(y, &time, &dt);
      }else{
        printf( "Algorithm\t%s\t not coded \n", TimeScheme);
        printf( "EXITING \n " );
        exit(1);
      }
      // cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<Tmax<<"\n";
    }
    MM.dt=dt;
    MM.time=time;
  }
  else if(SysType == "discrete"){
    eval_rhs(y);
  }
  else{
    cout << "ERROR: I don't know what is the system type." << endl;
    exit(1);
  }
}
/*----------------------------------------------- */
// function to find out whether given initial condition is a periodic orbit or not?
bool IsOrbit(double y[]){
  int mapdim = MM.mapdim; 
  double fy[mapdim];
  memcpy(fy,y,mapdim*sizeof(double));
  map_multiple_iter(fy);
  if (SqEr(fy,y,mapdim)/norm(y,mapdim)<err_tol){
    return 1;
  }
  else{
    return 0;
  }
}
/*----------------------------------------------- */
void coordinate_transform(double *y_trans, double *y){y_trans = y;}
void inv_coordinate_transform(double *y,double *y_trans){y=y_trans;}
/*----------------------------------------------- */