#include <iostream>
#include "model.h"
#include "mapdyn.h"
#include <cmath>
#include <sys/stat.h>
#include "constant.h"
#include <memory.h>
#include <cstdlib>
#include <algorithm>
#include "newton-krylov.h"
#include <string>
#include <fstream>
#include "misc.h"
#include "ode.h"
#include <time.h>
/* -----------------------------------------------*/
using namespace std;
MV MM;                        // Global variable MM
/* -----------------------------------------------*/
// void map_multiple_iter(double y[][ndim],double tAll);
void write_map_param(string fname);
void rnkt4(double *y,double *add_time, double* add_dt) __attribute__((weak));
void rnkf45(double *y,double *add_time, double* add_dt) __attribute__((weak));
void eval_rhs(double *y) __attribute__((weak));
void pre_next_iter(double *y,double *ytrans) __attribute__((weak));
// void map_trans_one_iter(double ytrans[]);
// Define coordinate transform in model.cpp file, if you wish to.
/*-----------------------------------------------*/
/*1) Prefix a means address to that particular varibale (avar -> address to var).
  2) Idea is not to use Eigen/Spectra namespace anywhere in these files.
  all of the converting function needs to be defined in utils.
/*-----------------------------------------------*/
void assign_map_param(){
  // In principle, this function converts set of ODEs to a poincare map by taking poincare section.
  // For a non-automous flow such as the particular problem of Elastic string in periodically driven
  // Stokes flow, we take Poincare section sin(omega*t)=0 which is 3*N-1 dimensional hyper-plane, 
  // where 'N' is number of ODE to solve. 
  // This is also the other way of saying that for "n" iteration of the map, we integrate the ODEs 
  // for "n" time period.
  /* Set up parameters for map iteration */
  MM.time = 0.;    // Ignore for discrete map
  MM.dt = 1.e-4;   // Ignore for discrete map
  MM.period = 4.;
  MM.iorbit = 1.;     // 0 if you already have the orbit,
                      // 1 for calculating the orbit using Newton-Krylov
                      // 2 for letting the simulation evolve to a stable orbit.
  // 0 for no stability analysis, 1 for yes.
  // It computes the eigenvalues and save them in a file.
  MM.istab = 0.;
  // MM.irel_orb = 0.;     // Do you want to search for relative periodic orbits?
                           // 0 -> no, 1-> yes. Symmetry needs to be defined in model.cpp file.
  MM.mapdim = Np;
  MM.guess_space = "Real";  // take two values "Real" or "Transformed"
  // (*MM).iter_method=1;   // 1 -> Newton-Raphson
  // I am commenting things for diagnostics.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve 
  // crosses Poincare section.
  // (*MM).tdiag = 0.;
  // (*MM).ldiag = 0.;
  // (*MM).ndiag = 1.;
  int suffix = 1;
  // string filename = "map_intials";
  // string filename_temp = filename+".txt";
  // while(IsPathExist(filename_temp)){
  //   suffix++;
  //   filename_temp = filename + to_string(suffix) + ".txt";
  // }
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
  // pout << "irel_orb = "<< MM.irel_orb << endl;
  pout << "mapdim = "<< MM.mapdim << endl;
  pout << "Guess space = " << MM.guess_space << endl;
}
/* -----------------------------------------------*/
void get_yall(double ytrans_all[], double yall[], double time[]){
  int mapdim = MM.mapdim;
  int period = MM.period;
  double fytrans[mapdim],fy[ndim];
  memcpy(fy,yall,ndim*sizeof(double));
  MM.time=0;
  cout << "# Starting map iteration " << endl;
  for (int iter = 0; iter < period-1; ++iter){
    map_one_iter(fy);
    coordinate_transform(fytrans,fy);
    memcpy(ytrans_all+mapdim*(iter+1),fytrans,mapdim*sizeof(double));
    memcpy(yall+ndim*(iter+1),fy,ndim*sizeof(double));
    time[iter+1]=MM.time;
  }
}
/*-----------------------------------------------*/
void get_yall(double ytrans_all[], double fytrans[], 
                    double yall[], double fy[], double time[]){
  int mapdim=MM.mapdim;
  int period=MM.period;
  get_yall(ytrans_all,yall,time);
  memcpy(fy,yall+ndim*(period-1),ndim*sizeof(double));
  map_one_iter(fy);
  coordinate_transform(fytrans,fy);
  time[period]=MM.time;
}
/*-----------------------------------------------*/
bool periodic_orbit(double ytrans_all[], double fytrans[],
                    double yall[], double fy[], double time[]){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method 
  // to go to nearby guess.
  int iorbit = MM.iorbit;
  int period = MM.period;
  int mapdim = MM.mapdim;
  //
  double ytrans[mapdim],y[ndim];
  memcpy(ytrans,ytrans_all,mapdim*sizeof(double));
  int MaxTry = (int) MaxIter/period;
  // int MaxTry = 1;
  bool success = 0;
  double Error = 0;
  //
  switch(iorbit){
    case 0:
      break ;         // The code should never come here.
    case 1:
      get_yall(ytrans_all,fytrans,yall,fy,time);
      if(norm(fytrans,mapdim)<err_tol*mapdim*10){Error = SqEr(fytrans,ytrans_all,mapdim)/(err_tol*mapdim);}
      else{Error = SqEr(fytrans,ytrans_all,mapdim)/norm(fytrans,mapdim);}
      // cout << "# Error = " << SqEr(fytrans,ytrans_all,mapdim)/mapdim << endl;
      if(Error < err_tol){
      // if(1){
        success = 1;
        return success;
      }else{
        bool success = newton_krylov(GG,ytrans,fytrans,mapdim,err_tol);
        add(fytrans,ytrans,fytrans,mapdim);
        inv_coordinate_transform(fy,fytrans);
        time[period]=MM.time;
        //
        inv_coordinate_transform(y,ytrans);
        memcpy(ytrans_all,ytrans,mapdim*sizeof(double));
        memcpy(yall,y,ndim*sizeof(double));
        get_yall(ytrans_all,yall,time);
        return success;
      }
      break;
    case 2:
      for (int itry = 0; itry < MaxTry; ++itry){
        cout << "# Starting next try:" << itry << endl;
        get_yall(ytrans_all,fytrans,yall,fy,time);
        // print(yall,period*ndim);
        if(norm(fytrans,mapdim)<err_tol*mapdim*10){Error = SqEr(fytrans,ytrans_all,mapdim)/(err_tol*mapdim);}
        else{Error = SqEr(fytrans,ytrans_all,mapdim)/norm(fytrans,mapdim);}
        cout << "# Error= " << Error << endl;
        if(Error < err_tol){
          // print(fytrans,mapdim);
          // print(ytrans_all,ndim);
          success = 1;
          return success;
        }
        // Preparing for next iteration
        memcpy(yall,fy,ndim*sizeof(double));
        memcpy(ytrans_all,fytrans,mapdim*sizeof(double));
        // print(yall,ndim*period);
        // An additional function in case you want to change anything in model before starting the next iteration.
        // A weak function is defined at the bottom of this file.
        pre_next_iter(fy,fytrans);
      }
      break;
  }
  return success;
}
/* ----------------------------------------------- */
// void Jacobian(double DerM[][MM.mapdim], double x[]){
//   // Take a small step in every direction and calculate the differences.
//   // dy/dx = (f(x+dx) - f(x-dx))/2dx
//   int period = MM.period;
//   double ypos[ndim],yneg[ndim];
//   // MatrixXd DerM(ndim,ndim);
//   for (int idim = 0; idim < ndim; ++idim){
//     cout << "Starting calculation for row " << idim+1 << endl;
//   // Calculation for f(x+dx)
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
void Jacobian(Matd *DerM, double x[]){
  // Take a small step in every direction and calculate the differences.
  // dy/dx = (f(x+dx) - f(x-dx))/2dx
  double delta = 1.e-2;
  cout << "# delta = " << delta <<endl;
  int period = MM.period;
  int mapdim = MM.mapdim;
  //
  double ypos[mapdim],yneg[mapdim];
  // MatrixXd DerM(mapdim,mapdim);
  for (int idim = 0; idim < mapdim; ++idim){
    cout << "# Starting calculation for row " << idim << endl;
    // Calculation for f(x+dx)
    memcpy(ypos,x,mapdim*sizeof(double));
    ypos[idim] = x[idim]+delta;
    map_multiple_iter(ypos);
    // Calculation for f(x-dx)
    memcpy(yneg,x,mapdim*sizeof(double));
    yneg[idim]=x[idim]-delta;
    map_multiple_iter(yneg);
    // Calculation of the derivative
    for (int jdim = 0; jdim < mapdim; ++jdim){
      (*DerM)(idim,jdim) = (ypos[jdim]-yneg[jdim])/(2*delta);
    }
  }
  // for (int idim = 0; idim < mapdim-2; ++idim){
  //   cout << "# Starting calculation for row " << idim+1 << endl;
  // // Calculation for f(x+dx)
  //   memcpy(ypos,x,mapdim*sizeof(double));
  //   ypos[idim+1] = x[idim+1]+delta;
  //   map_multiple_iter(ypos);
  //   // Calculation for f(x-dx)
  //   memcpy(yneg,x,mapdim*sizeof(double));
  //   yneg[idim+1]=x[idim+1]-delta;
  //   map_multiple_iter(yneg);
  //   // Calculation of the derivative
  //   for (int jdim = 0; jdim < mapdim-2; ++jdim){
  //     (*DerM)(idim,jdim) = (ypos[jdim+1]-yneg[jdim+1]);
  //   }
  // }
  // (*DerM) = 1/(2*delta)*(*DerM);
}
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
void GG(double ytrans[]){
  MM.time = 0;
  int mapdim = MM.mapdim;
  int period = MM.period;
  double fytrans[mapdim],fy[ndim];
  inv_coordinate_transform(fy,ytrans);
  cout << "# Starting map iteration" << endl;
  for (int iter = 0; iter < period; ++iter){
    map_one_iter(fy);
  }
  coordinate_transform(fytrans,fy);
  pre_next_iter(fy,fytrans);
  for (int idim = 0; idim < mapdim; ++idim){
    ytrans[idim] = fytrans[idim]-ytrans[idim];
  }
}
/* ----------------------------------------------- */
// it takes only one dimensional array as an input.
void map_multiple_iter(double ytrans[]){
  // int mapdim = MM.mapdim;
  // double fytrans[mapdim];
  // memcpy(fytrans,ytrans,mapdim*sizeof(double));
  // GG(fytrans);
  // for (int idim = 0; idim < mapdim; ++idim){
  //   ytrans[idim] = fytrans[idim]+ytrans[idim];
  // }
  MM.time=0;
  int period = MM.period;
  double y[ndim];
  cout << "# Starting map iteration" << endl;
  inv_coordinate_transform(y,ytrans);
  for (int iter = 0; iter < period; ++iter){
    map_one_iter(y);
  }
  coordinate_transform(ytrans,y);
}
/* ----------------------------------------------- */
void map_one_iter(double *y){
  if (SysType == "continuous"){
    // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
    // intermediate points as well.
    double time = MM.time;
    double Tmax = period + time;    // We typically get this by defining Poincare section. 
                                    // which can depend on the initial condition, but not in the case
                                    // Elastic string.
    // The function for Poincare section should be defined in model.cpp file.
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
    // print(y,ndim);
  }
  else if(SysType == "discrete"){
    eval_rhs(y);
  }
  else{
    cout << "# ERROR: I don't know what is the system type." << endl;
    exit(1);
  }
}
/*-----------------------------------------------*/
// function to find out whether given initial condition is a periodic orbit or not?
bool IsOrbit(double y[]){
  int mapdim = MM.mapdim;
  double fy[mapdim];
  memcpy(fy,y,mapdim*sizeof(double));
  map_multiple_iter(fy);
  double Error=0;
  if(norm(fy,mapdim)<err_tol*mapdim*10){Error = SqEr(fy,y,mapdim)/(err_tol*mapdim);}
  else{Error = SqEr(fy,y,mapdim)/norm(fy,mapdim);}
  cout << "# Error= " << Error << endl;
  return 1;
  if (Error<err_tol){
    return 1;
  }
  else{
    return 0;
  }
}
/*----------------------------------------------- */
void __attribute__((weak)) coordinate_transform(double *ytrans, double *y){
  int mapdim = MM.mapdim;
  memcpy(ytrans,y,mapdim*sizeof(double));
}
/*----------------------------------------------- */
void __attribute__((weak)) inv_coordinate_transform(double *y,double *ytrans){
  memcpy(y,ytrans,ndim*sizeof(double));
}
/*----------------------------------------------- */
void pre_next_iter(double *y,double *ytrans){}
/*----------------------------------------------- */