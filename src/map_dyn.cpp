#include <iostream>
#include <fstream>
#include "ode.h"
#include "model.h"
#include "map_dyn.h"
#include "math.h"
#include <Eigen/Eigenvalues>
#include<sys/stat.h>
/* -----------------------------------------------*/
using namespace std;
using namespace Eigen;
/* -----------------------------------------------*/
#define TimeScheme "rnkf45"
#define time_tol 1.e-9
#define err_tol 1.e-9
#define delta 1.e-4
#define wDataMeth 2
/* -----------------------------------------------*/
void map_one_iter(double *y, double *vel, MV* MM);
void pre_ode2map(MV* MM);
void periodic_orbit(double y0[],double vel[], MV* MM, int fnum);
template<typename T>  // This type of function definition can take any variable type. 
T abs(T value);       // Quite interesting, isn't it ?
void wData(ofstream *fptr, ofstream *fptr_vel, double y[], double vel[]) __attribute__((weak));
void map_multiple_iter(double y[],double vel[], MV *MM);
void map_multiple_iter(double y[],double vel[], MV *MM, ofstream *outfile, ofstream *outfile_vel);
double SqEr(double Arr1[], double Arr2[], int nn);
MatrixXd Jacobian(double x[], double vel[],MV *MM);
int lastfilenum();
/*-----------------------------------------------*/
/*1) Prefix a means address to that particular varibale (avar -> address to var).
*/ 
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
  (*MM).period = 1;
  (*MM).iorbit = 1; // 0 if you already have the orbit, 1 for calculating the orbit.
  // 0 for no stability analysis, 1 for yes.
  // It computes the eigenvalues and save them in the folder.
  (*MM).istab = 1;  
  // (*MM).iter_method=1;   // 1 -> Newton-Raphson
  // I am commenting things for diagnostics for the time being.
  // Since I am converting ODE to a map, I will just save things whenever the dynamical curve crosses 
  // Poincare section.
  // (*MM).tdiag = 0.;
  // (*MM).ldiag = 0.;
  // (*MM).ndiag = 1;
}
/* -----------------------------------------------*/
void periodic_orbit(double y0[],double vel[], MV* aMM, int fnum){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method 
  // to go to nearby guess.
  double *y=y0;
  double ynew[ndim];
  int period=(*aMM).period;
  
  ofstream outfile;
  ofstream outfile_vel;

  outfile.open("output/PSI0.txt");
  outfile_vel.open("output/VEL0.txt");
  wData(&outfile,&outfile_vel,y0,vel);                   // Code it in your model.cpp file

  // The multiple iter function will take care of one iteration running again and again.
  // and saving the intermediate data as well. It is like overdo but needed,
  // so that we do not commit silly mistakes of not initializing time to zero etc.
  map_multiple_iter(y,vel,aMM,&outfile,&outfile_vel);
  // Now check, did we hit the periodic orbit? If not, give a new guess.
  if(SqEr(y0,y,ndim)<err_tol){
    cout << "Voila! you got the periodic orbit with period " << period << endl;
    // Just write the output.
    wData(&outfile,&outfile_vel,y,vel);
    outfile.close();
    outfile_vel.close();
  }
  else{
    outfile.close();
    outfile_vel.close();
    // Somehow get a new guess now and call periodic_orbit function again.
    // First step is to get the derivative of the map.
    // dy/dx = (f(x+dx) - f(x-dx))/2dx
    MatrixXd DerM(ndim,ndim);
    DerM = Jacobian(y0,vel,aMM);
  }
}
/* ----------------------------------------------- */
MatrixXd Jacobian(double x[], double vel[], MV *aMM){
  // Take a small step in every direction and calculate the differences.
  // dy/dx = (f(x+dx) - f(x-dx))/2dx
  int period = (*aMM).period;
  MatrixXd DerM(ndim,ndim);
  for (int idim = 0; idim < ndim; ++idim){
    // Calculation for f(x+dx)
    double *yp=x;
    yp[idim]=x[idim]+delta;
    map_multiple_iter(yp,vel,aMM);
    // Calculation for f(x-dx)
    double *yn = x;
    yn[idim]=x[idim]-delta;
    map_multiple_iter(yn,vel,aMM);
    //Calculation of the derivative
    for (int jdim = 0; jdim < ndim; ++jdim){
      DerM(jdim,idim) = (yp[jdim]-yn[jdim])/(2*delta);
    }
  }
  return DerM;
}
/* ----------------------------------------------- */
// Wrapper function for map_multiple_iter.
// In case, we do not want to write data.
void map_multiple_iter(double y[], double vel[], MV *aMM){  
  (*aMM).time=0;
  int period = (*aMM).period;
  for (int iter = 0; iter < period; ++iter){
    cout << "Iteration " << iter << " started" << endl;
    map_one_iter(&y[0],&vel[0],aMM); 
  }
}
/* ----------------------------------------------- */
void map_multiple_iter(double y[],double vel[], MV *aMM, ofstream *aoutfile, ofstream *aoutfile_vel){
  (*aMM).time=0;                                  // Set initial time to zero.
  int period = (*aMM).period;
  for (int iter = 0; iter < period; ++iter){
    cout << "# -------- iteration: " << iter << " started" << endl;
    map_one_iter(&y[0],&vel[0],aMM);
    wData(aoutfile,aoutfile_vel,y,vel);
  }
  cout << "# ------- Completed iterations.---------- "<< endl;
}
/* ----------------------------------------------- */
void map_one_iter(double *y, double *vel, MV* MM){
  // This function convert ODE to map for 1 iteration. It also has flexibility to save a few 
  // intermediate points as well.
  double Tmax = 1*2*M_PI/omega;  // We typically get this by defining Poincare section. 
                                 // which can depend on the initial condition, but not in this case.
  // The function for Poincare section should be defined in model.cpp file.
  double time = (*MM).time;
  double dt = (*MM).dt;
  double ldiag = 0;               //Theoretically it should be bool but both works.
  while(abs(time-Tmax)>time_tol){
    // Time-stepper should exactly stop at Tmax.
    if (Tmax<time+dt){
      dt=Tmax-time;             
    }
    if( TimeScheme == "rnkt4"){
      rnkt4(&y[0], &vel[0], &time, &dt);
    }else if ( TimeScheme == "rnkf45" ){
      rnkf45(&y[0], &vel[0], &time, &dt);
    }else{
      printf( "Algorithm\t%s\t not coded \n", TimeScheme);
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
/*-----------------------------------------------*/
int lastfilenum(){
  // Returns the last file number
  int fnum=1;
  string filename = "output/PSI1";
  struct stat buffer;
  while(!stat(filename.c_str(), &buffer)){
    fnum=fnum+1;
    string filename = "output/PSI";
    filename.append(to_string(fnum));
  }
  return fnum;
}
/*-----------------------------------------------*/
// Move it to utilities
void wData(ofstream *fptr, ofstream *fptr_vel, double y[], double vel[]){
  switch(wDataMeth){
    case 1:
      for(int ip = 0; ip < Np; ++ip){
        for (int jp = 0; jp < pp; ++jp){
          *fptr << y[pp*ip+jp];
          *fptr << '\t';
        }
        // Now just throw away next three numbers as they contain values of velocity.
        for (int jp = 0; jp < pp; ++jp){
          *fptr << vel[pp*ip+jp];
          *fptr << '\t';
        }
        *fptr << '\n';
      }
      break;

    case 2:
      for (int idim = 0; idim < ndim; ++idim){
        *fptr << y[idim] << "\t";
        *fptr_vel << vel[idim] << "\t";
      }
      *fptr <<  "\n";
      *fptr << "#------------------------------------------------#\n";

      *fptr_vel <<  "\n";
      *fptr_vel << "#------------------------------------------------#\n";
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
// Move it to utilities
void rData(ofstream *fptr, ofstream *fptr_vel, double y[], double vel[]){
  switch(wDataMeth){
    case 1:
      for(int ip = 0; ip < Np; ++ip){
        for (int jp = 0; jp < pp; ++jp){
          *fptr << y[pp*ip+jp];
          *fptr << '\t';
        }
        // Now just throw away next three numbers as they contain values of velocity.
        for (int jp = 0; jp < pp; ++jp){
          *fptr << vel[pp*ip+jp];
          *fptr << '\t';
        }
        *fptr << '\n';
      }
      break;

    case 2:
      for (int idim = 0; idim < ndim; ++idim){
        *fptr << y[idim] << "\t";
        *fptr_vel << vel[idim] << "\t";
      }
      *fptr <<  "\n";
      *fptr << "#------------------------------------------------#\n";

      *fptr_vel <<  "\n";
      *fptr_vel << "#------------------------------------------------#\n";
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
// Move it to utilities
double SqEr(double Arr1[], double Arr2[], int nn){
  double error;
  for (int ii = 0; ii < nn; ++ii){
    error=error+(Arr2[ii]-Arr1[ii])*(Arr2[ii]-Arr1[ii]);
  }
  error = sqrt(error)/ndim;
  return error;
}
/* -----------------------------------------------*/
int main(){
  // Here I will define whether I am calling the code for fixed point or periodic orbits.
  // Idea is to use the same code for periodic orbit and fixed point both.
  // Fixed point is a periodic orbit with time-period 1 for a map. 
  MV MM;
  double y0[ndim],vel[ndim]; 
  MatrixXd DerM(ndim,ndim);               // This is better to allocate matrix.
  int fnum=0;
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
  fnum = lastfilenum();
  if (MM.iorbit){
    // The function use Newton-Raphson and calculate the nearest periodic orbit.
    periodic_orbit(&y0[0],&vel[0],&MM,fnum);
  }
  if(MM.istab){
    double *y=y0;
    // First check whether you actually have the periodic orbit?
    cout << "Is it a periodic orbit?"
            " (if you don't want this, please comment it out in map_dyn.cpp) " << endl;
    map_multiple_iter(y,vel,&MM);
    if(SqEr(y0,y,ndim)<err_tol){
      cout << "Yes!!! It is a periodic orbit. I shall calculate stability now." << endl;
      DerM = Jacobian(y0,vel,&MM);
      VectorXcd eivals = DerM.eigenvalues();
      cout << "The eigenvalues  are: " << endl << eivals << endl;
    }else{
      cout << "The guess is not periodic orbit." << endl << 
      " Did you cross-check the data or time-period? " << endl <<
      "If yes, use our periodic orbit solver to get nearest orbit to the guess." 
      "Set istab=0, iorbit=1" << endl;
    }
  }
}