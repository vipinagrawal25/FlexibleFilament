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
void map_one_iter(double *y, double *vel, MV* MM);
void pre_ode2map(MV* MM);
void periodic_orbit(double y0[],double vel[], MV* MM);
template<typename T>  // This type of function definition can take any variable type. 
T abs(T value);       // Quite interesting, isn't it ?
void wData(ofstream *fptr, ofstream *fptr_vel, double y[], double vel[]) __attribute__((weak));
void wData(double fy[][ndim], double vel_all[][ndim], int row) __attribute__((weak));
void map_multiple_iter(double y[],double vel[], MV *MM);
void map_multiple_iter(double y[][ndim],double vel[][ndim], MV *MM);
double SqEr(double Arr1[], double Arr2[],int ndim);
MatrixXd Jacobian(double x[], double vel[],MV *MM);
bool IsPathExist(const std::string &s)__attribute__((weak));
VectorXd arrtoVecXd(double arr[], int ndim);
double* VecXdtoArr(VectorXd* arrvec, int ndim);
bool IsOrbit(double y[],MV *aMM);
void write_map_param(MV *MM, string fname);
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
void periodic_orbit(double y[],double vel[], MV* aMM){
  // This function decide whether given initial condition is a periodic orbit or not.
  // If the initial point is not a periodic orbit, it uses the newton-raphson method 
  // to go to nearby guess.
  int iorbit = (*aMM).iorbit;
  int period = (*aMM).period;
  //
  double fy[period+1][ndim], vel_all[period+1][ndim],fy_old[period+1][ndim];
  VectorXd fyvec(ndim),step(ndim);
  MatrixXd GradF(ndim,ndim);
  int MaxTry = (int)MaxIter/period; 
  //
  switch(iorbit){
    case 0:
      break ;    // the code should never come here.
    case 1:
      for (int itry = 0; itry < MaxTry; ++itry){
        memcpy(fy[0],y,ndim*sizeof(double));      // Filling first row of fy[][] by initial guess.
        // The multiple iter function will take care of one iteration running again and again.
        // and saving the intermediate data as well. It is like overdo but needed,
        // so that we do not commit silly mistakes of not initializing time to zero etc.
        map_multiple_iter(fy,vel_all,aMM);
        // write data iff we found an orbit
        if(SqEr(fy[period+1],fy[0],ndim)<err_tol){
          cout << "Voila! you got the periodic orbit with period " << period << endl;
          wData(fy,vel_all,period+1);                  // Code it in your model.cpp file
          return;
        }else{
          //get a new guess now and call periodic_orbit function again.
          MatrixXd GradF(ndim,ndim);
          cout << "\n---- Jacobian Calculation for newton-raphson iteration: " << itry 
               << "----" << endl;
          GradF = Jacobian(y,vel,aMM);
          cout << "------------- Jacobian Calculation Completed ------------- " 
               << endl << endl;
          // Convert fy[iter+1][] to vecXd to use Eigen.
          fyvec = arrtoVecXd(fy[period],ndim);
          // Take a new guess.
          step = GradF.inverse() * fyvec;
          // Now take step: y = y - step;
          transform(y, y + ndim, VecXdtoArr(&step,ndim), y, minus<double>()); 
        }
      }
      break;
    case 2:
      memcpy(fy[0],y,ndim*sizeof(double));      // filling first row of fy[][] by initial guess.
      map_multiple_iter(fy,vel_all,aMM);      
      if (SqEr(fy[0],fy[period],ndim) < err_tol){
        cout << "Voila!!! You got the periodic orbit with period " << period << endl;
        wData(fy,vel_all,period+1);                  // Code it in your model.cpp file 
        return;
      }
      //
      for (int itry = 0; itry < MaxTry; ++itry){
        cout << "\n---------- Starting next try: " << itry << "----------" << endl;
        memcpy(fy_old,fy,ndim*(period+1)*sizeof(double));
        // Preparing fy for next iteration. 0th row is the starting point.
        memcpy(fy[0],fy[period],ndim*sizeof(double));
        map_multiple_iter(fy,vel_all,aMM);
        // Now take the difference of every row of fy_old and fy. Stop the simulation
        // if you find any periodic orbit.
        for (int iter = 1; iter < period+1; ++iter){

          cout << "Mean square distance: " << SqEr(fy_old[iter],fy[iter],ndim) << endl;
          if (SqEr(fy_old[iter],fy[iter],ndim)<err_tol){
            cout << "Voila!!! You got the periodic orbit with period " << period << endl;
            wData(fy,vel_all,period+1);                  // Code it in your model.cpp file 
            return;
          }
        }
      }
      break;
  }
}
/* ----------------------------------------------- */
MatrixXd Jacobian(double x[], double vel[], MV *aMM){
  // Take a small step in every direction and calculate the differences.
  // dy/dx = (f(x+dx) - f(x-dx))/2dx
  int period = (*aMM).period;
  double yp[ndim],yn[ndim];
  MatrixXd DerM(ndim,ndim);
  for (int idim = 0; idim < ndim; ++idim){
    cout << "Starting calculation for row " << idim+1 << endl;
    // Calculation for f(x+dx)
    memcpy(yp,x,ndim*sizeof(double));
    yp[idim] = x[idim]+delta;
    map_multiple_iter(yp,vel,aMM);
    // Calculation for f(x-dx)
    memcpy(yn,x,ndim*sizeof(double));
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
// it takes only one dimensional array as an input.
void map_multiple_iter(double y[], double vel[], MV *aMM){
  // Dignosis //
  (*aMM).time=0;
  int period = (*aMM).period;
  cout << "Starting map iteration" << endl;
  for (int iter = 0; iter < period; ++iter){
    map_one_iter(&y[0],&vel[0],aMM); 
  }
}
/* ----------------------------------------------- */
void map_multiple_iter(double fy[][ndim],double vel_all[][ndim], MV *aMM){
  double y[ndim],vel[ndim];
  aMM->time=0;                                  // Set initial time to zero.
  int period = (*aMM).period;
  memcpy(y,fy[0],ndim*sizeof(double));
  for (int iter = 0; iter < period; ++iter){
    cout << "# -------- iteration: " << iter+1 << " started" << endl;
    map_one_iter(&y[0],&vel[0],aMM);
    memcpy(fy[iter+1],y,ndim*sizeof(double));
    memcpy(vel_all[iter+1],vel,ndim*sizeof(double));
  }
  cout << "# ------- Completed iterations. "<< endl;
}
/* ----------------------------------------------- */
void map_one_iter(double *y, double *vel, MV* MM){
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
      rnkt4(&y[0], &vel[0], &time, &dt);
    }else if ( TimeScheme == "rnkf45" ){
      rnkf45(&y[0], &vel[0], &time, &dt);
    }else{
      printf( "Algorithm\t%s\t not coded \n", TimeScheme);
      printf( "EXITING \n " );
      exit(1);
    }
    // cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<Tmax<<"\n";
  }
  (*MM).dt=dt;
  (*MM).time=time;
}
/* ----------------------------------------------- */
// function to find out whether given initial condition is a periodic orbit or not?
bool IsOrbit(double y[],MV *aMM){
  double vel[ndim],fy[ndim];
  memcpy(fy,y,ndim*sizeof(double));
  map_multiple_iter(fy,vel,aMM);
  if (SqEr(fy,y,ndim)<err_tol){
    return 1;
  }
  else{
    return 0;
  }
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
// Move it to utilities
bool IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}
/*-----------------------------------------------*/
// Move it to utilities
// int lastfilenum(){
//   // Returns the last file number
//  =1;
//   string filename = "data/orbit1";
//   while(IsPathExist(filename)){
//     fnum=fnum+1;
//     string filename = "data/orbit" + to_string(fnum);
//   }
//   return fnum;
// }
/*-----------------------------------------------*/
// Move it to utilities
void wData(ofstream *fptr, ofstream *fptr_vel, double *y, double *vel, MV* MM){
  switch(wDataMeth){
    case 1:
      for(int ip = 0; ip < Np; ++ip){
        for (int jp = 0; jp < pp; ++jp){
          *fptr << y[pp*ip+jp];
          *fptr << '\t';
        }
        // Now next three numbers contain values of velocity.
        for (int jp = 0; jp < pp; ++jp){
          *fptr << vel[pp*ip+jp];
          *fptr << '\t';
        }
        *fptr << '\n';
      }
      break;

    case 2:
      (*fptr) << MM->time << "\t";
      (*fptr_vel) << MM->time << "\t";
      for (int idim = 0; idim < ndim; ++idim){
        *fptr << y[idim] << "\t";
        *fptr_vel << vel[idim] << "\t";
      }
      *fptr <<  endl;
      *fptr << "#------------------------------------------------#" << endl;

      *fptr_vel <<  endl;
      *fptr_vel << "#------------------------------------------------#" << endl;
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
// Wrapper function for wData, in case y and vel contains multiple rows.
// We do not pass the file pointer.
void wData(double y[][ndim], double vel[][ndim], MV* aMM,int nrow){
  switch(wDataMeth){
    case 1:{
      // for (int irow = 0; irow < nrow; ++irow){
      //   ofstream outfile("output/var"+to_string()+".txt",fstream::out);
      //   for(int ip = 0; ip < Np; ++ip){
      //     for (int jp = 0; jp < pp; ++jp){
      //       file << y[irow][pp*ip+jp];
      //       file << '\t';
      //     }
      //     // Now just throw away next three numbers as they contain values of velocity.
      //     for (int jp = 0; jp < pp; ++jp){
      //       file << vel[irow][pp*ip+jp];
      //       file << '\t';
      //     }
      //     file << '\n';
      //   }

      // }
      cout << "This is not necessary. I am EXITING from function wData." << endl;
      break;
    }
    case 2:{
      ofstream outfile("PSI",ofstream::out);
      ofstream outfile_vel("VEL",ofstream::out);
      for (int irow = 0; irow < nrow; ++irow){
        wData(&outfile,&outfile_vel,y[irow],vel[irow],aMM);
      }
      outfile.close();
      outfile_vel.close();
      break;
      // for (int idim = 0; idim < ndim; ++idim){
      //   *fptr << y[idim] << "\t";
      //   *fptr_vel << vel[idim] << "\t";
      // }
      // *fptr <<  "\n";
      // *fptr << "#------------------------------------------------#\n";

      // *fptr_vel <<  "\n";
      // *fptr_vel << "#------------------------------------------------#\n";
      // break;
    }
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
// Move it to utilities
double SqEr(double Arr1[], double Arr2[],int nn){
  double error=0;
  for (int cdim=0; cdim < nn; cdim++){
    error=error+(Arr2[cdim]-Arr1[cdim])*(Arr2[cdim]-Arr1[cdim]);
  }
  error = sqrt(error)/nn;
  return error;
}
/*-----------------------------------------------*/
// Move it to utilities
VectorXd arrtoVecXd(double arr[],int nn){
  VectorXd arrvec(nn);
  for (int idim = 0; idim < nn; ++idim){
    arrvec << arr[idim];
  }
  return arrvec;
}
/*-----------------------------------------------*/
// Move it to utilities
double* VecXdtoArr(VectorXd* arrvec, int nn){
  static double *arr {new double[nn]{}};
  for (int idim = 0; idim < nn; ++idim){
    arr[idim] = (*arrvec)(idim);
  }
  return arr;
}
/*-----------------------------------------------*/
int main(){
  // Here I will define whether I am calling the code for fixed point or periodic orbits.
  // Idea is to use the same code for periodic orbit and fixed point both.
  // Fixed point is a periodic orbit with time-period 1 for a map. 
  MV MM;
  double y[ndim],vel[ndim];
  MatrixXd DerM(ndim,ndim);                 // This is better to allocate matrix.
  //Initialize them to zero.
  for (int idim = 0; idim < ndim; ++idim){
    vel[idim]=0;
  }
  // First define all the parameters.
  pre_ode2map(&MM);       
  // Now get the initial configuration (t=0) of the system.
  iniconf(y);   // Implement a program for variable number of arguments.
  //
  /********** Diagnosis start **********/
  // for (int i = 0; i < ndim; ++i){
  //   cout << y[i] << " ";
  // }
  /********** Diagnosis end  **********/
  //
  write_param("wparam.txt");
  if(MM.iorbit){
    if (IsPathExist("VEL")){
      cout << "There is already some data in the folder. I can not replace that.\n"
           << "Please remove the PSI & VEL and run the exec again. ";
      exit(1);
    }
    // The function use Newton-Raphson and calculate the nearest periodic orbit.
    periodic_orbit(&y[0],&vel[0],&MM);
  }
  if(MM.istab){
    if (IsPathExist("eig") || IsPathExist("Eig")){
      cout << "I have already calculated the stability of this orbit."
              " There are already some eigenvalues in the folder. I can not replace that.\n"
              "Please remove the eig file and run the exec again. ";
      exit(1);
    }

    if(MM.iorbit){
      cout << "I shall calculate the stability of the orbit." << endl;
      DerM = Jacobian(y,vel,&MM);
      VectorXcd eivals = DerM.eigenvalues();
      //
      ofstream eigenfile;
      eigenfile.open( "eig" );
      eigenfile << "#------------- Jacobian matris is: -------------#\n" << DerM << endl
                << "\n#------------- The eigenvalues are: -------------#\n" << eivals << endl;
      eigenfile.close();
      //
    }
    else{
      // First check whether you actually have the periodic orbit?
      cout << "Is it a periodic orbit?"
              " (if you don't want this, please comment it out in map_dyn.cpp) " << endl;
      if(IsOrbit(y,&MM)){
        cout << "Yes!!! It is a periodic orbit. I shall calculate stability now." << endl;
        DerM = Jacobian(y,vel,&MM);
        //
        ofstream eigenfile;
        eigenfile.open( "eig" );
        eigenfile << "#------------- Jacobian matris is: -------------#" << endl << DerM << endl;
        //
        VectorXcd eivals = DerM.eigenvalues();
        eigenfile << "\n#------------- The eigenvalues are: -------------#" << endl << eivals << endl;
        eigenfile.close();
      }else{
        cout << "The guess is not periodic orbit." << endl << 
        " Did you cross-check the data or time-period? " << endl <<
        "If yes, use our periodic orbit solver to get nearest orbit to the guess." 
        "Set istab=0, iorbit=1" << endl;
      }
    } 
  }

}