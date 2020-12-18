#include <iostream>
#include <fstream>
#include "MapDyn.h"
#include "math.h"
#include <memory.h>
#include <cstdlib>
#include "misc.h"
#include "IO.h"
#include "JacobEigs.h"
#include <Spectra/GenEigsSolver.h>
/*----------------------------------------------*/
typedef Eigen::VectorXcd Vecc;
typedef Eigen::MatrixXcd Matc;
typedef Spectra::DenseGenMatProd<double> MatProd;
typedef Eigen::Map<MatrixXd> MapMat;
int const selectionRule = Spectra::LARGEST_REAL;
/*-----------------------------------------------*/
void calc_stab(double *y, int const neigs = 6);
// void wData_map(string fname_trans, double ytrans[], double fytrans[]);
/*-----------------------------------------------*/
int main(){
  // Here I will define whether I am calling the code for fixed point or periodic orbits.
  // Idea is to use the same code for periodic orbit and fixed point both.
  // Fixed point is a periodic orbit with time-period 1 for a map.
  bool success=1;
  assign_map_param();
  int mapdim = MM.mapdim;
  int period = MM.period;
  double y[ndim],fy[ndim],yall[ndim*period],time[period+1];
  double ytrans[mapdim],fytrans[mapdim],ytrans_all[mapdim*period];
  time[0]=0;
  // First define all the parameters.
  // Now get the initial configuration (t=0) of the system.
  if (MM.guess_space == "Real" || MM.guess_space == "real"){
    iniconf(y);
    coordinate_transform(ytrans,y);
    print(ytrans,mapdim);
  }
  else if(MM.guess_space== "Transformed" || MM.guess_space=="transformed"){
    iniconf(ytrans);
    inv_coordinate_transform(y,ytrans);
    print(y,ndim);
  }
  else{
    cout << "# -- guess_space is not mentioned."
            " Reading input as real space --# " << endl;
    iniconf(y);
    coordinate_transform(ytrans,y);
  }
  memcpy(ytrans_all,ytrans,mapdim*sizeof(double));
  memcpy(yall,y,ndim*sizeof(double));
  // Save parameters in starting.
  write_param("wparam.txt");
  if(MM.iorbit){
    if (IsPathExist("PSI")){
      cout << "# -- There is already some data in the folder. I can not replace that.-- #\n"
           << "# -- Please run make clean or remove the PSI and then run the exec again.-- # ";
      exit(1);
    }else{
      success = periodic_orbit(ytrans_all,fytrans,yall,fy,time);
      if(success){
        // Diagnosis //
        print(fytrans,mapdim);
        cout << "#-- Voila! you found the periodic orbit with period --#" << MM.period << endl;
        cout << "#-- I would go ahead and save it -:) --#" << endl;
        //
        ofstream outfile_trans("PSI_trans",ofstream::out);
        ofstream outfile("PSI",ofstream::out);
        for (int iter = 0; iter < period; ++iter){
          wData(&outfile_trans,ytrans_all+iter*mapdim,time[iter],MM.mapdim,1);
          wData(&outfile,yall+iter*ndim,time[iter]);
        }
        wData(&outfile_trans,fytrans,time[period],MM.mapdim,1);
        wData(&outfile,fy,time[period]);
        outfile.close();
        outfile_trans.close();
        //
      }else{
        cout << "#-- The code to Newton Krylov did not converge. Here are the options: --#\n"
             << "#-- 1) Change the initial guess. 2) Increase the number of trials --#" << endl;
      }
    }
  }
  // Calculate stability
  if(MM.istab){
    if (MM.iorbit){
      if (success){calc_stab(ytrans);}
      else{cout << "#-- Sorry!!! This is not a periodic orbit so it does not make sense. --#\n";}
    }
    else{
     // First check whether you actually have the periodic orbit?
      cout << "Is it a periodic orbit?"
             " (if you don't want this, please comment it out in MapSolve.cpp) " << endl;
      if(IsOrbit(ytrans)){calc_stab(ytrans);}
      // if (1){calc_stab(ytrans);}
      else{
        cout << "The guess is not periodic orbit." << endl << 
        " Did you cross-check the data or time-period? " << endl <<
        "If yes, use our periodic orbit solver to get nearest orbit to the guess." 
        "Set istab=0, iorbit=1" << endl;
      }
    }
  }
}
/*-----------------------------------------------*/
void calc_stab(double *y, int const neigs){
  int mapdim = MM.mapdim;
  Vecc eigval(neigs);
  Matc eigvec(neigs,mapdim);
  if(IsPathExist("eig") || IsPathExist("Eig")){
    cout << "I have already calculated the stability of this orbit."
            "There are already some eigenvalues in the folder. I can not replace that.\n"
            "Please remove the eig file and run the exec again. ";
    exit(1);
  }
 // bool success = jacob_eigval<selectionRule>(&eigval,neigs,mapdim,map_multiple_iter,y);
 bool success = jacob_eigvalvec<selectionRule>(&eigval,&eigvec,neigs,mapdim,map_multiple_iter,y);
 if(!success){
    cout << "Sorry, there is something wrong."
            "I can not converge to perform stability analysis in a matrix free way."
            "I will calculate the complete Jacobian using finite-difference method with the error o(eps^2)"
         << endl;
    // double DerM[ndim][ndim];
    // Jacobian(DerM,y);
    // MapMat Mat_Jacob(DerM,ndim,ndim);
    // MatProd op(Mat_Jacob);
    // GenEigsSolver< double, selectionRule, MatProd > eigs(&op, neigs, 2*neigs+1);
    // eigval = eigs.eigenvalues();
  }else{
    ofstream eigenfile;
    eigenfile.open( "eig" );
    eigenfile << "\n#------------- The eigenvalues are: -------------#\n" << eigval << endl;
    eigenfile << "\n#------------- The eigenvectors are: -------------#\n" << eigvec << endl;
    eigenfile.close();
  }
}
/*-----------------------------------------------*/
// This function should not be needed or should look differently when we implement everything
// as in overshooting method also.
// This function is written in very bad way.
// void wData_map(string fname_trans, double y[]){
//   int time = MM.time;
//   int period = MM.period;
//   MM.time = 0;
//   ofstream outfile(fname_trans,ofstream::out);
//   wData(&outfile,ytrans,0,MM.mapdim,1);
//   for (int iter = 0; iter < period-1; ++iter){
//     map_one_iter(y);

//     wData(&outfile,ytrans,MM.time,MM.mapdim,1);
//   }
//   wData(&outfile,fytrans,time,MM.mapdim,1);
//   outfile.close();
//   MM.time = time;
// }