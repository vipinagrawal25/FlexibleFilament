#include <iostream>
#include <fstream>
#include "MapDyn.h"
#include "math.h"
#include <memory.h>
#include <cstdlib>
#include "misc.h"
#include "IO.h"
/*-----------------------------------------------*/
int main(){
  // Here I will define whether I am calling the code for fixed point or periodic orbits.
  // Idea is to use the same code for periodic orbit and fixed point both.
  // Fixed point is a periodic orbit with time-period 1 for a map. 
  double y[ndim],fy[ndim];
  assign_map_param();
  // First define all the parameters.
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
    if (IsPathExist("PSI")){
      cout << "There is already some data in the folder. I can not replace that.\n"
           << "Please run make clean or remove the PSI and then run the exec again. ";
      exit(1);
    }else{ 
      periodic_orbit(y,fy);
      // test(y,fy);
      print(y,ndim);
      GG(y);
      print(y,ndim);
      // Now save the data;
      // ofstream outfile("PSI",ofstream::out);
      // wData(&outfile,y);
      // wData(&outfile,fy,MM.time);
    }
    // The function use Newton-Raphson and calculate the nearest periodic orbit.
  }
  // if(MM.istab){
  //   if (IsPathExist("eig") || IsPathExist("Eig")){
  //     cout << "I have already calculated the stability of this orbit."
  //             " There are already some eigenvalues in the folder. I can not replace that.\n"
  //             "Please remove the eig file and run the exec again. ";
  //     exit(1);
  //   }
  //   if(MM.iorbit){
  //     double DerM[ndim][ndim];
  //     cout << "I shall calculate the stability of the orbit." << endl;
  //     Jacobian(DerM[][ndim],y);
  //     VectorXcd eivals = DerM.eigenvalues();
  //     //
  //     ofstream eigenfile;
  //     eigenfile.open( "eig" );
  //     eigenfile << "#------------- Jacobian matris is: -------------#\n" << DerM << endl
  //               << "\n#------------- The eigenvalues are: -------------#\n" << eivals << endl;
  //     eigenfile.close();
  //     //
  //   }
  //   else{
  //     // First check whether you actually have the periodic orbit?
  //     cout << "Is it a periodic orbit?"
  //             " (if you don't want this, please comment it out in map_dyn.cpp) " << endl;
  //     if(IsOrbit(y)){
  //       cout << "Yes!!! It is a periodic orbit. I shall calculate stability now." << endl;
  //       Jacobian(DerM[][ndim],y);
  //       //
  //       ofstream eigenfile;
  //       eigenfile.open( "eig" );
  //       eigenfile << "#------------- Jacobian matris is: -------------#" << endl << DerM << endl;
  //       //
  //       VectorXcd eivals = DerM.eigenvalues();
  //       eigenfile << "\n#------------- The eigenvalues are: -------------#" << endl << eivals << endl;
  //       eigenfile.close();
  //     }else{
  //       cout << "The guess is not periodic orbit." << endl << 
  //       " Did you cross-check the data or time-period? " << endl <<
  //       "If yes, use our periodic orbit solver to get nearest orbit to the guess." 
  //       "Set istab=0, iorbit=1" << endl;
  //     }
  //   } 
  // }
}
