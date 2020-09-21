#include <iostream>
#include <fstream>
#include "map_dyn.h"
#include "math.h"
#include <memory.h>
#include <cstdlib>
/*-----------------------------------------------*/
int main(){
  // Here I will define whether I am calling the code for fixed point or periodic orbits.
  // Idea is to use the same code for periodic orbit and fixed point both.
  // Fixed point is a periodic orbit with time-period 1 for a map. 
  MV MM;
  double y[ndim];
  assign_map_param(&MM);
  // First define all the parameters.
  // Now get the initial configuration (t=0) of the system.
  // iniconf(y);   // Implement a program for variable number of arguments.
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
      double DerM[ndim][ndim];
      cout << "I shall calculate the stability of the orbit." << endl;
      Jacobian(DerM[][ndim],y,vel,&MM);
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
