#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string>
#include "ode.h"
#include "model.h"
#include "input.h"
//#include "cuda.h"
using namespace std;
void diagnos(int p);
/**************************/
unsigned int const ndim=Nensemble*pdim;
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double CurvSqr[Np], SS[Np];
  for (int ip = 0; ip < Np; ++ip)
  {
    CurvSqr[ip] = 0;
    SS[ip] = 0;
  }
  double time=0.;
  int k;
  int ldiagnos=0;
  iniconf(y);
  // cout << HH << endl;
  //  setup_cuda();
//----------------------------
  int itn=1;  
  // int temp = log10(TMAX/dt);
  // int temp2 = pow(10,temp+1);
  
  // Deleting contents of the folder and creating folder again.
  system("exec rmtrash output");    
  system("exec mkdir output");
  ofstream outfile_time;
  outfile_time.open("output/time.txt", ios::out); 

  ofstream outfile_curvature;
  outfile_curvature.open("output/curvature.txt", ios::out);

  ofstream outfile_SS;
  outfile_SS.open("output/material_point.txt", ios::out);  

  while(time < TMAX){
    ldiagnos=itn%idiag;

    // for(int ibody=0;ibody<Nensemble;ibody++){
    time = time+dt;
    // int irb=pdim*ibody;
    //euler(pdim,&y[irb],time-dt,dt);
    //rnkt2(pdim,&y[irb],time-dt,dt);
    // rnkt4(pdim,&y[irb],time-dt,dt);

    rnkf45(pdim, &y[0], time-dt, &dt, &CurvSqr[0], &SS[0], ldiagnos);
    // cout << dt << endl;
    // }
    if (ldiagnos == 0) {
      outfile_time << time << '\t' ;
       // cout << time << '\t' << y[0] << '\t' << (sin(2*time+10*sin(0.1*time)))/sqrt(6+3*cos(0.1*time)) << '\t' << 1/sqrt(6+3 *cos(0.1*time))<<endl;
      // cout << dt << endl;
      ofstream outfile;
      string l = "output/position";
      l.append(to_string(itn/idiag));
      l.append(".txt");
      outfile.open(l, ios::out);

      outfile_curvature << time << '\t' ;

      for (int ip = 0; ip < Np; ++ip)
      {
        outfile << y[3*ip] << '\t' << y[3*ip+1] << '\t' << y[3*ip+2] << endl ;
        outfile_curvature << CurvSqr[ip] << '\t' ;
        outfile_SS << SS[ip] << '\t';
        // cout << CurvSqr[ip] << endl;
      }

      outfile_curvature << endl;
      outfile_SS << endl;
      outfile.close();

      // for (int ip = 0; ip < Np; ++ip)
      // {
      //   outfile << y[3*ip+2] << '\t';
      // }
      // outfile << endl;
  }

  cout<<"Done, time="<<time<<"\t TMAX="<<TMAX<<"\n";
  itn=itn+1;

}
outfile_time.close();
outfile_curvature.close();
outfile_SS.close();
cout << itn/idiag << endl;
//----------------------------
}
/* ----------------------------------------*/
