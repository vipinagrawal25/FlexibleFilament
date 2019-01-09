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
  int tdiagnos = 0;
  iniconf(y, conf_number);
  double dt_min = 10;
  // cout << HH << endl;
  //  setup_cuda();
//----------------------------
  int itn=1;  
  int filenumber=1;
  clock_t timer;
  int timer_global;
  // int temp = log10(TMAX/dt);
  // int temp2 = pow(10,temp+1);
  
  // Deleting contents of the folder and creating folder again.
  system("exec rm -r output");    
  system("exec mkdir output");
  ofstream outfile_time;
  outfile_time.open("output/time.txt", ios::out);

  ofstream outfile_curvature;
  outfile_curvature.open("output/curvature.txt", ios::out);

  ofstream outfile_SS;
  outfile_SS.open("output/material_point.txt", ios::out);  

  // timer = clock();
  while(time < TMAX){
    // ldiagnos=itn%idiag;
    // tdiagnos=time%tdiag;
    // time = time + dt;
    // for(int ibody=0;ibody<Nensemble;ibody++){
    // int irb=pdim*ibody;
    //euler(pdim,&y[irb],time-dt,dt);
    //rnkt2(pdim,&y[irb],time-dt,dt);
    // rnkt4(pdim,&y[irb],time-dt,dt);

  	timer = clock();
    rnkf45(pdim, &y[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos);
    timer = clock() - timer;
    timer_global = timer_global + timer;
    // cout << timer << endl;
    if (dt<dt_min)
    {
        dt_min = dt;
    }
    // cout << dt << endl;
    // }
    tdiagnos = 0;
    if (time<=tdiag*filenumber && time+dt>=tdiag*filenumber) {
      outfile_time << time << '\t';
       // cout << time << '\t' << y[0] << '\t' << (sin(2*time+10*sin(0.1*time)))/sqrt(6+3*cos(0.1*time)) << '\t' << 1/sqrt(6+3 *cos(0.1*time))<<endl;
      // cout << dt << endl;
      ofstream outfile;
      string l = "output/position";
      l.append(to_string(filenumber));
      l.append(".txt");
      outfile.open(l, ios::out);

      outfile_curvature << time << '\t' ;

      for (int ip = 0; ip < Np; ++ip)
      {
        outfile << y[3*ip]/height << '\t' << y[3*ip+1]/height << '\t' << y[3*ip+2]/height << endl ;   
        /* Non-dimensionalizing the co-ordinate with respect to the height of the rod*/
        
        outfile_curvature << CurvSqr[ip]*aa*aa << '\t' ;  /* Square of curvature is non-dimensionalized with the multiplication of square of 
                                                             bead distance */ 
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
      filenumber = filenumber+1;
      cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<TMAX<<"\n";
      tdiagnos =1;
  }
  itn=itn+1;
}
// timer = clock()-timer;
outfile_time.close();
outfile_curvature.close();
outfile_SS.close();

ofstream outfile_information;
outfile_information.open("info.txt", ios::out | ios::app);
outfile_information << itn << '\t' << '\t' <<  ((double)timer_global)/CLOCKS_PER_SEC << '\t' << dt_min << '\t' << TMAX << '\t' << '\t' <<
viscosity << endl;

cout << "Total number of iteration: " << itn << endl;
// cout << "Total time elapsed: " << ((double)timer)/CLOCKS_PER_SEC << "s" << endl;
cout << "Total time elapsed: " << (double)timer_global/CLOCKS_PER_SEC << "s" << endl;
cout << "Minimum value of dt: " << dt_min << endl;
// cout << filenumber-1 << endl;
//----------------------------
}
/* ----------------------------------------*/
