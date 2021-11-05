#include <iostream>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include "ode.h"
#include "model.h"
#include "input.h"
#include <sys/stat.h>
#include "constant.h"
#include "io.h"
#include "misc.h"
#include <fstream>
#include <unistd.h>
#include <memory.h>
/********************************************/
using namespace std;
/* ----------------------------------------*/
int main(){
  pid_t pid = getpid();
  cout << "# ID for this process is: " << pid << endl;
  double y[ndim],y_prev[ndim],vel[ndim],EForceArr[ndim];
  double time,time_prev,timer_global;
  int filenumber=1;
  // But these things are only for diagnostic purpose and it's too much hassle for someting that is not even 
  // important. So we wil keep it in the wrong way, as most of time this would give right result.
  double dt_min = 10;
  double dt = start_dt;
  //----------------------------
  int itn=1;
  clock_t timer;
  // For storing the Mean square displacement of the rod with timer, every row would have different MSD wrt 
  // time for a particular value of AA.
  check_param();
  //
  time = start_time;
  if (time){
    // read data from output/time.txt
    // Let's just say that we are saving data as the same interval as earlier.
    iniconf(y,&time,tdiag);
    filenumber = (int) (time/tdiag);
  }
  if (time==0){
    iniconf(y);
  }
  // Deleting contents of the folder and creating folder again.
  // exit(1);
  eval_rhs();
  //
  system("exec mkdir output");
  //
  outfile.open("output/var0.txt");
  wData(&outfile,&outfile,y,vel);                                     // Code it in your model.cpp
  outfile.close();
  /*Opening every file again in mode. This thing does not depend on configuration number and that's why 
    it is outside the loop */
  fstream outfile_time("output/time.txt", ios::app);
  timer = clock();
  timer_global = timer/CLOCKS_PER_SEC;
  if (time>0){
    outfile_time << endl;
  }
  while(time < TMAX){
    memcpy(y_prev,y,ndim*sizeof(double));
    time_prev = time;

    rnkf45(ndim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], endlForceArr,tdiagnos);
    if (dt<dt_min){
      dt_min = dt;
    }
    if (time+dt>=tdiag*filenumber && time<tdiag*filenumber){tdiagnos=1;}
    else{tdiagnos=0;}
    if (time>=tdiag*filenumber){
      string l = "output/var" + to_string(filenumber) + ".txt";
      outfile.open(l, ios::out);
      wData(&outfile,&outfile,y,vel);
      outfile.close();
      /* Call a function to write both diagnostic variable.*/
      outfile_time << time;
      filenumber = filenumber+1;
      cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<TMAX<<"\n";
      tdiagnos=0;
      itn=itn+1;
    }
  }
  timer_global = clock()/CLOCKS_PER_SEC - timer_global;
  outfile_time.close();
  //
  cout << "Total number of iteration: " << itn << endl;
  cout << "Total time elapsed: " << timer_global << "s" << endl;
  cout << "Minimum value of dt: " << dt_min << endl;  
}
/********************************************/
void __attribute__((weak)) check_param(){
  cout << "You have implemented nothing for checking the parameters" << endl;
}
/********************************************/