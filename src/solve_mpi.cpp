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
#include <mpi.h>
/********************************************/
using namespace std;
/* ----------------------------------------*/
void read_param(string param_name);
void read_evolve(string fname);
void write_param(string fname);
/* ----------------------------------------*/
double period, AA, HH, TMAX;
int TotalFiles;
string datafile;
/* ----------------------------------------*/
int main(int argc, char** argv){
  /*-------------MPI part starts-------------------------*/
  int first_folder=1;
  MPI_Init(&argc, &argv);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  /*-------------MPI part ends-------------------------*/
  string run_dir = "run" + to_string(world_rank+first_folder) + "/";
  cout << "# I have access to " + run_dir << endl;
  //
  system( ("exec mkdir "+run_dir+"output").c_str() );
  fstream outfile_terminal(run_dir+"terminal.out", ios::app);
  fstream outfile_time(run_dir+"output/time.txt", ios::app);
  fstream outfile_curvature(run_dir+"output/curvature.txt", ios::app);
  fstream outfile_SS(run_dir+"output/material_point.txt", ios::app);
  //
  pid_t pid = getpid();
  outfile_terminal << "# ID for this process is: " << pid << endl;
  int ndim_tr = np_tracer*pp_tracer;
  //
  double y[ndim],y_prev[ndim],vel[ndim],EForceArr[ndim],sigma,facAA;
  //
  read_param(run_dir+"rparam.txt");
  read_evolve(run_dir+"revolve.txt");
  write_param(run_dir+"wparam.txt");
  //
  double CurvSqr[Np],SS[Np];
  double time,time_prev,timer_global;
  double y_tr[ndim_tr],vel_tr[ndim_tr];
  double yzero[2],yone[2];
  double velzero[2] = {0.,0.};
  int filenumber=1;
  string lastline;
  int ldiagnos=0;
  int tdiagnos = 1;
  double tdiag = TMAX/TotalFiles;
  // But these things are only for diagnostic purpose and it's too much hassle for someting that is not even 
  // important. So we wil keep it in the wrong way, as most of time this would give right result.
  ofstream outfile;
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
  }
  if (time==0){
    iniconf(y);
  }
  // Deleting contents of the folder and creating folder again.
  // exit(1);
  eval_rhs(time,y,vel,tdiagnos,CurvSqr,SS,EForceArr);
  //
  if (ievolve_save){
      outfile.open(run_dir+"output/var0.txt");
      wData(&outfile,&outfile,y,vel);                                     // Code it in your model.cpp
      outfile.close();
  }
  //
  if(itracer){
    iniconf_tr(y_tr);
    // memcpy(&y_tr[ndim_tr],&y_tr[0],ndim_tr*sizeof(double));
    eval_rhs_tr(time,EForceArr,y,y_tr,vel_tr);
    //
    outfile.open(run_dir+"output/tracer0.txt");
    wData(&outfile,&outfile,y_tr,vel_tr,time,np_tracer,pp_tracer);         // Code it in your model.cpp
    outfile.close();
  }
  // Initializing curv square and ss. It won't depend on the configuration number.
  // Call a diagnostic function for all these things. Save names in model file.
  for (int ip = 0; ip < Np; ++ip){
    CurvSqr[ip] = 0;
    SS[ip] = 0;
  }
  /* Opening every file again in mode. This thing does not depend on configuration number and that's why 
     it is outside the loop */
  // fstream outfile_MSD("MSD.txt", ios::app);
  timer = clock();
  timer_global = timer/CLOCKS_PER_SEC;
  while(time < TMAX){
    memcpy(y_prev,y,ndim*sizeof(double));
    time_prev = time;
    rnkf45(ndim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], EForceArr,tdiagnos);
    if (itracer){
      euler_tr(ndim_tr,y_prev,y_tr,vel_tr,time_prev,dt,EForceArr);
    }
    if (dt<dt_min){
      dt_min = dt;
    }
    if (time-start_time+dt>=tdiag*filenumber && time-start_time<tdiag*filenumber){tdiagnos=1;}
    else{tdiagnos=0;}
    //
    if (time-start_time>=tdiag*filenumber){
      //
      if (ievolve_save){
        string l = run_dir+"output/var" + to_string(filenumber) + ".txt";
        outfile.open(l, ios::out);
        wData(&outfile,&outfile,y,vel);
        outfile.close();
      }
      //
      if(itracer){
        string l_tr = run_dir+"output/tracer" + to_string(filenumber) + ".txt";
        outfile.open(l_tr, ios::out);
        wData(&outfile,&outfile,y_tr,vel_tr,time,np_tracer,pp_tracer);
        outfile.close();
      }
      //
      /* Call a function to write both diagnostic variable.*/
      outfile_curvature << time << '\t' ;
      outfile_time << time;
      for (int ip = 0; ip < Np; ++ip){
        /* Non-dimensionalizing the co-ordinate with respect to the height of the rod */
        outfile_curvature << CurvSqr[ip]*aa*aa << '\t';
       /*Square of curvature is non-dimensionalized with the multiplication
          of square of bead distance */
        outfile_SS << SS[ip]/SS[Np-1] << '\t';
      }
      /*---------------------------------- */
      outfile_curvature << endl;
      outfile_time << endl;
      outfile_SS << endl;
      filenumber = filenumber+1;
      outfile_terminal<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<TMAX<<"\n";
      tdiagnos=0;
      itn=itn+1;
    }
  }
  // Make a file to let the world know that you are done with the code //
  fstream outfile_done("done", ios::out);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  //
  outfile_terminal << "Total number of iteration: " << itn << endl;
  outfile_terminal << "Total time elapsed: " << timer_global << "s" << endl;
  outfile_terminal << "Minimum value of dt: " << dt_min << endl; 
  //----------------------------
  timer_global = clock()/CLOCKS_PER_SEC - timer_global;
  outfile_time.close();
  outfile_curvature.close();
  outfile_SS.close();
  outfile_terminal.close();
  MPI_Finalize();
}
/********************************************/
void read_param(string fname){
  double sigma, facAA;
  // string datafile;
  ifstream myfile;
  myfile.open(fname);
  myfile>>sigma;
  myfile>>facAA;
  myfile>>datafile;
  // setting up the parameters //
  period = 2*M_PI/(sigma*ShearRate);
  AA = 1.5*pow(10,-5)*pow(height,4)*facAA;
  HH = 16*AA/(dd*dd);
  //
  cout << "sigma = " << sigma << ";\tfacAA = " << facAA << ";\tdatafile = "<< datafile <<endl;
  myfile.close();
}
/********************************************/
void read_evolve(string fname){
  ifstream myfile;
  myfile.open(fname);
  myfile>>TMAX;
  myfile>>TotalFiles;
  cout << "TMAX = " << TMAX << "; TotalFiles = " << TotalFiles << endl;
  myfile.close();
}
/********************************************/
void  __attribute__((weak)) write_param( string fname ){
  cout << "The user has not implmented way to write data!!!" << endl;
}