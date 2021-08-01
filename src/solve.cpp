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
#include "IO.h"
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
  int ndim_tr = ntracer*ptracer;
  //
  double y[ndim],y_prev[ndim],vel[ndim],EForceArr[ndim];
  double CurvSqr[Np],SS[Np];
  double time,time_prev,timer_global;
  double y_tr[ndim_tr],vel_tr[ndim_tr];
  double yzero[2],yone[2];
  double velzero[2] = {0.,0.};
  int filenumber;
  string lastline;
  int ldiagnos=0;
  int tdiagnos = 1;
  // But these things are only for diagnostic purpose and it's too much hassle for someting that is not even 
  // important. So we wil keep it in the wrong way, as most of time this would give right result.
  double dt_min = 10;
  //----------------------------
  int itn=1; 
  clock_t timer;
  // For storing the Mean square displacement of the rod with timer, every row would have different MSD wrt 
  // time for a particular value of AA.
  check_param();
  //
  if (time){
    // read data from output/time.txt
    // Let's just say that we are saving data as the same interval as earlier.
    iniconf(y,&time);
  }
  if (time==0){
    iniconf(y);
  }
  // Deleting contents of the folder and creating folder again.
  // exit(1);
  eval_rhs(time,y,vel,tdiagnos,CurvSqr,SS,EForceArr);
  //
  system("exec mkdir output");
  ofstream outfile;
  outfile.open("output/var0.txt");
  if (bcb==2){
    calc_yone(yone,time);
    calc_yzero(yzero,time);
    //
    wData(&outfile,&outfile,yzero,velzero,time,1,pp);
    wData(&outfile,&outfile,yone,velzero,time,1,pp);
  }
  wData(&outfile,&outfile,y,vel);                                     // Code it in your model.cpp
  outfile.close();
  //
  if(itracer){
    iniconf_tr(y_tr);
    // memcpy(&y_tr[ndim_tr],&y_tr[0],ndim_tr*sizeof(double));
    eval_rhs_tr(time,EForceArr,y,y_tr,vel_tr);
    //
    outfile.open("output/tracer0.txt");
    wData(&outfile,&outfile,y_tr,vel_tr,time,ntracer,ptracer);         // Code it in your model.cpp
    outfile.close();
  }
  //
  // Initializing curv square and ss. It won't depend on the configuration number.
  // Call a diagnostic function for all these things. Save names in model file.
  for (int ip = 0; ip < Np; ++ip){
    CurvSqr[ip] = 0;
    SS[ip] = 0;
  }
  /* Opening every file again in mode. This thing does not depend on configuration number and that's why 
     it is outside the loop */
  // fstream outfile_MSD("MSD.txt", ios::app);
  fstream outfile_time("output/time.txt", ios::app);
  fstream outfile_curvature("output/curvature.txt", ios::app);
  fstream outfile_SS("output/material_point.txt", ios::app);
  timer = clock();
  timer_global = timer/CLOCKS_PER_SEC;
  while(time < TMAX){
    //euler(pdim,&y[irb],time-dt,dt);
    //rnkt2(pdim,&y[irb],time-dt,dt);
    // rnkt4(pdim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos);
    memcpy(y_prev,y,ndim*sizeof(double));
    time_prev = time;
    rnkf45(ndim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], EForceArr,tdiagnos);
    if (itracer){
      euler_tr(ndim_tr,y_prev,y_tr,vel_tr,time_prev,dt,EForceArr);
    }
    // DP54(pdim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos);
    // cout << time << endl;
    if (dt<dt_min){
      dt_min = dt;
    }
    if (time+dt>=tdiag*filenumber && time<tdiag*filenumber){tdiagnos=1;}
    else{tdiagnos=0;}
    // cout << time << endl;
    if (time>=tdiag*filenumber){
      //
      string l = "output/var" + to_string(filenumber) + ".txt";
      outfile.open(l, ios::out);
      if (bcb==2){
        calc_yzero(yzero,time);
        calc_yone(yone,time);
        wData(&outfile,&outfile,yzero,velzero,time,1,pp);
        wData(&outfile,&outfile,yone,velzero,time,1,pp);
      }
      wData(&outfile,&outfile,y,vel);
      outfile.close();
      //
      if(itracer){
        string l_tr = "output/tracer" + to_string(filenumber) + ".txt";
        outfile.open(l_tr, ios::out);
        wData(&outfile,&outfile,y_tr,vel_tr,time,ntracer,ptracer);
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
      cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<TMAX<<"\n";
      tdiagnos=0;
      itn=itn+1;
    }
  }
  timer_global = clock()/CLOCKS_PER_SEC - timer_global;
  outfile_time.close();
  outfile_curvature.close();
  outfile_SS.close();
  //
  cout << "Total number of iteration: " << itn << endl;
  // cout << "Total time elapsed: " << ((double)timer)/CLOCKS_PER_SEC << "s" << endl;
  cout << "Total time elapsed: " << timer_global << "s" << endl;
  cout << "Minimum value of dt: " << dt_min << endl;  
  // cout << "Difference between max length and Minimum length of the rod: " << lengthmax-lengthmin << endl;
  // cout << "Max Length: " << lengthmax << '\t' << "Min Length: " <<lengthmin << endl;
  // cout << "The average change in the length of the rod is: " << sqrt(MSElen)/itn << endl;
// cout << filenumber-1 << endl;
//----------------------------
}
/********************************************/
// if(niniconf==-1){
//   filenumber = lastfile+1;
//   // -----------------------------------------------
//   double num = 0.0;
//   string line;
//   ifstream outfile_time("output/time.txt");
//   fstream outfile_time_new("output/time_new.txt", ios::out);   
//   Just a new time file which would be renamed later anyway.
//   // ifstream outfile_MSD("MSD.txt");
//   // fstream outfile_MSD_new("MSD_new.txt",ios::out);
//   ifstream outfile_curvature("output/curvature.txt");
//   fstream outfile_curvature_new("output/curvature_new.txt",ios::out) ;
//   ifstream outfile_SS("output/material_point.txt");
//   fstream outfile_SS_new("output/material_point_new.txt",ios::out) ;
//   int ifile = 1;
//   // This loop just reads data from existing time file and dump it into the new file till it is allowed.
//   // This is necessary as you can ask to delete data after a certain iteration.
//   while(ifile<lastfile){
//     if (outfile_time >> num){
//       time = num;
//       outfile_time_new << num << endl;
//       //
//       getline(outfile_curvature,line,'\n');
//       outfile_curvature_new << line << endl;
//       //
//       getline(outfile_SS,line,'\n');
//       outfile_SS_new << line << endl;
//       ifile = ifile+1;
//     }
//     else{
//       cout << "ERROR: The last code has not been run till the file mentioned.";
//       cout <<  "This does not make sense." << endl;
//       return 0;
//     }
//   }
//   // Closing all the files which we have opened.
//   outfile_time.close();
//   // outfile_MSD.close();
//   outfile_curvature.close();
//   outfile_SS.close();
//   // Now close the new files as well
//   outfile_time_new.close();
//   // outfile_MSD_new.close();
//   outfile_curvature_new.close();
//   outfile_SS_new.close();
//   // Closing files done
//   system("exec rm -f output/time.txt");
//   system("exec mv output/time_new.txt output/time.txt");
//   // fstream outfile_time("output/time.txt", ios::app);    // Now opening file again in append mode.
//   // system("exec rm -f MSD.txt");
//   // system("exec mv MSD_new.txt MSD.txt");
//   // fstream outfile_time("MSD.txt", ios::app);            // Now opening file again in append mode. 
//   system("exec rm -f output/curvature.txt");
//   system("exec mv output/curvature_new.txt output/curvature.txt");
//   // fstream outfile_curvature("output/curvature.txt", ios::app);    
    // Now opening file again in append mode.
//   system("exec rm -f output/material_point.txt");
//   system("exec mv output/material_point_new.txt output/material_point.txt");
//   // fstream outfile_SS("output/material_point.txt", ios::app);    
     // Now opening file again in append mode.      
//   // ------------------------------------------------------------------------------------------------
//   /*Code to remove the contents of the output folder after the last file mentioned. 
//   The code would have already returned an error message if lastfile is more than the 
    // total number of files present*/
//   // The idea is that I check whether the file exists or not, if not come out of loop immediately, 
//   // else remove the file if it has index more than lastfile mentioned in '.h' file.
//   ifile = lastfile+1;
//   // cout << ifile << endl;
//   string filename = "output/var" + to_string(ifile) + ".txt";
  
//   struct stat buffer;
//   while(!stat(filename.c_str(), &buffer)){
//     string removefile = "exec rm -f " + filename ;
//     system(removefile.c_str());
//     string filename = "output/var" + to_string(ifile) + ".txt";
//     ifile = ifile+1;
//   }
//   //Store the values of initial file into y0.
//   // ifstream initialfile("output/var0.txt", ios::in);
//   //keep storing values from the text file so long as data exists:
//   // rData(&initialfile,&y0[0]);     // code it in your model.cpp file.
//   // initialfile.close();
// }
// else{