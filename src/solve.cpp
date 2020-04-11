#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string>
#include<vector>
#include "ode.h"
#include "model.h"
#include "input.h"
#include<sys/stat.h>
//#include "cuda.h"
using namespace std;
/**************************/
unsigned int const ndim=Nensemble*pdim;
/* ----------------------------------------*/
int main()
{
  double y[ndim], y0[ndim], CurvSqr[Np], SS[Np],time,MeanSqDis,timer_global;
  double vel[ndim];
  int filenumber;
  string lastline;
  int ldiagnos=0;
  int tdiagnos = 1;
  // remember that the lengthmin and lengthmax might be different when the code is getting started in the middle i.e. conf_number==-1
  // But these things are only for diagnostic purpose and it's too much hassle for someting that is not even important. So we wil keep 
  // it in the wrong way, as most of time this would give right result.
  double lengthmin=height;
  double lengthmax=height;
  double dt_min = 10;
  double gamma = 8*M_PI*viscosity*ShearRate*pow(height,4)/(AA);
  double MSElen = 0;
  // This gamma is non-dimensional number which describes the competition between viscous forces and elastic forces.
  // cout << HH << endl;
  //  setup_cuda();
//----------------------------
  int itn=1; 
  clock_t timer;
  // For storing the Mean square displacement of the rod with timer, every row would have different MSD wrt time 
  // for a particular value of AA.
  if (dd>aa){
      cout << "ERROR: The diameter of a particle should be less than the distance between two particles." << endl;
      return 0;
  }
  
  if(lastfile){
    filenumber = lastfile+1;
    // -----------------------------------------------------------------------------------------------------
    double num = 0.0;
    string line;
    ifstream outfile_time("output/time.txt");
    fstream outfile_time_new("output/time_new.txt", ios::out);   // Just a new time file which would be renamed later anyway.
    ifstream outfile_MSD("MSD.txt");
    fstream outfile_MSD_new("MSD_new.txt",ios::out);
    ifstream outfile_curvature("output/curvature.txt");
    fstream outfile_curvature_new("output/curvature_new.txt",ios::out) ;
    ifstream outfile_SS("output/material_point.txt");
    fstream outfile_SS_new("output/material_point_new.txt",ios::out) ;
    int ifile = 1;
    // This loop just make reads data from existing time file and dump it into the new file till it is allowed.
    // This is necessary as you can ask to delete data after a certain iteration.
    while(ifile<lastfile){
      if (outfile_time >> num){
        time = num;
        outfile_time_new << num << endl;
        // 
        outfile_MSD >> num;
        outfile_MSD_new << num << endl;
        //
        getline(outfile_curvature,line,'\n');
        outfile_curvature_new << line << endl;
        //
        getline(outfile_SS,line,'\n');
        outfile_SS_new << line << endl;
        ifile = ifile+1;
      }
      else{
        cout << "ERROR: The last code has not been run till the file which has been mentioned. This does not make sense." << endl;
        return 0;
      }
    }
    // Closing all the files which we have opened.
    outfile_time.close();
    outfile_MSD.close();
    outfile_curvature.close();
    outfile_SS.close();
    // Now close the new files as well
    outfile_time_new.close();
    outfile_MSD_new.close();
    outfile_curvature_new.close();
    outfile_SS_new.close();
    // Closing files done
    system("exec rm -f output/time.txt");
    system("exec mv output/time_new.txt output/time.txt");
    // fstream outfile_time("output/time.txt", ios::app);    // Now opening file again in append mode.
    system("exec rm -f MSD.txt");
    system("exec mv MSD_new.txt MSD.txt");
    // fstream outfile_time("MSD.txt", ios::app);    // Now opening file again in append mode. 
    system("exec rm -f output/curvature.txt");
    system("exec mv output/curvature_new.txt output/curvature.txt");
    // fstream outfile_curvature("output/curvature.txt", ios::app);    // Now opening file again in append mode.
    system("exec rm -f output/material_point.txt");
    system("exec mv output/material_point_new.txt output/material_point.txt");
    // fstream outfile_SS("output/material_point.txt", ios::app);    // Now opening file again in append mode.      
    // -----------------------------------------------------------------------------------------------------------
    /*Code to remove the contents of the output folder after the last file mentioned. The code would have already returned an error message
    if lastfile is more than the total number of files present*/
    // The idea is that I check whether the file exists or not, if not come out of loop immediately, else remove the file if it has index more than
    // lastfile mentioned in '.h' file.
    ifile = lastfile+1;
    // cout << ifile << endl;
    string filename = "output/var";
    filename.append(to_string(ifile));
    filename.append(".txt");
    struct stat buffer;
    while(!stat(filename.c_str(), &buffer)){
      string removefile = "exec rm -f ";
      removefile.append(filename);
      system(removefile.c_str());
      string filename = "output/var";
      filename.append(to_string(ifile));
      filename.append(".txt");
      ifile = ifile+1;
    }
    //Store the values of initial file into y0.
    ifstream initialfile("output/var0.txt", ios::in);
    //keep storing values from the text file so long as data exists:
    for (int ip = 0; ip < Np; ++ip){
      initialfile >> y0[3*ip];
      initialfile >> y0[3*ip+1];
      initialfile >> y0[3*ip+2];
   // Now just throw away next three numbers as they contain values of velocity.
      initialfile >> num;
      initialfile >> num;
      initialfile >> num;
    }
  }
  else{
    filenumber = 1;
    time = 0;
    // Deleting contents of the folder and creating folder again.
    system("exec rm -rf output");
    system("exec mkdir output");
    system("exec rm -f MSD.txt");
    iniconf(y0, vel, conf_number);
    ofstream outfile;
    outfile.open("output/var0.txt"); 
    for (int ip = 0; ip < Np; ++ip){
      outfile << y0[3*ip] << '\t' << y0[3*ip+1] << '\t' << y0[3*ip+2] << '\t' 
              << vel[3*ip] << '\t' << vel[3*ip+1]  << '\t' << vel[3*ip+2] << endl ; 
    }
  } 
  // Initializing curv square and ss. It won't depend on the configuration number.
  for (int ip = 0; ip < Np; ++ip){
    CurvSqr[ip] = 0;
    SS[ip] = 0;
  }
  /*Opening every file again in append mode. This thing does not depend on configuration number and that's why it is outside the loop*/
  fstream outfile_MSD("MSD.txt", ios::app);
  fstream outfile_time("output/time.txt", ios::app);
  fstream outfile_curvature("output/curvature.txt", ios::app);
  fstream outfile_SS("output/material_point.txt", ios::app);
  iniconf(y, vel, conf_number);
  timer = clock();
  timer_global = timer/CLOCKS_PER_SEC;
  while(time < TMAX){
    //euler(pdim,&y[irb],time-dt,dt);
    //rnkt2(pdim,&y[irb],time-dt,dt);
    // rnkt4(pdim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos);
  	rnkf45(pdim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos); 
    // DP54(pdim, &y[0], &vel[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos); 
    // cout << time << endl;
    if (dt<dt_min){
      dt_min = dt;
    }
    // MSElen = MSElen+(height-SS[Np-1])*(height-SS[Np-1]);
    if (SS[Np-1]>lengthmax){
        lengthmax=SS[Np-1];
    }
    if (SS[Np-1]<lengthmin){
        lengthmin=SS[Np-1];
        // cout << lengthmin << endl;
    }
    tdiagnos = 0;
    // cout << time << endl;
    if (time>=tdiag*filenumber){
      // cout << dt << endl;
      ofstream outfile;
      string l = "output/var";
      l.append(to_string(filenumber));
      l.append(".txt");
      outfile.open(l, ios::out);
      outfile_curvature << time << '\t' ;
      outfile_time << time;
      for (int ip = 0; ip < Np; ++ip){
        outfile << y[3*ip] << '\t' << y[3*ip+1] << '\t' << y[3*ip+2] << '\t' 
                << vel[3*ip] << '\t' << vel[3*ip+1] << '\t' << vel[3*ip+2] << endl;
        /* Non-dimensionalizing the co-ordinate with respect to the height of the rod*/
        outfile_curvature << CurvSqr[ip]*aa*aa << '\t';  /*Square of curvature is non-dimensionalized with the multiplication of square of 
                                                             bead distance */   
        outfile_SS << SS[ip]/SS[Np-1] << '\t';
      }
      MeanSqDis = 0;
      for (int idim = 0; idim < ndim; ++idim){
          MeanSqDis = MeanSqDis+(y[idim]-y0[idim])*(y[idim]-y0[idim]);
      }
      MeanSqDis = sqrt(MeanSqDis);
      outfile_MSD << MeanSqDis;
      outfile_curvature << endl;
      outfile_MSD << endl;
      outfile_time << endl;
      outfile_SS << endl;
      outfile.close(); 

      filenumber = filenumber+1;
      cout<<"Done, time="<<time << "\t dt=" << dt <<"\t TMAX="<<TMAX<<"\n";
      tdiagnos =1;
    }
    itn=itn+1;
  }
// timer = clock()-timer;
  timer_global = clock()/CLOCKS_PER_SEC - timer_global;  
  outfile_time.close();
  outfile_curvature.close();
  outfile_SS.close();
  // outfile_MSD << endl;
  outfile_MSD.close();

  ofstream outfile_information;
  outfile_information.open("../info.csv", ios::out | ios::app);
  outfile_information << itn << "," <<  timer_global << "," << TMAX << ',' << dt_min << "," << viscosity << ','
  << ShearRate << ',' <<  omega << "," << Np << "," << AA << "," << HH << "," << dd << "," << height << "," << sigma << "," 
  << gamma << "," << HH*aa*aa/AA << endl;
  outfile_information.close();

  outfile_information.open("info.txt", ios::out | ios::app);
  outfile_information << itn << '\t' <<  timer_global << '\t' << TMAX << '\t' << dt_min << '\t' << viscosity << '\t'
  << ShearRate << '\t' <<  omega << '\t' << Np << '\t' << AA << '\t' << HH << '\t' << dd << '\t' << height << '\t' << sigma << '\t' 
  << gamma << '\t' << HH*aa*aa/AA << endl;
  outfile_information.close();

  cout << "Total number of iteration: " << itn << endl;
  // cout << "Total time elapsed: " << ((double)timer)/CLOCKS_PER_SEC << "s" << endl;
  cout << "Total time elapsed: " << timer_global << "s" << endl;
  cout << "Minimum value of dt: " << dt_min << endl;  
  // cout << "Difference between max length and Minimum length of the rod: " << lengthmax-lengthmin << endl;
  cout << "Max Length: " << lengthmax << '\t' << "Min Length: " <<lengthmin << endl;
  // cout << "The average change in the length of the rod is: " << sqrt(MSElen)/itn << endl;
  cout << "Gamma: " << gamma << '\t' << "K: " << HH*aa*aa/AA << endl;
// cout << filenumber-1 << endl;
//----------------------------
}
/* ----------------------------------------*/
