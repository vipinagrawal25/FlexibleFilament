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
/**************************/
unsigned int const ndim=Nensemble*pdim;
/* ----------------------------------------*/
int main(){
  double y[ndim], y0[ndim], CurvSqr[Np], SS[Np],time,MeanSqDis,timer_global;
  int filenumber;

  int ldiagnos=0;
  int tdiagnos = 1;

  // remember that the lengthmin and lengthmax might be different when the code is getting started in the middle i.e. conf_number==-1
  // But these things are only for diagnostic purpose and it's too much hassle for someting that is not even important. So we wil keep 
  // it in the wrong way, as most of time this would give right result.
  
  double lengthmin=height;
  double lengthmax=height;
  
  double dt_min = 10;

  double gamma = 8*M_PI*viscosity*aa*aa*aa*ShearRate*height/AA ;
  double MSElen = 0;
  // This gamma is non-dimensional number which describes the competition between viscous forces and elastic forces.
  
  // cout << HH << endl;
  //  setup_cuda();
//----------------------------
  int itn=1; 

  clock_t timer;

  // For storing the Mean square displacement of the rod with timer, every row would have different MSD wrt time 
  // for a particular value of AA.

  fstream outfile_MSD("MSD.txt", ios::app);

  if (dd>aa){
      cout << "ERROR: The diameter of a particle should be less than the distance between two particles." << endl;
      return 0;
  }

  iniconf(y, conf_number);

  if (conf_number == -1)
  {
      filenumber = lastfile+1;
      time = 40;

      // This whole thing has been written to calculate the time at which the last code was stopped.
      ifstream input_time ("output/time.txt",ios::in);
      // input_time.open();
      if(input_time.is_open()) 
      {
        input_time.seekg(-1,ios_base::end);                // go to one spot before the EOF

        bool keepLooping = true;
        while(keepLooping) 
        {
            char ch;
            input_time.get(ch);                            // Get current byte's data

            if((int)input_time.tellg() <= 1) {             // If the data was at or before the 0th byte
                input_time.seekg(0);                       // The first line is the last line
                keepLooping = false;                // So stop there
            }
            else if(ch == '\n') {                   // If the data was a newline
                keepLooping = false;                // Stop at the current position.
            }
            else {                                  // If the data was neither a newline nor at the 0 byte
                input_time.seekg(-2,ios_base::cur);        // Move to the front of that data, then to the front of the data before it
            }
        }

        string lastline;
        getline(input_time,lastline);                      // Read the current line
        time = stod(lastline);

        input_time.close();   // close the file and open it in append mode later.
      // }

      ifstream myfile("output/position1.txt");
      string line;
      // getline(myfile,line);
      // cout << "Aisa ho hi nahi sakta ki ye yaha na aaye" << endl;
      int ip =0;
      while ( getline (myfile,line,'\t') )
      {
          // cout << "Yadi ye yaha nahi aa raha hai to iska matlab ye while loop ne gandagi faila rakhi hai" <<endl;
          // cout << line << '\n';
          y0[ip] = stod(line);
      }
      myfile.close();
  }
  else
  {
      filenumber = 1;
      time = 0;
      // Deleting contents of the folder and creating folder again.
      system("exec rm -rf output");
      system("exec mkdir output");
      system("exec rm -f MSD.txt");

      iniconf(y0, conf_number);

      ofstream outfile;
      outfile.open("output/position0.txt"); 

      for (int ip = 0; ip < Np; ++ip)
      {
          outfile << y0[3*ip] << '\t' << y0[3*ip+1] << '\t' << y0[3*ip+2] << endl ; 
      }

      for (int ip = 0; ip < Np; ++ip)
      {
        CurvSqr[ip] = 0;
        SS[ip] = 0;
      }
  }


  fstream outfile_time("output/time.txt", ios::app);
  fstream outfile_curvature("output/curvature.txt", ios::app);
  fstream outfile_SS("output/material_point.txt", ios::app);
  ifstream myfile ("output/position1.txt", ios::in);

  timer = clock();
  timer_global = timer/CLOCKS_PER_SEC;



  while(time < TMAX)
  {
    // cout << time << endl;
    // ldiagnos=itn%idiag;
    // tdiagnos=time%tdiag;
    // time = time + dt;
    // for(int ibody=0;ibody<Nensemble;ibody++){
    // int irb=pdim*ibody;
    //euler(pdim,&y[irb],time-dt,dt);
    //rnkt2(pdim,&y[irb],time-dt,dt);
    // rnkt4(pdim,&y[irb],time-dt,dt);

  	// timer = clock();

    rnkf45(pdim, &y[0], &time, &dt, &CurvSqr[0], &SS[0], tdiagnos);           

    // timer = clock() - timer;
    // timer_global = timer_global + timer;
    // cout << timer << endl;
    if (dt<dt_min)
    {
      dt_min = dt;
    }

    MSElen = MSElen+(height-SS[Np-1])*(height-SS[Np-1]);

    if (SS[Np-1]>lengthmax)
    {
        lengthmax=SS[Np-1];
    }

    if (SS[Np-1]<lengthmin)
    {
        lengthmin=SS[Np-1];
        // cout << lengthmin << endl;
    }

    // cout << dt << endl;
    // }
    tdiagnos = 0;
    // cout << "Yaar ye code chal kyu nahi raha hai " << endl;
    if (time<=tdiag*filenumber && time+dt>=tdiag*filenumber) 
    {
       // cout << time << '\t' << y[0] << '\t' << (sin(2*time+10*sin(0.1*time)))/sqrt(6+3*cos(0.1*time)) << '\t' << 1/sqrt(6+3 *cos(0.1*time))<<endl;
      // cout << dt << endl;
      ofstream outfile;
      string l = "output/position";
      l.append(to_string(filenumber));
      l.append(".txt");
      outfile.open(l, ios::out);

      outfile_curvature << time*omega/M_PI << '\t' ;
      outfile_time << time*omega/M_PI;

      for (int ip = 0; ip < Np; ++ip)
      {
        outfile << y[3*ip] << '\t' << y[3*ip+1] << '\t' << y[3*ip+2] << endl ; 

        /* Non-dimensionalizing the co-ordinate with respect to the height of the rod*/
        
        outfile_curvature << CurvSqr[ip]*aa*aa << '\t' ;  /*Square of curvature is non-dimensionalized with the multiplication of square of 
                                                             bead distance */   
        outfile_SS << SS[ip] << '\t';        
        
    
        // MeanSqDis = MeanSqDis+(y[3*idim]-y0[idim])*(y[idim]-y0[idim]);

        // outfile_MSD << MeanSqDis << ';' ; 
        // cout << CurvSqr[ip] << endl;
      }

      if (SaveInfo == 'Y')
      {
          MeanSqDis = 0;
          for (int idim = 0; idim < ndim; ++idim)
          {
              MeanSqDis = MeanSqDis+(y[idim]-y0[idim])*(y[idim]-y0[idim]);
          }
          MeanSqDis = sqrt(MeanSqDis);
          outfile_MSD << MeanSqDis;
      }
       
      outfile_curvature << endl;
      outfile_MSD << endl;
      outfile_time << endl;
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
  timer_global = clock()/CLOCKS_PER_SEC - timer_global;  
  outfile_time.close();
  outfile_curvature.close();
  outfile_SS.close();

  if (SaveInfo=='Y')
  {
    outfile_MSD << endl;
    outfile_MSD.close();

    ofstream outfile_information;
    outfile_information.open("../info.csv", ios::out | ios::app);
    outfile_information << itn << ";" <<  ((double)timer_global)/CLOCKS_PER_SEC << ";" << TMAX << ';' << dt_min << ";" << viscosity << ';'
    << ShearRate << ';' <<  omega << ";" << Np << ";" << AA << ";" << HH << ";" << dd << ";" << height << ";" << sigma << ";" 
    << gamma << ";" << HH*aa*aa/AA << endl;
  }


  cout << "Total number of iteration: " << itn << endl;
  // cout << "Total time elapsed: " << ((double)timer)/CLOCKS_PER_SEC << "s" << endl;
  cout << "Total time elapsed: " << timer_global << "s" << endl;
  cout << "Minimum value of dt: " << dt_min << endl;  
  // cout << "Difference between max length and Minimum length of the rod: " << lengthmax-lengthmin << endl;
  cout << "Max Length: " << lengthmax << '\t' << "Min Length: " <<lengthmin << endl;
  cout << "The average change in the length of the rod is: " << sqrt(MSElen)/itn << endl;

// cout << filenumber-1 << endl;
//----------------------------
}
/* ----------------------------------------*/
