#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "misc.h"
#include "model.h"
#include <iomanip>
#include "input.h"
using namespace std;
/*-----------------------------------------------*/
void wData(ofstream *fptr, double *y, double time, int Npt, int ppt){
  // If we do not want to write velocity in a file.
  int ndimt = Npt*ppt;
  switch(wDataMeth){
    case 1:
      for(int ip = 0; ip < Npt; ++ip){
        for (int jp = 0; jp < ppt; ++jp){
          *fptr << setprecision(15) << y[ppt*ip+jp];
          *fptr << '\t';
        }
        *fptr << '\n';
      }
      break;
    case 2:
      (*fptr) << time << "\t";
      for (int idim = 0; idim < ndimt; ++idim){
        *fptr << setprecision(15) << y[idim] << "\t";
      }
      *fptr <<  endl;
      *fptr << "#------------------------------------------------#" << endl;
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
void wData(ofstream *fptr, ofstream *fptr_vel, double *y, double *vel, double time,
           int Npt, int ppt){
  int ndimt = Npt*ppt;
  switch(wDataMeth){
    case 1:
      for(int ip = 0; ip < Npt; ++ip){
        // *fptr << setprecision(15) << 0;
        // *fptr << setprecision(15) << '\t';
        for (int jp = 0; jp < ppt; ++jp){
          *fptr << setprecision(15) << y[ppt*ip+jp];
          *fptr << setprecision(15) << '\t';
        }
        // Now next three numbers contain values of velocity.
        // *fptr << setprecision(15) << 0;
        // *fptr << setprecision(15) << '\t';
        for (int jp = 0; jp < ppt; ++jp){
          *fptr << setprecision(15) << vel[ppt*ip+jp];
          *fptr << '\t';
        }
        *fptr << '\n';
      }
      break;

    case 2:
      (*fptr) << setprecision(15) << time << "\t";
      (*fptr_vel) << setprecision(15) << time << "\t";
      for (int idim = 0; idim < ndimt; ++idim){
        *fptr << setprecision(15) << y[idim] << "\t";
        *fptr_vel << setprecision(15) << vel[idim] << "\t";
      }
      *fptr <<  endl;
      *fptr << "#------------------------------------------------#" << endl;

      *fptr_vel <<  endl;
      *fptr_vel << "#------------------------------------------------#" << endl;
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/********************************************/
void rData(double *y,string fname, int Npt, int ppt){
  double ch;
  int cnt=0;
  ifstream myfile;
  //
  double ndimt = Npt*ppt;
  switch(rDataMeth){
    case 1:
      myfile.open(fname);
      while(myfile >> ch){cnt++;}
      myfile.close();
      cout << "# Number of elements read: " << cnt << endl;
      if(cnt==Npt){}
      else if(cnt==ndimt){
        cnt=0;
        myfile.open(fname);
        while(myfile >> ch){y[cnt]=ch;cnt++;}
        myfile.close();
      }
      else if(cnt==pp*ndim){
        myfile.open(fname);
        for(int idim = 0; idim < Npt; ++idim){
          for (int ip = 0; ip < ppt; ++ip){
            myfile >> y[idim*ppt+ip];
          }
          for (int ip = 0; ip < pp; ++ip){
            myfile >> ch;
          }
        }
        myfile.close();
      }
      break;
    // case 2:
    //   if(){
    //     l = "data/";
    //     l.append(filename);
    //   }else{
    //     l=filename;
    //   }
    //   while( getline(myfile,token) ){
    //     line=token;
    //     getline(myfile,token);          // Dumping this: #---------
    //   }
    //   // Now convert all the tab separated entries to array.
    //   iss.str(line);
    //   // getline(iss, token, '\t');     // First entry is time, delete that.
    //   for (int idim = 0; idim < ndimt; ++idim){
    //     getline(iss, token, '\t');      // Get next token.
    //     y[idim]=stod(token);            // Convert to double and store it in y.
    //   }
    //   myfile.close();
    //   break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData( double y[], string filename) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/*-----------------------------------------------*/
// Wrapper function for wData, in case y and vel contains multiple rows.
// We do not pass the file pointer, but the fname.
// void wData(string fname, double y[][ndim], double time[], int nrow){
//   switch(wDataMeth){
//     case 1:{
//       // for (int irow = 0; irow < nrow; ++irow){
//       //   ofstream outfile("output/var"+to_string()+".txt",fstream::out);
//       //   for(int ip = 0; ip < Np; ++ip){
//       //     for (int jp = 0; jp < pp; ++jp){
//       //       file << y[irow][pp*ip+jp];
//       //       file << '\t';
//       //     }
//       //     // Now just throw away next three numbers as they contain values of velocity.
//       //     for (int jp = 0; jp < pp; ++jp){
//       //       file << vel[irow][pp*ip+jp];
//       //       file << '\t';
//       //     }
//       //     file << '\n';
//       //   }
//       // }
//       cout << "This is not necessary. I am EXITING from function wData." << endl;
//       break;
//     }
//     case 2:{
//       ofstream outfile(fname,ofstream::out);
//       for (int irow = 0; irow < nrow; ++irow){
//         wData(&outfile,y[irow],time[irow]);
//       }
//       outfile.close();
//       break;
//       // for (int idim = 0; idim < ndim; ++idim){
//       //   *fptr << y[idim] << "\t";
//       //   *fptr_vel << vel[idim] << "\t";
//       // }
//       // *fptr <<  "\n";
//       // *fptr << "#------------------------------------------------#\n";

//       // *fptr_vel <<  "\n";
//       // *fptr_vel << "#------------------------------------------------#\n";
//       // break;
//     }
//     default:
//       cout << "Hey, your choice of writing data does not exist. "
//               "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
//               "in model.cpp file." << endl;
//       exit(1);
//   }
// }
// /*-----------------------------------------------*/
// void rData(ifstream *fptr, double *y){
//   double num=0.0;
//   for(int ip = 0; ip < Np; ++ip){
//     *fptr >> y[3*ip];
//     *fptr >> y[3*ip+1];
//     *fptr >> y[3*ip+2];
//     // Now just throw away next three numbers as they contain values of velocity.
//     *fptr >> num;
//     *fptr >> num;
//     *fptr >> num;
//   }
// }