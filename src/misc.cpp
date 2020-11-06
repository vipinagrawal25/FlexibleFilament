#include <iostream>
#include <sys/stat.h>
#include "math.h"
#include <string>
using namespace std;
/*-----------------------------------------------*/
double SqEr(double Arr1[], double Arr2[],int nn){
  double error=0;
  for (int cdim=0; cdim < nn; cdim++){
    error=error+(Arr2[cdim]-Arr1[cdim])*(Arr2[cdim]-Arr1[cdim]);
  }
  error = sqrt(error);
  return error;
}
/*-----------------------------------------------*/
bool IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}
/*-----------------------------------------------*/
// template<typename T>
// T abs(T value){
//   if (value>=0){
//     return value;
//   }else{
//     return value*(-1);
//   }
// }
/******************************************************/
void __attribute__((weak)) check_param(){
  cout << "I believe all your model parameters are physical. Otherwise, define function: "
          "void check_param() in model.cpp file" << endl;
}
// -----------------------------------------------
void __attribute__((weak)) write_param(string fname){
  cout << "I can not find any implementation to write model parameters." 
          "Hence I will not write anything." << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int start, int skip, int end){
  for (int in = 0; in < int(end/skip); ++in){
      cout << arr[skip*in+start] << " ";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int nn){
  for (int in = 0; in < nn; ++in){
      cout << arr[in] << "\t";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void add(double *added, double *arr1, double *arr2, int nn ){
  for (int ii = 0; ii < nn; ++ii){
    added[ii] = arr1[ii] + arr2[ii];
  }
}
/*-----------------------------------------------*/
void substract(double *subsed, double *arr1, double *arr2, int nn ){
  for (int ii = 0; ii < nn; ++ii){
    subsed[ii] = arr1[ii] - arr2[ii];
  }
}
/*-----------------------------------------------*/
double norm(double *y, int nn ){
  double norr =0;
  for (int ii = 0; ii < nn; ++ii){
    norr = norr+y[ii]*y[ii];
  }
  return sqrt(norr);
}
/*-----------------------------------------------*/
// template<typename T, typename ... Args>
// T First(T first,Args ... args ){
//   return first;
// }
// /*-----------------------------------------------*/
// template<typename T, typename ... Args>
// T Last(T first, T second, Args ... args ){
//   return Last(args...);
// }
// /*-----------------------------------------------*/
// double upScale(double *y, int factor, int nn ){
//   for (int in = 0; in < nn; ++in){
//     y[in] = y[in]*factor;
//   }
// }
/*-----------------------------------------------*/
double downScale(double *yscaled, double *y, int factor, int nn ){
  for (int in = 0; in < nn/factor; ++in){
    yscaled[in] = y[in*factor];
  }
}