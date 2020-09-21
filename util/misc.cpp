#include <Eigen/Dense>
#include <sys/stat.h>
using namespace Eigen;
/* -----------------------------------------------*/
template<typename T>  // This type of function definition can take any variable type. 
T abs(T value);       // Quite interesting, isn't it ?
double SqEr(double Arr1[], double Arr2[],int ndim);
bool IsPathExist(const std::string &s);
/*-----------------------------------------------*/
double SqEr(double Arr1[], double Arr2[],int nn){
  double error=0;
  for (int cdim=0; cdim < nn; cdim++){
    error=error+(Arr2[cdim]-Arr1[cdim])*(Arr2[cdim]-Arr1[cdim]);
  }
  error = sqrt(error)/nn;
  return error;
}
/*-----------------------------------------------*/
bool IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}
/*-----------------------------------------------*/
template<typename T>
T abs(T value){
  if (value>=0){
    return value;
  }else{
    return value*(-1);
  }
}