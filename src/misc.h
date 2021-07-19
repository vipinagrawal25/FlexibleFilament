#include <string>
#ifndef FILE_MISC_SEEN
#define FILE_MISC_SEEN
using namespace std;
/* -----------------------------------------------*/
// template<typename T>  // This type of function definition can take any variable type. 
// T abs(T value);       // Quite interesting, isn't it ?
double SqEr(double Arr1[], double Arr2[],int ndim);
bool IsPathExist(const std::string &s);
void check_param() ;
void write_param(string fname);
void print(double*, int);
void print(double*, int start, int skip,int end);
void add(double *added, double *arr1, double *arr2, int nn);
void substract(double *subsed, double *arr1, double *arr2, int nn);
double norm(double *y, int nn );
double norm(const double *y, int nn );
void downScale(double *yscaled, double *y, int factor, int nn );
void zeros(double *yzero, int ndim);
// void print(vec2 *arr, int nn);
/* -----------------------------------------------*/
#endif