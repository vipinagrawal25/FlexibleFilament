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
void add(double *added, double *arr1, double *arr2, int nn );
double norm(double *y, int nn );
/* -----------------------------------------------*/
#endif