#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN
#include "math.h"
/**************************/
using namespace std;
/* ----------------------------------------*/
unsigned int const Np= 1;					// total number of points on the rod.
unsigned int const pp= 2;
unsigned int const ndim=pp*Np;
double const aa = 1.4;
double const bb = 0.3;	
double const omega = 1;			// It is not used anywhere.	
/* ----------------------------------------*/
void iniconf(double *y); 	// The configuration number is defined for different
							// Initial configuration into the system.
void eval_rhs(double *y);
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */