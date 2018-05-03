#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN
/*---------------------------------------*/
#include<math.h>
using namespace std;
/* ----------------------------------------*/
unsigned int const Np=100;// total number of points on the rod.
unsigned int const pdim=3*Np; 
double const AA = 0.001;		  // A is bending rigidity constant (A = EI), unit -> Pa.m^4
double const aa =  1.0/(double)(Np-1); 	// distance between two nodes,
                                                 //which should be constant accodrding to the constraint.
double const HH = 16*AA/(aa*aa);			// H should be atleast 1500 times of A. Follow: bit.ly/2r23lmA unit -> Pa.m^2
double const OneByGamma=1.;
double const Z0=0.;
double const FF0 = 0;			// Force Value on the ends
/* ----------------------------------------*/
void iniconf(double *y);
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */
