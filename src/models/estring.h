#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN

/**************************/
using namespace std;
/* ----------------------------------------*/
unsigned int const Np= 30;			// total number of points on the rod.
unsigned int const pdim=3*Np;
double const height = 1;					// height of the box we are doing simulations in.
double const aa =  height/(double)(Np); 	// distance between two nodes.
// AA = AA/aa ;							// See the equation number 35 in DOI: 10.1063/1.1848511

double const factor1 = 	2*pow(10,-5);	// unit -> Pa.m^2, AA = factor1*N*l^2
/* This factor depends only on material and not on anything else. So for a particular material, irrespective of
it's length factor1 would be constant iff mass per unit length is kept constant.			
checout file https://drive.google.com/file/d/1_Ig_sJLZF71APnuuX74-eByuhI4O6atD/view?usp=sharing */

double const AA = factor1*Np*height*height;		 	// AA is bending rigidity constant (AA = EI), unit -> Pa.m^4
double const dd = 0.1*aa;				// This is also not needed in the case, if back reaction from the fluid is neglected.
double const viscosity = 5.;				// Equivalent to kinematic viscosity of glycerin
double const HH = 16*AA*(1./(aa*aa));		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
// The below line is commented because of some error. So I included it in last formulae itself.
// HH =  HH/aa; 			// This is done to model the rod perfectly (Also see: https://goo.gl/nbELSx Week 27/2018 DS)
// double const HH = 0.;

double const OneByGamma=1.;				// Mobility is just the inverse of Gamma. So OneByGamma is actually the mobility of the system if no
										// back reaction is taken into the account.
double const Z0=0.;						// If we want the bottom point of the rod to be fixed.
// double const FFZ0 = 1*height;		// Force Value on the ends
double const FFZ0 = 0;					// Different force for different configuration.

double const omega = 0.4;
char const UseRP = 'Y';					// Y is for the new one, O is for old one and N is for not using 
										// Rotne Pragor tensor. 
int const conf_number = 1;
double const ShearRate = 0.5;

/* ----------------------------------------*/
void iniconf(double *y, int conf_number); // The § number is defined for different
											// Initial configuration into the system.
void diagnos(double time, int itn, double y[]);
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */
