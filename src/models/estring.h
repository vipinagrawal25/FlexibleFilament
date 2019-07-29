#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN
// #include <../input.h>
/**************************/
using namespace std;
/* ----------------------------------------*/
unsigned int const Np=101;			// total number of points on the rod.
unsigned int const pdim=3*Np;
double const height = 1;					// height of the box we are doing simulations in.
double const aa = height/(double)(Np-1); 	// distance between two nodes.

double const dd = 0.005*height;				// r/l ratio for the rod has been kept constant. It should be noted that the particles would also have same diameter. 

double const viscosity = 10;				// Equivalent to kinematic viscosity of glycerin

// AA is bending rigidity constant (AA = EI), unit -> Pa.m^4, 3.14/(16*4) = 0.05

double const Z0=0.;						// If we want the bottom point of the rod to be fixed.
// double const FFZ0 = 1*height;		// Force Value on the ends
double const FFZ0 = 0;					// Different force for different configuration.

// Sigma is a dimensionless number, which is described as frequency parameter.
double const sigma=1.5;					
double const ShearRate = 1;
double const omega = ShearRate*sigma;

double const AA = 1.5*pow(10,-4)*0.1;
double const HH = 64*AA/(aa*aa);		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2

char const UseRP = 'Y';					// Y is for the new one, O is for old one and N is for not using 
										// Rotne Pragor tensor. 

// This sets up the initial configuration of the system. For more information go through the Readme file.
int const conf_number = 2;	
int const lastfile=0;

char const SaveInfo = 'Y';	// This decides whether the simulation is important enough to save the information.
								// Set N for not saving.

/* ----------------------------------------*/
void iniconf(double *y, int conf_number); 	// The configuration number is defined for different
											// Initial configuration into the system.
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */
