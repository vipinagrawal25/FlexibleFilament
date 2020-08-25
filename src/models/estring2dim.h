#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN
#include "math.h"
/**************************/
using namespace std;
/* ----------------------------------------*/
double const height = 1.28;					// height of the box we are doing simulations in.
unsigned int const Np= 8;					// total number of points on the rod.
unsigned int const pp=2;
unsigned int const ndim=pp*Np;				
double const aa = height/(double)(Np-1); 	// distance between two nodes.
double const dd = 0.005; 					// diameter of the filament.
											// The particles would also have same diameter. 
double const viscosity = 10;				// Equivalent to kinematic viscosity of glycerin
double const Z0=0.;		 					// If we want the bottom point of the rod to be fixed.
double const FFY0 = 0;	  					// Amplitude of the force.
// Sigma is a dimensionless number, which is described as frequency parameter.
double const sigma=0.75;					
double const ShearRate = 2;
double const omega = ShearRate*sigma;
// double const AA = 1.5*pow(10,-5)*pow(height,4)*100;
double const AA = 1.5*pow(10,-5)*pow(height,4)*100;
double const HH = 16*AA/(dd*dd);			// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
char const UseRP = 'Y';						// Y is for the new one, O is for old one and N is for not using 
											// Rotne Pragor tensor. 
// This sets up the initial configuration of the system. For more information go through the Readme file.
bool const bcb = 1 ;      					// Boundary condition at bottom, 0 means fixed, 1 means free.
bool const bct = 1 ;       					// Boundary condition at top.  
bool const iext_force = 0; 					// Whether to apply the external force or not
int const floc=Np-1;						// External force location
int const iext_flow = 3;   					/* External flow: 1->shear rate changes as step function.
						  					3-> Sine type shear rate.*/
int const niniconf = 1;    	 				/* Configuration of the system at t = 0. 0 -> sine perturbation in the filament.
											1 -> Straight filament, 2-> */
int const lastfile= 0;						// If you want to read input from some file, Mention the file number.
/* ----------------------------------------*/
void iniconf(double *y); 	// The configuration number is defined for different
											// Initial configuration into the system.
void rData(ifstream *fptr,double *y);
void check_param();
void write_param(string fname);
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */