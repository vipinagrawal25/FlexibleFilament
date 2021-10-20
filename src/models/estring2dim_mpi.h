#ifndef FILE_ESTRING2DIM_SEEN
#define FILE_ESTRING2DIM_SEEN
#include "math.h"
#include <string>
/**************************/
extern double sigma;
extern double facAA;
/**************************/
double const height = 1.28;					// height of the box we are doing simulations in.
unsigned int const Np= 256;					// total number of points on the rod.
unsigned int const pp= 2;
unsigned int const ndim=pp*Np;				
double const aa = height/(double)(Np-1); 	// distance between two nodes.
double const dd = 0.005; 					// diameter of the filament.
double const viscosity = 10;				// Equivalent to kinematic viscosity of glycerin
double const ShearRate = 2;
double const angularVel = 3;
double const period = 2*M_PI/(sigma*ShearRate);
double const AA = 1.5*pow(10,-5)*pow(height,4)*facAA;
double const HH = 16*AA/(dd*dd);			// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
char const UseRP = 'Y';						// N -- Local drag, Y -- Global drag
/* BC Not needed in Peskin way*/
double const bcb = 1;      					// Boundary condition at bottom, 0 - fixed, 1 - free, 2 - clamped.
double const bct = 1;       				// Boundary condition at top. 1 - free.
//
int const icons = 0;						// 0 -- No constrain, 1 -- circle, 2 -- '8'.
int const np_cons = 2;						// Number of points to constrain?
int const loc_con[np_cons] = {0,1};
//
int const iext_force = 0; 					// 0 -- No force, 1 means Peskin way to clamp points
int const fnp = 2;
int const floc[fnp] = {0,1};					// External force location
//
int const iext_flow = 3;   /* External flow: 1->shear rate changes as step function. 3-> Sine type shear rate.*/
int const niniconf = 1;    /* Configuration of the system at t = 0. 0 -> sine perturbation in  the  filament. 
							1 -> Straight filament, 2-> , -1 -> read from a file*/
int const ievolve_save=0;
string const datafile = "sol1";				// If you want to read input from some file, Mention the file name.
int const itracer=2;						// 0 -> No tracers, 1 -> tracers on square lattice, 3 -> tracer on circular lattice? 
int const np_tracer=1024;
int const pp_tracer=3;
/**************************/
void iniconf(double *y); 	// The configuration number is defined for different
							// Initial configuration into the system.
void iniconf(double *y, double *aTime, double tdiag); 	
void iniconf_tr(double *y_tr);
void check_param();
void write_param(string fname);
void eval_rhs(double time, double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]);
void eval_rhs(double time, double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[],
			 				double EForceArr[], bool iEforceArr=1);
void coordinate_transform(double *y_trans, double *y);
void inv_coordinate_transform(double *y, double *y_trans);
void pre_next_iter(double *y, double *y_trans);
void eval_rhs_tr(double time,double EForceArr[],double y[],double y_tr[],double rhs[]);
void calc_yone(double *yone, double time);
void calc_yzero(double *yzero, double time);
void set_param();
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */
/* ----------------------------------------*/