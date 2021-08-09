#ifndef FILE_ESTRING2DIM_SEEN
#define FILE_ESTRING2DIM_SEEN
#include "math.h"
#include <string>
/**************************/
struct MParam{
	double height;		// height of the box we are doing simulations in.
	unsigned int Np;	// total number of points on the rod.
	unsigned int pp;
	unsigned int ndim;				
	double dd; 			// diameter of the filament.
						// The particles would also have same diameter. 
	double viscosity ;	// Equivalent to kinematic viscosity of glycerin
	double ShearRate ;
	double angularVel ;   
	double period;
	double AA;
	double HH ;			// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
	char UseRP ;		// N -- Local drag, Y -- Global drag
	/* BC Not needed in Peskin way*/
	int bcb ;      		// Boundary condition at bottom, 0 - fixed, 1 - free, 2 - clamped.
	int bct;       		// Boundary condition at top. 1 - free.
	//
	int icons;			// 0 -- Norain, 1 -- circle, 2 -- '8'.
	int np_cons;		// Number of points torain?
	int loc_con[np_cons];
	//
	int iext_force; 	// 0 -- No force, 1 means Peskin way to clamp points
	int fnp ;
	int floc[fnp] ;		// External force location
	//
	int iext_flow ;   /* External flow: 1->shear rate changes as step function. 3-> Sine type shear rate.*/
	int niniconf ;    /* Configuration of the system at t = 0. 0 -> sine perturbation in  the  filament. 
						1 -> Straight filament, 2-> , -1 -> read from a file*/
	int ievolve_save ;
	string datafile ;	// If you want to read input from some file, Mention the file name.
	int itracer ;		// 0 -> No tracers, 1 -> tracers on square lattice, 3 -> tracer on circular lattice? 
	int np_tracer ;
	int pp_tracer ;
}
/**************************/
extern MParam param;
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