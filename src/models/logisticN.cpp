#include <iostream>
#include "model.h"
#include "misc.h"
#include <memory.h>
using namespace std;
/* -----------------------------------------------*/
double lambda[Np];
string cc;
/*-----------------------------------------------*/
void eval_rhs(double *y){
	double ytemp[ndim];
	memcpy(ytemp,y,ndim*sizeof(double));
	for (int ip = 0; ip < Np; ++ip){
		y[ip] = lambda[ip]*ytemp[ip]*(1-ytemp[ip]);
	}
}
/*-----------------------------------------------*/
void iniconf(double *y){
	definelambda();
	y[0] = 0.611260467;
	y[1] = 0.950484434;
	y[2] = 0.188255099;
	y[3] = 0.611260467;
	y[4] = 0.950484434;
	y[5] = 0.188255099;

	for (int ip = 6; ip < Np; ++ip){
		y[ip]= 0;
	}
	// y[0] = (lambda[0]-1)/lambda[0];
}
/*-----------------------------------------------*/
void definelambda(){
	for (int ip = 0; ip < Np; ++ip){
		lambda[ip]=4;
	}
	lambda[0]=4;
}
/* ----------------------------------------------- */
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
	eval_rhs(y);
}
/* ----------------------------------------------- */