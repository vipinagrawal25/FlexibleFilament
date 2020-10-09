#include <iostream>
#include "model.h"
#include "misc.h"
#include <memory.h>
using namespace std;
/* -----------------------------------------------*/
double lambda[Np];
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
	for (int ip = 0; ip < Np; ++ip){
		y[ip]=0;
	}
}
/*-----------------------------------------------*/
void definelambda(){
	for (int ip = 0; ip < Np; ++ip){
		lambda[ip]=0.5;
	}
	lambda[0]=1.5;
}
/* -----------------------------------------------*/