#include <iostream>
#include "model.h"
#include "misc.h"
#include <memory.h>
using namespace std;
/* ----------------------------------------------- */
void eval_rhs(double *y){
	double ytemp[ndim];
	memcpy(ytemp,y,ndim*sizeof(double));
	y[0] = 1.0-aa*ytemp[0]*ytemp[0] + ytemp[1];
	y[1] = bb*ytemp[0];
}
void iniconf(double *y){
	y[0]=1;
	y[1]=0;
}
/* ----------------------------------------------- */