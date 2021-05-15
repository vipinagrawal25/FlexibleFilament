#include <iostream>
#include "model.h"
#include "misc.h"
#include <memory.h>
#include <fstream>
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
	// y[0] = 0.859517625;
	// y[1] = 0.41718115212;
	// y[2] = 0.84005228779;
	// y[3] = 0.4642291456;
	// y[0] = 0.86;
	// y[1] = 0.43;
	// y[2] = 0.84;
	// y[3] = 0.46;
	ifstream myfile;
	myfile.open("PSI");
	double ch;
 	myfile >> ch;
 	int cnt =0;
  	while(myfile >> ch){y[cnt]=ch;cnt++;}
  	cout << "# Number of elements read: " << cnt << endl;
	// y[0] = (lambda[0]-1)/lambda[0];
}
/*-----------------------------------------------*/
void definelambda(){
	for (int ip = 0; ip < Np; ++ip){
		lambda[ip]=3.455;
	}
}
/* ----------------------------------------------- */
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
	eval_rhs(y);
}
/* ----------------------------------------------- */