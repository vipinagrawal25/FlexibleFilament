#include <iostream>
#include "model.h"
#include "misc.h"
#include <fstream>
using namespace std;
int main(){
	double y_trans[Np],y[ndim];
	// iniconf(y_trans);
	// print(y_trans,Np);
	iniconf(y);
	print(y,0,2,ndim);
	print(y,1,2,ndim);
	coordinate_transform(y_trans,y);
	// print(y_trans,Np);
	inv_coordinate_transform(y,y_trans);
	print(y,0,2,ndim);
	print(y,1,2,ndim);
	// inv_coordinate_transform(y,y_trans);
	// print(y,ndim);
	// substract(kappa2, kappa2, kappa, Np);
	// print(kappa2,Np);
	return 0;
}