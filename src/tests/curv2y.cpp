#include <iostream>
#include "model.h"
#include "misc.h"
using namespace std;
int main(){
	double y[ndim],kappa[Np],y2[ndim],kappa2[Np];
	iniconf(y);
	print(y,0,2,ndim);
	print(y,1,2,ndim);
	coordinate_transform(kappa,y);
	print(kappa,Np);
	inv_coordinate_transform(y2,kappa);
	print(y2,0,2,ndim);
	print(y2,1,2,ndim);
	coordinate_transform(kappa2,y2);
	print(kappa2,Np);
	return 0;
}