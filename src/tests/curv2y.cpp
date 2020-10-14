#include <iostream>
#include "model.h"
#include "misc.h"
using namespace std;
int main(){
	double y[ndim],y_trans[Np],y2[ndim];
	iniconf(y);
	print(y,ndim);
	coordinate_transform(y_trans, y);
	print(y_trans,Np);
	inv_coordinate_transform(y2, y_trans);
	print(y2,ndim);
	return 0;
}