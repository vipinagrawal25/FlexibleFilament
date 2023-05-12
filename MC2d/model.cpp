#include "model.h"
/* -----------------------------------------------*/
void set_param(MESH *mesh){
	cout << "Reading and setting up the parameters." << endl;
	int Np=16;
	int pp=3;
	(*mesh).Np=Np;
	(*mesh).pp=pp;
	(*mesh).ndim=Np*pp;
}
/* -----------------------------------------------*/
void iniconf(MESH *mesh, int size_mesh){
	if((*mesh).ndim==0){set_param(mesh);}
	
}