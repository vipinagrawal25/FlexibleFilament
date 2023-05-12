#ifndef UTIL_H
#define UTIL_H
#include "TriTri_intersection/src/tri_tri_intersect.c"
#include "../include/global.h"
int tri_tri_intersect(POSITION *first_tri, POSITION *second_tri){
	double uu0[3],uu1[3],uu2[3],vv0[3],vv1[3],vv2[3];
	Pos2y(uu0,first_tri[0],0);
	Pos2y(uu1,first_tri[1],0);
	Pos2y(uu2,first_tri[2],0);
	//
	Pos2y(vv0,second_tri[0],0);
	Pos2y(vv1,second_tri[1],0);
	Pos2y(vv2,second_tri[2],0);
	return tri_tri_intersection_test_3d(uu0,uu1,uu2,vv0,vv1,vv2);
}
void Pos2y(double *y, POSITION XX, int ip){
  y[3*ip] = XX.x;
  y[3*ip+1] = XX.y;
  y[3*ip+2] = XX.z;
}
#endif