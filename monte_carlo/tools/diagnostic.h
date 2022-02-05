#ifndef DIAGNOSTIC
#define DIAGNOSTIC
#include "../utils/util.h"
#include "../include/global.h"
//
void mesh_intersection(POSITION *pos, int *triangles, int N ){
	// The function	will print the triangle numbers which are intersecting.
	POSITION first_tri[3], second_tri[3];
	for (int it = 0; it< 3*N; it=it+3){
        first_tri[0] = Pos[triangles[it]];
        first_tri[1] = Pos[triangles[it+1]];
        first_tri[2] = Pos[triangles[it+2]];
        //
        for (int jt = it+3; jt < 3*N; jt=jt+3){
        	second_tri[0] = Pos[triangles[jt]];
        	second_tri[1] = Pos[triangles[jt+1]];
        	second_tri[2] = Pos[triangles[jt+2]];
        }
        //
        if (tri_tri_intersection(first_tri,second_tri) == 1){
        	cout << it << "\t" << jt << endl;
        }
    }
}
#endif