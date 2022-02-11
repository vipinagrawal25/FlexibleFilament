#include<math.h>
#include <iostream>
#include <string.h>
using namespace std;
/*******************************************/
void visit_vtk_io(double *points, 
        int *triangles, 
        int Np, char filename[], string dataname){

    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, mydata\n ASCII \n DATASET POLYDATA \n";
    int num_triangles, i;
    FILE *fid;

    num_triangles = 2*Np - 4;
    fid = fopen(filename, "w");
    fprintf(fid, "%s", first_headers);
    fprintf(fid, "%s %d %s", "POINTS", Np, "float \n");
    for(i = 0; i< 3*Np; i=i+3){
        fprintf(fid, "%g %g %g \n", points[i], points[i+1], points[i+2]);
    }
    fprintf(fid, "%s %d  %d %s", "POLYGONS", num_triangles, 4*num_triangles, " \n");
    for(i = 0; i< 3*num_triangles; i=i+3){
        fprintf(fid, "%d %d %d %d\n", 3, triangles[i], triangles[i+1], triangles[i+2]);
    }
    fclose(fid);
}

void visit_vtk_io_point_data(bool *data, 
        int Np, char filename[], 
        string dataname){

    int num_triangles, i;
    FILE *fid;

    num_triangles = 2*Np - 4;
    fid = fopen(filename, "a");
    fprintf(fid, "%s %d \n", "POINT_DATA", Np);
    fprintf(fid, "%s %s %s \n", "SCALARS", dataname, "float 1");
    fprintf(fid, "%s \n", "LOOKUP_TABLE default ");
    for(i = 0; i< Np; i=i+1){
        if(data[i]){
            fprintf(fid, "%g \n", 1.0);
        } else {
            fprintf(fid, "%g \n", 0.0);
        }
    }
    fclose(fid);
}

void visit_vtk_io_cell_data(double *data, 
        int Np, char filename[], 
        string dataname){

    int num_triangles, i;
    FILE *fid;

    num_triangles = 2*Np - 4;
    fid = fopen(filename, "a");
    /* if(fid != NULL)printf("I am null"); */
    fprintf(fid, "%s %d  \n", "CELL_DATA", Np);
    fprintf(fid, "%s %s %s \n", "SCALARS", dataname, "float 1");
    fprintf(fid, "%s \n", "LOOKUP_TABLE default ");
    for(i = 0; i< Np; i=i+1){
        fprintf(fid, "%g \n", data[i]);
    }
    fclose(fid);
}

void visit_vtk_io_afm_tip(double *data, 
        int Np, char filename[]){

    int i, j, here, est, nth, n_est;
    double x, y;
    FILE *fid;
    int iext, idx;
    double dx;

    iext = (int) sqrt(Np);
    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, afmtip\n ASCII \n DATASET POLYDATA \n";

    fid = fopen(filename, "w");
    fprintf(fid, "%s", first_headers);
    fprintf(fid, "%s %d %s", "POINTS", Np, "float \n");

    for(i=0; i < 3*Np; i=i+3){
        fprintf(fid, "%g %g %g\n", data[i], data[i+1], data[i+2]);
    }

    int npoly = (iext-1)*(iext-1);
    fprintf(fid, "%s %d  %d %s", "POLYGONS", 4*npoly, 5*4*npoly, " \n");
    for(j=0; j < iext-1; j++){
        for(i=0; i < iext-1; i++){
            here = j*iext + i;
            est = (i+1)%iext + ((iext + j)%iext) *iext;
            nth = ((iext + i)%iext) + (j + 1 + iext)%iext *iext;
            n_est = ((i + 1)%iext) + ((j + 1 + iext)%iext)*iext;
            fprintf(fid, "%d %d %d %d %d\n", 4, here, est, n_est, nth);
        }
    }
    fclose(fid);
}