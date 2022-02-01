#include<stdio.h>
void visit_vtk_io(double *points, 
        int *triangles, 
        int Np, char filename[], 
        char dataname[]){

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

void visit_vtk_io_points_data(double *data, 
        int Np, char filename[], 
        char dataname[]){

    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, mydata\n ASCII \n DATASET POLYDATA \n";
    int num_triangles, i;
    FILE *fid;

    num_triangles = 2*Np - 4;
    fid = fopen(filename, "a");
    fprintf(fid, "%s %d %s \n", "POINT_DATA", Np);
    fprintf(fid, "%s %s %s \n", "SCALARS", dataname, "float 1");
    fprintf(fid, "%s \n", "LOOKUP_TABLE default ");
    for(i = 0; i< Np; i=i+1){
        fprintf(fid, "%g \n", data[i]);
    }
    fclose(fid);
}

void visit_vtk_io_cell_data(double *data, 
        int Np, char filename[], 
        char dataname[]){

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


