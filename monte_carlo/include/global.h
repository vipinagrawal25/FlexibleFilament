// this is a global file
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#define pi 3.14159265358979
#define R_del 0.05

typedef struct{
    double  x,y,z;
}POSITION;  //not included the celid and particle id which i shall do in the cell linked list part
//
typedef struct{
    int list_ss[200];
    int cnt_ss;
}Neighbours;
//
typedef struct{
    double dfac;
    int mc_iter;
    double kBT;
    double delta; // increment of position
    char *metric;
}MCpara;
//
typedef struct{
    double coef_bend;
    double coef_str;
    double radius;
    double sigma, epsilon;
    double pos_bot_wall;
    double coef_vol;
    double av_bond_len;
    int N;
    int num_triangles;
    int num_nbr;
}MBRANE_para;
//
typedef struct{
    int i1, i2;
}int2;
//
typedef struct{
    int size_node_nbr;
    int *cmlist;
    int *node_nbr_list;
    int2 *bond_nbr_list;
}MESH;

//
typedef struct{
    int N;
    double sigma;
    double epsilon;
    double len;
    double r_cut;
}LJpara;

#endif
