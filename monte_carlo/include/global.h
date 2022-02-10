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
#include <iostream>
#include <cstdlib>

using namespace std;

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
    int one_mc_iter;
    int tot_mc_iter;
    int dump_skip;
    double kBT;
    double delta; // increment of position
    bool is_restart;
    char *metric;
}MCpara;
//
typedef struct{
    double coef_bend;  //coefficient bending
    double coef_str;  //coefficient stretching
    double coef_vol_expansion;   //coefficient of volume expansion
    double radius;  // radius of ball
    double sigma, epsilon; // sigma and epsilon for the bottom attractive wall
    double pos_bot_wall;  // position of the bottom attractive wall
    double av_bond_len; // average length of the bond
    double *tot_energy;
    double *volume; // these are updated after each monte carlo 
    int N;   // number of particles in mesh
    int num_triangles;  //number of triangles 
    int num_nbr; // sum  neighbours of all particles
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
int N; // number of points defining in afm tip curve
double extent[4]; // extent in x- x+ y- y+
POSITION *tip_curve;
double tip_pos_z; // position of tip in z
double tip_rad; // radius of the tip
double sigma, epsilon;
}AFM_para;
//
typedef struct{
    int N;
    double sigma;
    double epsilon;
    double len;
    double r_cut;
}LJpara;
#endif
