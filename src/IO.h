#include "model.h"
#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN
/*--------------------------------------------------*/
void rData(double *y);
void wData(ofstream *fptr, double *y, double time=0.0, int Npt=Np, int ppt=pp) __attribute__((weak));
void wData(ofstream *fptr, ofstream *fptr_vel, double *y, double *vel, double time=0,
		   int Npt=Np, int ppt=pp) __attribute__((weak));
// void rData(ifstream *fptr, double *y) __attribute__((weak));
// void wData(string fname, double y[][ndim], double time[], int nrow)__attribute__((weak));;
/*--------------------------------------------------*/
#endif