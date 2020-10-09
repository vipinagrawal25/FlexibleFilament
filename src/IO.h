#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN
/*--------------------------------------------------*/
void wData(ofstream *fptr, double *y, double time=0.0) __attribute__((weak));
void wData(ofstream *fptr, ofstream *fptr_vel, double *y, double *vel, double time=0) __attribute__((weak));
// void rData(ifstream *fptr, double *y) __attribute__((weak));
void rData(double *y, string filename) __attribute__((weak));
// void wData(string fname, double y[][ndim], double time[], int nrow)__attribute__((weak));;
/*--------------------------------------------------*/
#endif
