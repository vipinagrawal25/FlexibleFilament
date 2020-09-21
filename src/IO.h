void wData(ofstream *fptr, double y[], double time) __attribute__((weak));
void wData(double fy[][ndim], double vel_all[][ndim], MV *aMM, int row) __attribute__((weak));
void wData(ofstream *fptr, ofstream *fptr_vel, double *y, double *vel, double time)__attribute__((weak));