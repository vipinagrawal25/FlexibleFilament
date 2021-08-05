#ifndef FILE_INPUT_SEEN
#define FILE_INPUT_SEEN
#include <cmath>
/* ----------------------------------------*/
unsigned int const Nensemble=1;
double const start_dt = 0.0001;
double const TMAX=960;
double const TotalFiles = 48000;
double const tdiag=TMAX/TotalFiles;					/*This decides after how much time interval the file should be saved.*/
<<<<<<< HEAD
double const start_time = 10*2*M_PI/3;
=======
double const start_time = 0;
>>>>>>> tracer
int const wDataMeth = 1;
int const rDataMeth = 1;
/* ----------------------------------------*/
#endif /* !FILE_INPUT_SEEN */