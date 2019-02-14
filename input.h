#ifndef FILE_INPUT_SEEN
#define FILE_INPUT_SEEN
/* ----------------------------------------*/
unsigned int const Nensemble=1;
double dt = pow(10,-4);
double TMAX=20/1;
int idiag=10;
double TotalFiles = 2000;
double tdiag=TMAX/TotalFiles;	/*This decides after how much time interval the file should be saved.*/
/* ----------------------------------------*/
#endif /* !FILE_INPUT_SEEN */
