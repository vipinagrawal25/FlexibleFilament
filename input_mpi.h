#ifndef FILE_INPUT_SEEN
#define FILE_INPUT_SEEN
/* ----------------------------------------*/
extern double TMAX;
extern int TotalFiles;
/* ----------------------------------------*/
unsigned int const Nensemble=1;
double const start_dt = 1.e-5;
double const tdiag=TMAX/TotalFiles;	/*This decides after how much time interval the file should be saved.*/
int const wDataMeth = 1;	 // Writind data methos -> 1 for saving each snap in different file
int const rDataMeth = 1;						//	  2 for PSI in one file and velocity in another file.
double const start_time = 0;
/* ----------------------------------------*/
#endif /* !FILE_INPUT_SEEN */