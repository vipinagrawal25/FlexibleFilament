#ifndef FILE_INPUT_SEEN
#define FILE_INPUT_SEEN
/* ----------------------------------------*/
unsigned int const Nensemble=1;
double dt = pow(10,-5);
double TMAX = 1.;
int idiag = 10;
double TotalFiles = 50;
double tdiag=TMAX/TotalFiles;	/*This decides after how much time interval the file should be saved.*/
int const wDataMeth = 1;	 //Writind data methos -> 1 for saving each snap in different file
int const rDataMeth = 1;						//	  2 for PSI in one file and velocity in another file.
double time = 0;
/* ----------------------------------------*/
#endif /* !FILE_INPUT_SEEN */