#ifndef FILE_ESTRING_SEEN
#define FILE_ESTRING_SEEN
/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, double* add_SS, bool flag_kappa);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec3 X[], int i, int j, double *Distance, vec3 *rij);
/**************************/
using namespace std;
/* ----------------------------------------*/
unsigned int const Np= 40;				// total number of points on the rod.
unsigned int const pdim=3*Np; 
double const aa =  2.0/(double)(Np); 	// distance between two nodes.
// AA = AA/aa ;							// See the equation number 35 in DOI: 10.1063/1.1848511
double const AA = 0.001;		 // 0.001 is bending rigidity constant (AA = EI/aa), unit -> Pa.m^4
double const dd = 0.1*aa;
double const viscosity = 1;
double const HH = AA*(1./(aa*aa));		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
// The below line is commented because of some error. So I included it in last formulae itself.
// HH =  HH/aa; 			// This is done to model the rod perfectly (Also see: https://goo.gl/nbELSx Week 27/2018 DS)
// double const HH = 0.;
double const OneByGamma=1.;				// Mobility is just the inverse of Gamma.
double const Z0=0.;
double const FFZ0 = 0.1;				// Force Value on the ends
double const omega = 0.02;
char const UseRP = 'N';
/* ----------------------------------------*/
void iniconf(double *y);
void diagnos(double time, int itn, double y[]);
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */
