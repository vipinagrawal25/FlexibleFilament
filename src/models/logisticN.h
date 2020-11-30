#ifndef FILE_LOGISTICN_SEEN
#define FILE_LOGISTICN_SEEN
/**************************/
using namespace std;
/* ----------------------------------------*/
unsigned int const Np=8;
unsigned int const pp= 1;
unsigned int const ndim=pp*Np;
extern double lambda[Np];
// double const *lambda = new double(Np);
double const omega=1;
/* ----------------------------------------*/
void iniconf(double *y); 	// The configuration number is defined for different
							// Initial configuration into the system.
void eval_rhs(double *y);
void definelambda();
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]);
extern string cc;
/* ----------------------------------------*/
#endif /* !FILE_ESTRING_SEEN */