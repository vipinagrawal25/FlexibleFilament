#ifndef FILE_ODEN_SEEN
#define FILE_ODEN_SEEN
/*--------------------------------------------------*/
/* We solve for NN number of ordinary differential equations. 
Each of these equations can run independent of every other one
except when we need to average over them. Each of these equations
can have a stochastic component, and each must have a realisation of 
the noise independent of the other equations. 
The number of variables in each ode is pp. For example, we are solving
NN copies of the one dimensional harmonic oscillator then pp = 2 
(x and p ) */ 
#define NN 16
#define pp  2
#define ndim NN*pp
void alloc_odeN( double **PSI, double **psi);
void free_odeN( double **PSI, double **psi );
void H2D(double psi[], double PSI[], int Nsize );
void D2H(double PSI[], double psi[], int Nsize );
void iniconf(  double PSI[], double psi[]);
#endif /* !ODEN_SEEN */
