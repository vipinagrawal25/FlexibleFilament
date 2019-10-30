#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN
#include "cuda.h"
/*--------------------------------------------------*/
/* We solve for a chain, a one dimensional system. 
Each node may be connected to every other node. 
The number of nodes is NN. 
The number of degrees of freedom at each node is pp */ 
#define NN 100
#define pp  3
#define ndim NN*pp
#define TimeScheme "rnkf45"
/*--------------------------------------------------*/
struct MPARAM {
  double height ;	// height of the box we are doing simulations in.
  double aa; 	// distance between two nodes.
  double Dbyell; // diameter/length of the filament.
  double dd ;	/* r/l ratio for the rod has been kept constant. 
                   It should be noted that the particles would also have same diameter. */
  double viscosity ;				
  double  Z0;	  // If we want the bottom point of the rod to be fixed.
  double Famp ; // Force Value on the ends
  // Sigma is a dimensionless number, which is described as frequency parameter.
  double sigma ;					
  double ShearRate ;
  double omega ; // frequency of external force.
  double  factorAA ; 
  double AA ;
  double HH ;	// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
  double KK;
  int qdiag;
  int bcb;
  int bct;
  int global_drag;
  int iext_force;
  int floc;
  int iext_flow ;
  int iniconf;
};
const int size_MPARAM = 14*sizeof( double ) + 8*sizeof( int );
const int size_psi = ndim*sizeof(double);
/* boundary condition :  
   0 => clamped
   1=> free
   2 => hinged
*/
/* iext_force : implements external force on the filament
   periodic forcing at  position floc */
/* external flow:
   0 => no external flow.
  1 => time-dependent shear U = ( ShearRate*z, 0, 0 ) * square_wave(omega*time) 
  2 => time-INdependent shear U = ( ShearRate*z, 0, 0 ) 
/*------------------------------------------------------------------*/
void alloc_chain( double **PSI, double **psi);
void free_chain( double **PSI, double **psi );
void H2D(double psi[], double PSI[], int Nsize );
void D2H(double PSI[], double psi[], int Nsize );
void set_param( MPARAM *PARAM, MPARAM **dev_param ) ;
void write_param( MPARAM *PARAM, char *fname );
int pre_diag( double **DIAG , double **dev_diag, MPARAM PARAM );
__device__ void model_rhs( double dpsi[], double psi[], int kelement, double tau,
                           MPARAM *param, double *diag, CRASH *bug, int ldiag  );
void initial_configuration( double PSI[], MPARAM PARAM );
void wPSI ( double PSI[], double tau ); 
void wDIAG( double DIAG[], double tau, MPARAM PARAM );
bool check_param(MPARAM PARAM);
/*------------------------------------------------------------------*/
#endif /* !MODEL_SEEN */
