#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN
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
  double HH ;		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
  double KK;
  int qdiag;
  int bcb;
  int bct;
  int global_drag;
  int iext_force;
  int floc;
  int iext_flow ;
 };
const int size_MPARAM = 14*sizeof( double ) + 7*sizeof( int );
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
*/
extern struct MPARAM host_param;
extern struct MPARAM *dev_param;
extern double *DIAG;
extern double *dev_diag;
extern int size_diag;
/* The two following structures are defined in cuda.h file */
extern struct CRASH BUG;
extern struct CRASH *dev_bug;
/*------------------------------------------------------------------*/
__host__ void set_param( void ) ;
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          MPARAM *dev_param, double *diag, CRASH *crash  );
__host__ void initial_configuration( double PSI[] );
/*------------------------------------------------------------------*/
#endif /* !MODEL_SEEN */
