#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN
/*--------------------------------------------------*/
struct MPARAM {
  double height ;	// height of the box we are doing simulations in.
  double aa; 	// distance between two nodes.
  double Dbyell // diameter/length of the filament.
  double dd ;	/* r/l ratio for the rod has been kept constant. 
                   It should be noted that the particles would also have same diameter. */
  double viscosity ;				
  double  Z0;	  // If we want the bottom point of the rod to be fixed.
  double FFZ0 ; // Force Value on the ends
// Sigma is a dimensionless number, which is described as frequency parameter.
  double sigma ;					
  double ShearRate ;
  double omega ;
  double  factorAA ; 
  double AA ;
  double HH ;		// Follow: bit.ly/2r23lmA unit -> Pa.m^4/m^2 -> Pa.m^2
  double KK;
  int qdiag;
};
const int size_MPARAM = 13*sizeof( double ) + sizeof( int );
extern struct MPARAM host_param;
extern struct MPARAM *dev_param;
extern double *DIAG;
extern double *dev_diag;
void set_param( void ) ;
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          MPARAM *dev_param, double *diag  );
__host__ void initial_configuration( double PSI[] );
void write_param( void );
#endif /* !MODEL_SEEN */
