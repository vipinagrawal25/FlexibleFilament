#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN
/*--------------------------------------------------*/
struct MPARAM {
  double OM;
};
const int size_MPARAM = sizeof( double ) + 0*sizeof( int );
extern struct MPARAM host_param;
extern struct MPARAM *dev_param;
void set_param( void ) ;
__device__ void eval_rhs( double dpsi[], double psi[], int kelement, double tau,
                          MPARAM *dev_param  );
__host__ void initial_configuration( double PSI[] );
void write_param( void );
#endif /* !MODEL_SEEN */
