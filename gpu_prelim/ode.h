// __global__ eval_rhs(double time,double *y,double *rhs);
//void rnkt4(unsigned int ndim, double *y, double time,double dt);
/* These are device pointers I need for the time-stepping code */
/* I always need one .. */
double *temp;
/* For higher-order time-steppers I need more. Here I go 
   up to 4th order (even adaptive) time-steppers. So I need .. */
double *k1,*k2,*k3,*k4, *k5;
__global__ void pre_euler( void );
__global__ void euler( double *psi, double time,double dt);
