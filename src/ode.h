// void euler(unsigned int ndim, double *y, double time,double dt);
// void rnkt2(unsigned int ndim, double *y, double time,double dt);
void rnkt4(unsigned int ndim, double *y, double *vel, double *time,double *dt, double* CurvSqr, 
			double* SS, double ldiagnos);
void rnkt4(double *y, double *add_time, double* add_dt);
void rnkf45(unsigned int ndim, double *y, double *vel, double *time, double *dt, double* CurvSqr, 
			double* SS, double ldiagnos);
void rnkf45(double *y, double *add_time, double* add_dt);
void DP54(unsigned int ndim, double *y, double *vel, double *time, double *dt, double* CurvSqr, 
		 double* SS, double ldiagnos);