void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[],double SS[]);
// void euler(unsigned int ndim, double *y, double time,double dt);
// void rnkt2(unsigned int ndim, double *y, double time,double dt);
void rnkt4(unsigned int ndim, double *y, double *time,double *dt, double* CurvSqr, double* SS, double ldiagnos);
void rnkf45(unsigned int ndim, double *y, double *time, double *dt, double* CurvSqr, double* SS, double ldiagnos);
void DP54(unsigned int ndim, double *y, double *time, double *dt, double* CurvSqr, double* SS, double ldiagnos);
