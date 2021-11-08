void euler_tr(unsigned int ndim, double *y, double *y_tr, double* vel_tr, double time, double dt,double *EForceArr);
void rnkt4(unsigned int ndim, double *y, double *vel, double *time,double *dt, double ldiagnos);
void rnkt4(double *y, double *add_time, double* add_dt);
void rnkf45(unsigned int ndim, double *y, double *vel, double *time, double *dt, double *EForceArr, double ldiagnos);
void rnkf45(double *y, double *add_time, double* add_dt);
void DP54(unsigned int ndim, double *y, double *vel, double *add_time, double* add_dt, double ldiagnos);
void eval_rhs(double rhs[], double y[], double time);