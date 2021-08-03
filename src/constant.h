#ifndef FILE_CONSTANT_SEEN
#define FILE_CONSTANT_SEEN
//---------------------------//
#define TimeScheme "rnkf45"
#define time_tol 1.e-9
#define err_tol 1.e-2
// #define delta 1.e-2
#define tiny 1.e-15
#define SysType "continuous"
//--------------------------//
#endif
/*  TimeScheme: You can choose Euler, rnkt2, rnkt4, rnkf45(adaptive runge-kutta scheme), 
			   DP54 (Dormand-Prince method)
	time_tol -> 1 iteration of time has error of a time_tol.
	err_tol -> 	Mainly for map_dynamics. An orbit/fixed point has this error
	delta	-> 	delta == deltax,deltay etc. dx or dy to calculate Jacobian. 
				It should not be very small otherwise numerical errors get significant in the calculation of Jacobian.
	wDataMeth(Writing/Reading data type): 
				1 -> to save every snap in different file
                2 -> to save all the snaps in single file called PSI, VEL
	SysType	-> Is it a map or ODE? */