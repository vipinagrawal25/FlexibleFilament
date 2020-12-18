#ifndef FILE_NEWTONKRYLOV_SEEN
#define FILE_NEWTONKRYLOV_SEEN
/*****************************************************************/
// template<int Maxtry=128>
bool newton_krylov(void func(double*), double yini[], int ndim,
                   double tol=1.e-3, double eps_temp=1.e-4,int Maxtry=128, bool verbose=1);

// template<int Maxtry=128>
bool newton_krylov(void func(double*), double yini[], double fy[], int ndim,
                   double tol=1.e-3, double eps_temp=1.e-4, int Maxtry=128, bool verbose=1);
/*****************************************************************/
#endif