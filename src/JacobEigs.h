#ifndef FILE_JACOBEIGS_SEEN
#define FILE_JACOBEIGS_SEEN
/*****************************************************************/
// template<int const SelectionRule, int neigs, int ndim, typename T>
// bool jacob_eigvalvec(double eigval[], double eigvec[][ndim], void func(double*, T*),
//                      double Xini[], T *TV);

// template<int const SelectionRule, int neigs, int ndim>
// bool jacob_eigvalvec(double eigval[], double eigvec[][ndim], void func(double*), double Xini[]);

// template<int const SelectionRule, int neigs, int ndim, typename T>
// bool jacob_eigval(double eigval[], void func(double*, T*), double Xini[], T *TV);

// template<int neigs, int ndim>
bool jacob_eigval(double eigval[], void func(double*), double Xini[]);
/*****************************************************************/
#endif