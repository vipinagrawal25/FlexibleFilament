#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
#include "misc.h"
#include <string>
#include <complex>
using namespace std;
/*----------------------------*/
#ifndef FILE_JACOBEIGS_SEEN
#define FILE_JACOBEIGS_SEEN
/*----------------------------*/
using namespace Spectra;
using namespace Eigen;
using namespace std;
// /*****************************************************************/
// // template<int const SelectionRule, int neigs, int ndim, typename T>
// // bool jacob_eigvalvec(double eigval[], double eigvec[][ndim], void func(double*, T*),
// //                      double Xini[], T *TV);

// template<int const SelectionRule>
// bool jacob_eigvalvec(VectorXcd *eigval, MatrixXcd *eigvec, int neigs, int ndim,
//                      void func(double*), double Xini[], double eps_temp=1.e-4);

// template<int const SelectionRule, typename T>
// bool jacob_eigval(VectorXcd *eigval, int neigs, int ndim, void func(double*, T*), 
//                   double Xini[], T *TV, double eps_temp=1.e-4);

// template<int const selectionRule>
// bool jacob_eigval(double eigval[], int neigs, int ndim, void func(double*), 
//                   double Xini[],double eps_temp=1.e-4);
/*****************************************************************/
/* Question: Given a function f(x), how would you calculate 'n'
eigen values of the derivative of this function. 
Things to do :
  1) We need to decide the folder for this file. 
  2) class JdotX should have private variables instead of public ones.
     use a function to set the values of private variables.
  3) Find out cleverer way for passing a structure or not. Right now we
      are just overloading it, but is there a way not to do that? */
/*--------------------------------------------------------------------------*/
template <typename T>
class T_JdotX{
public:
    int rows() { return nn; }
    int cols() { return nn; }
    /* Right now they are public variable but we should define them as private variable
       and use function to define them. */
    function<void(double*, T*)> f_display;
    int nn;
    T* aTV;     // Address of the variable TV. TV stand for varible with T type.
    double *Xstar;
    double eps;
    /*----------------*/
    void perform_op(const double *x_in, double *y_out){
      double Xpos[rows()],Xneg[rows()];
      for (int i = 0; i < rows(); ++i){
        Xpos[i] = Xstar[i] + eps*x_in[i];
        Xneg[i] = Xstar[i] - eps*x_in[i];
      }
      f_display(Xpos,aTV);
      f_display(Xneg,aTV);
      for (int i = 0; i < rows(); ++i){
        y_out[i] = (Xpos[i] - Xneg[i])/(2*eps);
      }
    }
private:
};
/*--------------------------------------------------------------------------*/
// Wrapper to JdotX class without template. 
class JdotX{
public:
    int rows() { return nn; }
    int cols() { return nn; }
    /* Right now they are public variable but we should define them as private variable
       and use function to define them. */
    // Map<VectorXd> *Xstar;
    function<void(double*)> f_display;
    int nn;
    double *Xstar;
    double eps;
    /*----------------*/
    void perform_op(const double *x_in, double *y_out){
      cout << "Arnoldi iteration " << iter << " for eigenvalues." << endl;
      double Xpos[rows()],Xneg[rows()];
      // Map<VectorXd> rhs(x_in,rows());   // rhs due to same terminology as Eigen.
      // VectorXd Xpos = Xstar+eps*rhs;
      // VectorXd Xneg = Xstar-eps*rhs;
      for (int i = 0; i < rows(); ++i){
        Xpos[i] = Xstar[i] + eps*x_in[i];
        Xneg[i] = Xstar[i] - eps*x_in[i];
      }
      f_display(Xpos);
      f_display(Xneg);
      for (int i = 0; i < rows(); ++i){
        y_out[i] = (Xpos[i] - Xneg[i])/(2*eps);
      }
      iter++;
      // VectorXd temp =  (Xpos - Xneg)/(2*eps);
      // y_out = temp.data();
    }
    /*----------------*/
private:
  int iter=1;
  // void my_func(VectorXd *vec){ 
  //   double *vecdoub = (*vec).data();
  //   f_display(vecdoub);
  //   new (vec) Map<VectorXd>(vecdoub,rows());
  // }
};
/*--------------------------------------------------------------------------*/
// /* The function expects a function with first argument to be double pointer, 
//   and the second arguement can be a structure. 
//   Typename of structure needs to be passed as template arguments. */
// template<int const SelectionRule, int neigs, int ndim, typename T>
// bool jacob_eigvalvec(double eigval[], double eigvec[][ndim], void func(double*, T*),
//                      double Xini[], T *TV){
//   JdotX<T> op;
//   Map<VectorXd>xstar_op(Xini,ndim);
//   op.Xstar = &xstar_op;
//   op.f_display = func;
//   op.aTV = TV;

//   GenEigsSolver< double, SelectionRule, JdotX<T> > eigs(&op, neigs, 2*neigs+1);
//   eigs.init();
//   eigs.compute(ndim,1.e-3);
//   if(eigs.info() == SUCCESSFUL){
//       return 1;
//       VectorXd eigval_vec = eigs.eigenvalues();
//       eigval = eigval_vec.data();
//       MatrixXd eigvec_mat = eigs.eigenvectors();
//       eigvec = eigvec_mat.data();
//   }else{
//     return 0;
//   }
// }
// /*--------------------------------------------------------------------------*/
template<int const SelectionRule>
bool jacob_eigvalvec(VectorXcd *eigval, MatrixXcd *eigvec, int neigs, int ndim,void func(double*), 
                     double Xini[], double eps_temp=1.e-4){
  JdotX op;
  op.Xstar = Xini;
  op.f_display = func;
  op.nn = ndim;
  op.eps = eps_temp;

  int ncv = min(2*neigs+1,ndim);
  GenEigsSolver< double, SelectionRule, JdotX > eigs(&op, neigs, ncv);
  eigs.init();
  eigs.compute(ndim,1.e-3);
  if(eigs.info() == SUCCESSFUL){
    *eigval = eigs.eigenvalues();
    *eigvec = eigs.eigenvectors();
    return 1;
  }else{
    return 0;
  }

}
/*--------------------------------------------------------------------------*/
template<int const selectionRule, typename T>
bool jacob_eigval(VectorXcd *eigval, int neigs, int ndim, void func(double*, T*),
                  double Xini[], T *TV, double eps_temp=1.e-4){
  T_JdotX<T> op;
  op.Xstar = Xini;
  op.f_display = func;
  op.nn = ndim;
  op.aTV = TV;
  op.eps = eps_temp;

  int ncv = min(2*neigs+1,ndim);
  GenEigsSolver< double, selectionRule, T_JdotX<T> > eigs(&op, neigs, ncv);
  eigs.init();
  eigs.compute(ndim,1.e-3);

  if(eigs.info() == SUCCESSFUL){
    *eigval = eigs.eigenvalues();
    cout << *eigval;
    return 1;
  }else{
    return 0;
  }
}
/*--------------------------------------------------------------------------*/
// Wrapper function for the function which does not have structure as the second argument.
// The only argument is the array.
template<int const selectionRule>
bool jacob_eigval(VectorXcd *eigval, int neigs, int ndim, void func(double*),
                  double Xini[], double eps_temp=1.e-4){
  // TempV tempV;
  // void func_temp = [func](double Xini[],TempV * tempV) {func(Xini)};
  // jacob_eigval<selectionRule>(eigval,neigs,ndim,func,Xini,&tempV,eps_temp);
  cout << selectionRule << endl;
  JdotX op;
  op.Xstar = Xini;
  op.f_display = func;
  op.nn = ndim;
  op.eps=eps_temp;

  int ncv = min(2*neigs+1,ndim);
  GenEigsSolver< double, selectionRule, JdotX > eigs(&op, neigs, ncv);
  eigs.init();
  int nconv = eigs.compute(ndim);
  if(eigs.info() == SUCCESSFUL){
    *eigval = eigs.eigenvalues();
    cout << *eigval;
    return 1;
  }else{
    return 0;
  }
}
/*****************************************************************/
#endif