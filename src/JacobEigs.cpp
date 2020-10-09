/* Question: Given a function f(x), how would you calculate 'n'
eigen values of the derivative of this function. 
Things to do :
  1) We need to decide the folder for this file. 
  2) class JdotX should have private variables instead of public ones.
     use a function to set the values of private variables.
  3) Find out cleverer way for passing a structure or not. Right now we
      are just overloading it, but is there a way not to do that? */
/***********************************************************/
#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
#include "JacobEigs.h"
#include "misc.h"
/*----------------------------*/
using namespace Spectra;
// using namespace Eigen;
// using namespace std; 
/***********************************************************/
// template <typename T>
// class JdotX{
// public:
//     int rows() { return Xstar -> rows(); }
//     int cols() { return Xstar -> cols(); }
//     /* Right now they are public variable but we should define them as private variable
//        and use function to define them. */
//     Map<VectorXd> *Xstar;
//     function<void(double*, T*)> f_display;
//     T* aTV;     // Address of the variable TV. TV stand for varible with T type.
//     /*----------------*/
//     void perform_op(const double *x_in, double *y_out){
//       Map<VectorXd> rhs(x_in,rows());   // rhs due to same terminology as Eigen.
//       VectorXd Xpos = Xstar+eps*rhs;
//       VectorXd Xneg = Xstar-eps*rhs;

//       VectorXd temp =  my_func(Xpos) - my_func(Xneg);
//       temp = temp/(2*eps);
//       yout = temp.data();
//     }
// private:
//   void my_func(VectorXd *vec){ 
//     double *vecdoub = (*vec).data();
//     f_display(vecdoub,aTV);
//     new (vec) Map<VectorXd>(vecdoub,rows());
//   }
// };
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
    /*----------------*/
    void perform_op(const double *x_in, double *y_out){
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
      // VectorXd temp =  (Xpos - Xneg)/(2*eps);
      // y_out = temp.data();
    }
    /*----------------*/
private:
  double eps=1.e-4;
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
//   eigs.compute();
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
// template<int const SelectionRule, int neigs, int ndim>
// bool jacob_eigvalvec(double eigval[], double eigvec[][ndim], void func(double*), double Xini[]){
//   JdotX op;
//   Map<VectorXd> xstar_op(Xini,ndim);
//   op.Xstar = &xstar_op;
//   op.f_display = func;
  
//   GenEigsSolver< double, SelectionRule, JdotX > eigs(&op, neigs, 2*neigs+1);
//   eigs.init();
//   eigs.compute();
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
// template<int const SelectionRule, int neigs, int ndim, typename T>
// bool jacob_eigval(double eigval[], void func(double*, T*), double Xini[], T *TV){
//   JdotX<T> op;
//   Map<VectorXd>xstar_op(Xini,ndim);
//   op.Xstar = &xstar_op;
//   op.f_display = func;
//   op.aTV = TV;

//   GenEigsSolver< double, SelectionRule, JdotX<T> > eigs(&op, neigs, 2*neigs+1);
//   eigs.init();
//   eigs.compute();
//   if(eigs.info() == SUCCESSFUL){
//       return 1;
//       VectorXd eigval_vec = eigs.eigenvalues();
//       eigval = eigval_vec.data();
//   }else{
//     return 0;
//   }
// }
/*--------------------------------------------------------------------------*/
// Wrapper function for the function which does not have structure as the second argument.
// The only argument is the array.
// template<int neigs, int ndim>
bool jacob_eigval(double eigval[], void func(double*), double Xini[]){
  int neigs = 3;
  int ndim = 8;
  JdotX op;
  // Map<VectorXd> xstar_op(Xini,ndim);
  op.Xstar = Xini;
  op.f_display = func;
  op.nn = ndim;

  GenEigsSolver< double, LARGEST_REAL, JdotX > eigs(&op, neigs, 2*neigs+1);
  eigs.init();
  eigs.compute();
  if(eigs.info() == SUCCESSFUL){
      Eigen::VectorXcd eigval_vec = eigs.eigenvalues();
      cout << eigval_vec << endl;
      // eigval = eigval_vec.data();
      return 1;
  }else{
    return 0;
  }
}
/*--------------------------------------------------------------------------*/