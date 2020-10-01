// This file implements Matrix-free Newton Krylov method of findign zeros using GMRES. 
// We use predefined Eigen library and define a wrapper class to use GMRES in a matrix-free way.
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "misc.h"
#include "NewtonKrylov.h"
/*---------------------------------------------------------------------------- */
double eps;
using namespace std;
using namespace Eigen;
/*---------------------------------------------------------------------------- */
class MatrixReplacement;
namespace Eigen {
namespace internal {
  template<>
  struct traits<MatrixReplacement> :  public internal::traits< SparseMatrix<double> >
  {};
}
}
class MatrixReplacement : public EigenBase<MatrixReplacement> {
public:
  VectorXd *my_vec;
  function<void(double*)>f_display;
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum{
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
  };
 
  Index rows() const { return my_vec->rows(); }
  Index cols() const { return my_vec->cols(); }
 
  template<typename Rhs>
  Product<MatrixReplacement,Rhs,AliasFreeProduct> operator*(const MatrixBase<Rhs>& x) const 
  {
    return Product<MatrixReplacement,Rhs,AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  // template<typename T>
  // void assign_var(T vec, void func(double*) ){
  //   mp_vec = vec;
  //   f_display = func;
  // }
  const void my_func(VectorXd *vec)const{ 
    double *vecdoub = (*vec).data();
    f_display(vecdoub);
    new (vec) Map<VectorXd>(vecdoub,rows());
  }
  // void my_vec(VectorXd *XX) const { XX =  mp_vec; }
private:
};
/*----------------------------------------------------------------------------*/
// Implementation of scaleAndAddTo function
namespace Eigen {
namespace internal {
  template<typename Rhs>
  struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> 
  // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement,Rhs, generic_product_impl<MatrixReplacement,Rhs> >
  {
    typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha){
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      // assert(alpha==Scalar(1) && "scaling is not implemented");
      // EIGEN_ONLY_USED_FOR_DEBUG(alpha);
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // VectorXd Xpos(lhs.rows());
      // lhs.my_vec(&Xpos);
      // VectorXd *Xpos = lhs.my_vec;
      // VectorXd *Xneg = lhs.my_vec;
      // &Xneg = &Xpos;
      VectorXd Xpos = *lhs.my_vec+eps*rhs;
      VectorXd Xneg = *lhs.my_vec-eps*rhs;
      // cout << "Xpos = " << Xpos << endl;
      // cout << "Xneg = " << Xneg << endl; 
      lhs.my_func(&Xpos);
      lhs.my_func(&Xneg);
      // cout << "Xpos new = " << Xpos << endl;
      // cout << "Xneg new= " << Xneg << endl; 
      VectorXd temp =  Xpos - Xneg;
      dst.noalias() =  temp/(2*eps);
      // for(Index i=0; i<lhs.cols(); ++i)
      //   dst += rhs(i) * lhs.my_matrix().col(i);
    }
  };
}
}
/*-----------------------------------------------------------------------------*/
// newton_krylov function should accept variable number of arguments.
// It should also accept guess as a function or a vector.
// function should also accept guess() as an argument.
// Find X* such that F(X*) = 0
// aMM will contain every other parameter in a structure form, it should be passed to a function.
// template<int Maxtry>
void newton_krylov(void func(double*), double yini[], int ndim,
                  int Maxtry, double tol, double eps_temp){
	// input for Xstar is the initial guess.
  eps=eps_temp;
  MatrixReplacement AA;
  AA.f_display = func;
  double Err=1;
  VectorXd Xstar = Map<VectorXd>( yini, ndim, 1 );
  // Map<VectorXd> Xstar(yini,ndim);           // Converting double to VectorXd
  VectorXd deltaX(ndim),bb(ndim);
  double *bbdoub;
  while(Err>tol){
    Xstar = Xstar+deltaX;
    AA.my_vec = &Xstar;
    GMRES<MatrixReplacement, IdentityPreconditioner> gmres;
    gmres.compute(AA);
    bb = Xstar;
    double *bbdoub = bb.data();
    func(bbdoub);
    bb = Map<VectorXd>(bbdoub,ndim,1);
    deltaX = gmres.solve(-bb);
    Err = deltaX.norm()/Xstar.norm();
  }
  yini = Xstar.data();
} 
// /*--------------------------------------*/
// // Wrapper function
// template <double Maxtry, double tol, double eps, typename T, typename ... Args>
// double *newton_krylov(double Xstar, int ndim, double *func(), double *guess(), Args ... args){
//   // If we get guess from a function.
// }
// --------------------------------------
// // Wrapper function
// template <double Maxtry, double tol, double eps, typename T, typename ... Args>
// double *newton_krylov(double Xstar, int ndim, double *func(), void Gradient(), Args ... args){
// 	// Xstar is the initial guess.
// 	// In case gradient is retuned directly from a function. 
// }
/*-----------------------------------------------------------------------------*/
// template<int Maxtry>
void newton_krylov(void func(double*), double yini[], double gy[], int ndim,
                  int Maxtry, double tol, double eps_temp){
  MatrixReplacement AA;
  AA.f_display = func;
  eps=eps_temp;
  double Err=1;

  VectorXd Xstar = Map<VectorXd>( yini, ndim,1 );
  VectorXd deltaX(ndim),bb(ndim);
  double *bbdoub;
  while(Err>tol){
    Xstar = Xstar+deltaX;
    AA.my_vec = &Xstar;
    GMRES<MatrixReplacement, IdentityPreconditioner> gmres;
    gmres.compute(AA);
    bb = Xstar;
    double *bbdoub = bb.data();
    func(bbdoub);
    bb = Map<VectorXd> (bbdoub,ndim,1);
    deltaX = gmres.solve(-bb);
    // cout << deltaX << endl;
    Err = deltaX.norm()/Xstar.norm();
    // cout << Err << endl;
  }
  // yini = Xstar.data();
  for (int idim = 0; idim < ndim; ++idim){
    gy[idim] = bb(idim);
    yini[idim] = Xstar(idim);
  }
  // fy = bb.data();
}
/*-----------------------------------------------------------------------------*/
void test(void func(double*), double yini[], double fy[], int ndim,
          int Maxtry, double tol, double eps_temp){
  VectorXd Xstar = Map<VectorXd>(yini, ndim,1);
  // VectorXd Xstar(ndim);
  // Xstar.data() = yini;
  for (int idim = 0; idim < ndim-1; ++idim){
    Xstar(idim) = Xstar(idim)+5;
  }
}
/*-----------------------------------------------------------------------------*/