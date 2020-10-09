// 1) Use template to pass a type of structure. Let's say that the function GG(y,MV)
// typename of second argument can be passed using typenamen and template.

// This file implements Matrix-free Newton Krylov method of findign zeros using GMRES. 
// We use predefined Eigen library and define a wrapper class to use GMRES in a matrix-free way.

// We should make templates like JacobEig function.
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "misc.h"
#include "NewtonKrylov.h"
#include <cmath>
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
      VectorXd Xpos = *lhs.my_vec+eps*rhs;
      VectorXd Xneg = *lhs.my_vec-eps*rhs;

      lhs.my_func(&Xpos);
      lhs.my_func(&Xneg);
      // cout << "Xpos new = " << Xpos << endl;
      // cout << "Xneg new= " << Xneg << endl; 
      VectorXd temp =  Xpos - Xneg;
      temp = temp/(2*eps);
      // cout << "Jdotdu norm = " << temp.norm() << endl;
      cout << "Jdotdu = " << temp << endl;
      dst.noalias() =  temp;
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
void newton_krylov(void func(double*), double Xini[], int ndim,
                  int Maxtry, double tol, double eps_temp){
	// input for Xstar is the initial guess.
  eps=eps_temp;
  MatrixReplacement AA;
  AA.f_display = func;
  double Err=1;
  VectorXd Xstar = Map<VectorXd>( Xini, ndim, 1 );
  // Map<VectorXd> Xstar(Xini,ndim);           // Converting double to VectorXd
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
    // Err = deltaX.norm()/Xstar.norm();
    // Err = max(deltaX.norm()/ndim, deltaX.maxCoeff());
    Err = (bb).norm()/Xstar.norm();
  }
  Xini = Xstar.data();
} 
/*-----------------------------------------------------------------------------*/
void newton_krylov(void func(double*), double Xini[], double gx[], int ndim,
                  int Maxtry, double tol, double eps_temp){
  MatrixReplacement AA;
  AA.f_display = func;
  eps=eps_temp;
  double Err=1;

  VectorXd Xstar = Map<VectorXd>( Xini, ndim,1 );
  VectorXd deltaX(ndim),bb(ndim);
  double *bbdoub;
  int itry=0;
  while(Err>tol){
    itry++;
    cout << "#Starting NewtonKrylov iteration: " << itry << endl;
    Xstar = Xstar+deltaX;
    AA.my_vec = &Xstar;
    GMRES<MatrixReplacement, IdentityPreconditioner> gmres;
    gmres.compute(AA);
    bb = Xstar;
    double *bbdoub = bb.data();
    func(bbdoub);
    bb = Map<VectorXd> (bbdoub,ndim,1);
    cout << "bb Norm = " << bb.norm() << endl;
    deltaX = gmres.solve(-bb);
    // Err = deltaX.norm()/Xstar.norm();
    // Err = max(deltaX.norm()/ndim, deltaX.maxCoeff());
    Err = (bb).norm()/Xstar.norm();
    cout << "#Finished NewtonKrylov iteration: " << itry << endl;
    cout << "#Error = " << Err << endl;
  }
  // Xini = Xstar.data();
  for (int idim = 0; idim < ndim; ++idim){
    gx[idim] = bb(idim);
    Xini[idim] = Xstar(idim);
  }
  // fy = bb.data();
}
/*-----------------------------------------------------------------------------*/