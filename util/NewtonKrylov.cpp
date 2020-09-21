// This file implements Matrix-free Newton Krylov method of findign zeros using GMRES. 
// We use predefined Eigen library and define a wrapper class to use GMRES in a matrix-free way.
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <../src/map_dyn.h>
#include <../src/constant.h>

using namespace std;
using namespace Eigen;

// Find X* such that F(X*) = 0
template <typename F>
void newton_krylov(VectorXd Xstar,  void func(),  ){
	
}