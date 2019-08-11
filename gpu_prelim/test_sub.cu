#include <stdio.h>
#include "sub.h"
  void allocA ( double **PSI, double **dev_psi) {
  double *temp;
  *PSI = (double *) malloc( ndim*sizeof( double ) );
  cudaMalloc( (void**)&temp, ndim*sizeof( double ) );
  *dev_psi = temp;
}
/*--------------------------------*/
void setA ( double A[], double dev_a[] ){
  for(int i = 0; i<ndim; i++){
    A[i] = (double) i ; 
}
  cudaMemcpy( dev_a,  A, ndim * sizeof(double),
              cudaMemcpyHostToDevice);
}
/*------------------------------------*/
void wA ( double A[] ){
  for( int i = 0; i<ndim; i++){
    printf ( "A[%d]=%f\t  ", i, A[i] );
  }
  printf ("\n" );
}
/*------------------------------------*/
int main ( void ) {
  double *A, *dev_a ;
//  double *dev_b;
  allocA ( &A, &dev_a );
  setA ( A, dev_a ) ;
  printf( " before \n " );
  wA ( A ) ;
 // dev_b = preAdd( ) ;
  preAdd( ) ;
  //AddB( A, dev_a, dev_b ) ;
  AddB( A, dev_a  ) ;
  printf( " after \n " );
  wA ( A ) ;
}
