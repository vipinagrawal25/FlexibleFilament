#ifndef FILE_2Tens_SEEN
#define FILE_2Tens_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "2vec.h"
using namespace std;
/*-------------------------*/
class Tens2{
public:
  double xx,xy,yx,yy;
  Tens2();
  Tens2(int,int,int,int);
  Tens2(float,float,float,float);
  Tens2(double,double,double,double);
    // definining operators 
  Tens2 inline operator+(Tens2);
  Tens2 inline operator-(Tens2);
  Tens2 inline operator*(double);
  Tens2 inline operator/(double);
private:
};
Tens2::Tens2(){
  xx = 0.0;
  xy = 0.0;
  
  yx = 0.0;
  yy = 0.0;
}
Tens2::Tens2(int a, int b, int c, int d){
  xx = (double) a;
  xy = (double) b;
  
  yx = (double) c;
  yy = (double) d;
}
Tens2::Tens2(float a, float b, float c, float d){
  xx = (double) a;
  xy = (double) b;
  
  yx = (double) c;
  yy = (double) d;
}
Tens2::Tens2(double a, double b, double c, double d){
  xx =  a;
  xy =  b;

  yx =  c;
  yy =  d;
}
Tens2 Tens2::operator+(Tens2 param){
  Tens2 temp;
  temp.xx = xx+param.xx;
  temp.xy = xy+param.xy;

  temp.yx = yx+param.yx;
  temp.yy = yy+param.yy;

  return(temp);
}
Tens2 Tens2::operator-(Tens2 param){
  Tens2 temp;
  temp.xx = xx-param.xx;
  temp.xy = xy-param.xy;

  temp.yx = yx-param.yx;
  temp.yy = yy-param.yy;

  return(temp);
}
Tens2 Tens2::operator*(double param){
  Tens2 temp;
  temp.xx=param*xx;
  temp.xy=param*xy;

  temp.yx=param*yx;
  temp.yy=param*yy;

  return(temp);
}
Tens2 operator*(const double param, Tens2 a){
    // This function is just here because multiplication between double and a Tensor needs to be commutative.
    return(a*param);
}
Tens2 Tens2::operator/(double param){
  Tens2 temp;
  temp.xx=xx/param;
  temp.xy=xy/param;

  temp.yx=yx/param;
  temp.yy=yy/param;

  return(temp);
}
void PTens2(Tens2 a){
  cout<<a.xx<<"\t"<<a.xy<<"\n"<<a.yx<<"\t"<<a.yy<<"\n";
}
// This defines the kronecker delta function or a unit tensor matrix. dab = \delta_{\alpha \beta}
Tens2 dab(1.,0.,0.,1.);
/*---------------------------------------*/
#endif /* !FILE_3vec_SEEN */