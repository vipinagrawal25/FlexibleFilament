#ifndef FILE_2Tens_SEEN
#define FILE_2Tens_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
using namespace std;
/*-------------------------*/
class Tens2{
public:
  double xx,xy,xz,yx,yy,yz,zx,zy,zz;
  Tens2();
  Tens2(int,int,int);
  Tens2(float,float,float);
  Tens2(double,double,double);
    // definining operators 
  Tens2 operator+(Tens2);
  Tens2 operator-(Tens2);
  Tens2 operator*(double);
  Tens2 operator/(double);
};

Tens2::Tens2(){
  xx = 0.0;
  xy = 0.0;
  xz = 0.0;
  
  yx = 0.0;
  yy = 0.0;
  yz = 0.0;
  
  zx = 0.0;
  zy = 0.0;
  zz = 0.0;
}

Tens2::Tens2(int a, int b, int c, int d, int e, int f, int g, int h, int i){
  xx = (double) a;
  xy = (double) b;
  xz = (double) c;
  
  yx = (double) d;
  yy = (double) e;
  yz = (double) f;
  
  zx = (double) g;
  zy = (double) h;
  zz = (double) i;
}

Tens2::Tens2(float a, float b, float c, float d, float e, float f, float g, float h, float i){
  xx = (double) a;
  xy = (double) b;
  xz = (double) c;
  
  yx = (double) d;
  yy = (double) e;
  yz = (double) f;
  
  zx = (double) g;
  zy = (double) h;
  zz = (double) i;
}

Tens2::Tens2(double a, doubleb, double c, double d, double e, double f, double g, double h, double i){
  xx =  a;
  xy =  b;
  xz =  c;

  yx =  d;
  yy =  e;
  yz =  f;
  
  zx =  g;
  zy =  h;
  zz =  i;
}

Tens2 Tens2::operator+(Tens2 param){
  Tens2 temp;
  temp.xx = xx+param.xx;
  temp.xy = xy+param.xy;
  temp.xz = xz+param.xz;

  temp.yx = yx+param.yx;
  temp.yy = yy+param.yy;
  temp.yz = yz+param.yz;
  
  temp.zx = zx+param.zx;
  temp.zy = zy+param.zy;
  temp.zz = zz+param.zz;

  return(temp);
}

Tens2 Tens2::operator-(Tens2 param){
  Tens2 temp;
  temp.xx = xx-param.xx;
  temp.xy = xy-param.xy;
  temp.xz = xz-param.xz;

  temp.yx = yx-param.yx;
  temp.yy = yy-param.yy;
  temp.yz = yz-param.yz;
  
  temp.zx = zx-param.zx;
  temp.zy = zy-param.zy;
  temp.zz = zz-param.zz;

  return(temp);
}

Tens2 Tens2::operator*(double param){
  Tens2 temp;
  temp.xx=param*xx;
  temp.xy=param*xy;
  temp.xz=param*xz;

  temp.yx=param*yx;
  temp.yy=param*yy;
  temp.yz=param*yz;

  temp.zx=param*zx;
  temp.zy=param*zy;
  temp.zz=param*zz;

  return(temp);
}

Tens2 operator*(const double param, const Tens2 a)
{
    // This function is just here because multiplication between double and a Tensor needs to be commutative.
    return a*param;
}

Tens2 Tens2::operator/(double param){
  Tens2 temp;
  temp.xx=param/xx;
  temp.xy=param/xy;
  temp.xz=param/xz;

  temp.yx=param/yx;
  temp.yy=param/yy;
  temp.yz=param/yz;

  temp.zx=param/zx;
  temp.zy=param/zy;
  temp.zz=param/zz;

  return(temp);
}

vec3 dot(Tens2 a, vec3 b)
{
    // Ci = Aij*Bj; (Einstein Convention)
    // It should be a simple matrix multiplication, where A is 3X3 matrix and B is 3X1.

    vec3 temp;

    temp.x = (a.xx)*(b.x)+(a.xy)*(b.y)+(a.xz)*(b.z);
    temp.y = (a.yx)*(b.x)+(a.yy)*(b.y)+(a.yz)*(b.z);
    temp.z = (a.zx)*(b.x)+(a.zy)*(b.y)+(a.zz)*(b.z);

    return (temp);
}

// This defines the kronecker delta function or a unit tensor matrix. dab = \delta_{\alpha \beta}
Tens2 dab(1,0,0,0,1,0,0,0,1);

/*---------------------------------------*/
#endif /* !FILE_3vec_SEEN */

