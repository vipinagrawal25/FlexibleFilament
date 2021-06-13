#ifndef FILE_3vec_SEEN
#define FILE_3vec_SEEN
/*---------------------------------------*/
#include <iostream>
#include <math.h>
#include "2Tens.h"
using namespace std;
/*-------------------------*/
class vec3{
public:
  double x,y,z;
  vec3();
  vec3(int,int,int);
  vec3(float,float,float);
  vec3(double,double,double);
    // definining operators 
  vec3 operator+(vec3);
  vec3 operator-(vec3);
  vec3 operator*(double);
  vec3 operator/(double);
  Tens2 operator*(vec3);
private:
};
vec3::vec3(){
  x = 0.0;
  y = 0.0;
  z = 0.0;
}
vec3::vec3(int a, int b, int c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
vec3::vec3(float a, float b, float c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
vec3::vec3(double a, double b, double c){
  x =  a;
  y =  b;
  z =  c;
}

vec3 vec3::operator+(vec3 param){
  vec3 temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  temp.z = z+param.z;
  return(temp);
}
vec3 vec3::operator-(vec3 param){
  vec3 temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  temp.z = z-param.z;
  return(temp);
}

vec3 vec3::operator*(double param){
  vec3 temp;
  temp.x=param*x;
  temp.y=param*y;
  temp.z=param*z;
  return(temp);
}

Tens2 vec3::operator*(vec3 param)
{   
    // This is a direct product of 2 vectors. Assuming vec3 a column vector with 3 entries
    // then this operator defines w = uv' (Where v' is the transpose of the matrix.)
    // It should be noted that this product is not commutative.
    Tens2 temp;
    temp.xx = x*(param.x);
    temp.xy = x*(param.y);
    temp.xz = x*(param.z);

    temp.yx = y*(param.x);
    temp.yy = y*(param.y);
    temp.yz = y*(param.z);

    temp.zx = z*(param.x);
    temp.zy = z*(param.y);
    temp.zz = z*(param.z);

    return(temp);
}

vec3 vec3::operator/(double param){
  vec3 temp;
  temp.x=x/param;
  temp.y=y/param;
  temp.z=z/param;
  return(temp);
}

double dot(vec3 a, vec3 b){
  double temp;
  temp = a.x*b.x+a.y*b.y+a.z*b.z;
  return(temp);
}

vec3 cross(vec3 a, vec3 b){
  vec3 temp;
  temp.x = a.y*b.z-a.z*b.y;
  temp.y = -(a.x*b.z-b.x*a.z);
  temp.z = a.x*b.y-a.y*b.x;
  return(temp);
}

double norm(vec3 a){
  return( sqrt( a.x*a.x+a.y*a.y+a.z*a.z) );}
/*---------------------------------------*/
double sqnorm(vec3 a){
  return( a.x*a.x+a.y*a.y+a.z*a.z) ;}
/*---------------------------------------*/
void PVec3(vec3 a){
  cout<<a.x<<"\t"<<a.y<<"\t"<<a.z<<"\n";
}
/*---------------------------------------*/
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
/*---------------------------------------*/
#endif /* !FILE_3vec_SEEN */

