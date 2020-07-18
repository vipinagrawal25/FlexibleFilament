#ifndef FILE_3vec_SEEN
#define FILE_3vec_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "2b2Tens.h"
using namespace std;
/*-------------------------*/
class vec2{
public:
  double x,y,z;
  vec2();
  vec2(int,int);
  vec2(float,float);
  vec2(double,double);
    // definining operators 
  vec2 operator+(vec2);
  vec2 operator-(vec2);
  vec2 operator*(double);
  vec2 operator/(double);
  Tens2b2 operator*(vec2);

private:
};
vec2::vec2(){
  x = 0.0;
  y = 0.0;
}
vec2::vec2(int a, int b){
  x = (double) a;
  y = (double) b;
}
vec2::vec2(float a, float b){
  x = (double) a;
  y = (double) b;
}
vec2::vec2(double a, double b){
  x =  a;
  y =  b;
}

vec2 vec2::operator+(vec2 param){
  vec2 temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  return(temp);
}
vec2 vec2::operator-(vec2 param){
  vec2 temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  return(temp);
}

vec2 vec2::operator*(double param){
  vec2 temp;
  temp.x=param*x;
  temp.y=param*y;
  return(temp);
}

Tens2b2 vec2::operator*(vec2 param)
{   
    // This is a direct product of 2 vectors. Assuming vec2 a column vector with 3 entries
    // then this operator defines w = uv' (Where v' is the transpose of the matrix.)
    // It should be noted that this product is not commutative.

    Tens2b2 temp;
    temp.xx = x*(param.x);
    temp.xy = x*(param.y);

    temp.yx = y*(param.x);
    temp.yy = y*(param.y);

    return(temp);
}

vec2 vec2::operator/(double param){
  vec2 temp;
  temp.x=x/param;
  temp.y=y/param;
  return(temp);
}

double dot(vec2 a, vec2 b){
  double temp;
  temp = a.x*b.x+a.y*b.y;
  return(temp);
}
double cross(vec2 a, vec2 b){
  double temp;
  temp = a.x*b.y-a.y*b.x;
  return(temp);
}
double norm(vec2 a){
  return( sqrt( a.x*a.x+a.y*a.y) );}
/*---------------------------------------*/
double sqnorm(vec2 a){
  return( a.x*a.x+a.y*a.y) ;}
/*---------------------------------------*/
void PVec2(vec2 a){
  cout<<a.x<<"\t"<<a.y<<"\n";
}

vec2 dot(Tens2b2 a, vec2 b){
    vec2 temp;
    temp.x = (a.xx)*(b.x)+(a.xy)*(b.y);
    temp.y = (a.yx)*(b.x)+(a.yy)*(b.y);

    return (temp);
}
/*---------------------------------------*/
#endif /* !FILE_3vec_SEEN */
