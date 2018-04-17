#ifndef FILE_3vec_SEEN
#define FILE_3vec_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
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
#endif /* !FILE_3vec_SEEN */

