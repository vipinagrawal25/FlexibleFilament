#ifndef FILE_vec3_SEEN
#define FILE_vec3_SEEN
/*---------------------------------------*/
#include<stdlib.h>
#include <iostream>
#include <math.h>
#include <array>
#include <iomanip>
/*---------------------------------------*/
class vec3{
 private:
	std::array<double,3> a {0., 0., 0.};
 public:
	// Constructors
	vec3() = default;
	vec3(const double& x, const double& y, const double& z){
		a[0]=x;
		a[1]=y;
		a[2]=z;
	}
	vec3(const std::array<double, 3>& A){
		for( size_t i {}; i<3; ++i){
		a[i]=A[i];
		}
	}
	// copy constructor
	vec3(const vec3& vv){
		for( size_t i {}; i<3; ++i){
		a[i]=vv.a[i];
		}
	}
	// Accessors
	double x() const{return a[0];}
	double y() const{return a[1];}
	double z() const{return a[2];}
	std::array<double,3> geta() const{return a;}
	/*----------ALGEBRA------------------------*/
	vec3 operator+(const vec3& A){
		return vec3 {a[0]+A.a[0],
				             a[1]+A.a[1],
				a[2]+A.a[2]};
	}
	vec3 operator-(const vec3& A){
		return vec3 {a[0]-A.a[0],
				             a[1]-A.a[1],
				a[2]-A.a[2]};
	}

	/*----------------------------------*/
};
vec3 operator/(const vec3& A, const double& phi){
	return vec3 {A.x()/phi, A.y()/phi, A.z()/phi};
}
vec3 operator*(const double& phi, const vec3& A){
	std::array<double,3> C;
	for(size_t i{}; i<3; ++i) {
		C[i] = phi*A.geta()[i]; 
	}
	return vec3(C);
}
vec3 operator*(const vec3& A, const double& phi){
	return phi*A;
}
double dot(const vec3& A, const vec3& B){
	double C = 0;
	for(size_t i{}; i<3; ++i) {
		C +=  A.geta()[i]*B.geta()[i];
	}
	return C;
}
inline double norm(const vec3& AA){
	return sqrt(dot(AA, AA)) ;
}
void normalize(vec3& AA){
	AA = AA/norm(AA);
}
vec3 normalize(const vec3& AA){
  return AA/norm(AA);
}
	/*-------------------------*/
vec3 cross(const vec3& A, const vec3& B){
	return vec3{ A.y()*B.z() - A.z()*B.y() , 
			-A.x()*B.z() + A.z()*B.x() , 
			A.x()*B.y() - A.y()*B.x()
			}; 
}
void printv(const vec3& A){
	std::cout<< std::fixed << std::setprecision(2)
					 << std::setw(8) << A.x()
					 << std::setw(8) << A.y()
					 << std::setw(8) << A.z()<< std::endl;
}
/*---------------------------------------------------*/
#endif /* !FILE_vec3_SEEN */
