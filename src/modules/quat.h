#ifndef FILE_quat_SEEN
#define FILE_quat_SEEN
/*---------------------------------------*/
#include<stdlib.h>
#include <iostream>
#include <math.h>
#include <array>
#include <iomanip>
#include "vec3.h"
/*---------------------------------------*/
class quat{
 private:
  double w {1};
	vec3 vv {0., 0., 0.};
 public:
	// Constructors
	quat() = default;
	quat(const double& r, const double& x, const double& y, const double& z){
		w = r;
		vv =vec3(x,y,z);
	}
	quat(const double& r, const std::array<double, 3>& uu){
		w = r;
		vv = vec3(uu);
	}
	quat(const double& r, const vec3& uu){
		w = r;
		vv = uu;
	}
#define qzero quat(0, 0, 0, 0);
	// copy constructor
	quat(const quat& Q){
		w=Q.w;
		vv=Q.vv;
	}
	double r() const{return w;}
	double x() const{return vv.x();}
	double y() const{return vv.y();}
	double z() const{return vv.z();}
	vec3 vec() const{return vv;}
	// functions
	inline quat operator*(const quat& Q){
		return quat{w*Q.w - dot(vv,Q.vv),
				w*Q.vv + Q.w*vv + cross(vv, Q.vv)};
	}
	inline quat operator+(const quat& Q){
		return quat{w+Q.w,
				             vv+Q.vv};
	}
	inline quat operator-(const quat& Q){
		return quat{w-Q.w,
				             vv-Q.vv};
	}
	void printq(){
		std::cout<< std::fixed << std::setprecision(2)
						 << std::setw(8) << w
						 << std::setw(8) << vv.x()
						 << std::setw(8) << vv.y()
						 << std::setw(8) << vv.z()<< std::endl;
	}

	/*----------------------------------*/
};
inline quat inv(const quat& Q){
	return quat(Q.r(), -Q.x(),  -Q.y(), -Q.z());
	}
inline quat operator*(const double& phi, const quat& Q){
	return quat{phi*Q.r(), phi*Q.vec()};
}
inline quat operator*(const vec3& uu, const quat& Q){
	return quat(0,uu)*Q;
}
inline quat operator*(quat Q, const vec3& uu){
	return Q*quat(0,uu);
}
quat rotation(const double& theta, const vec3& axis){
	normalize(axis);
	return quat( cos(theta/2.), axis*sin(theta/2.) );
}
/* COUNTERCLOCKWISE rotation */
vec3 rotate(vec3& uu, const double& theta, const vec3& axis){
	quat qt = rotation( theta, axis) ;
	quat qtemp = qt*uu*inv(qt);
	return qtemp.vec();
}
vec3 rotate(vec3& uu, const quat& q){
	quat qtemp = q*uu*inv(q);
	return qtemp.vec();
}
/* find out the quaternion that take uu to vv */
quat find_quat(const vec3& uu, const vec3& vv){
	double Cth = dot(uu,vv)/(norm(uu)*norm(vv));
	quat qq {1, 0, 0, 0};
	if (Cth != 1.){
		double Csqrthby2 = (1 + Cth)/2.;
		double Ssqrthby2 = 1 - Csqrthby2;
		vec3 axis = normalize(cross(uu,vv));
		qq= quat(sqrt(Csqrthby2), sqrt(Ssqrthby2)*axis);
	}
	return qq;
}
/*---------------------------------------------------*/
#endif /* !FILE_quat_SEEN */
