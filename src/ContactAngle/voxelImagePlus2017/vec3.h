#ifndef VEC3_H
#define VEC3_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>

#ifndef verySmall
	#define verySmall  1.0e-30
#endif

class vec3
{
public:

	double	x;
	double	y;
	double	z;


	vec3() {}
//	vec3(double v[3])  { x = v[0];  y = v[1];  z = v[2]; }
	vec3(double r, double s, double t)  { x = r;  y = s;  z = t; }

	vec3& Set(double r, double s, double t)  { x = r;  y = s;  z = t;  return (*this); }

	double&       operator [](long k) {	return ((&x)[k]); }

	const double& operator [](long k) const  { return ((&x)[k]); }

	vec3& operator +=(const vec3& v)  { x += v.x;  y += v.y;  z += v.z;  return (*this); }
	//vec3& operator +=(double t)  { x += t;  y += t;  z += t;  return (*this); }
	vec3& operator -=(const vec3& v)  { x -= v.x;  y -= v.y;  z -= v.z;  return (*this); }
	//vec3& operator -=(double t)  { x -= t;  y -= t;  z -= t;  return (*this); }
	vec3& operator *=(double t)  { x *= t;  y *= t;  z *= t;  return (*this); }
	vec3& operator /=(double t)  {  double f = 1.0 / t;  x *= f;  y *= f;  z *= f;  return (*this); }
	vec3& operator ^=(const vec3& v)  { double r, s;  r=y*v.z-z*v.y;  s=z*v.x-x*v.z;  z=x*v.y-y*v.x;  x=r; y=s; 	return (*this); }
	vec3& operator *=(const vec3& v)  { x *= v.x;  y *= v.y;  z *= v.z;  return (*this); }
	vec3  operator -(void) const  { return (vec3(-x, -y, -z)); }
	vec3  operator +(const vec3& v) const  { return (vec3(x+v.x, y+v.y, z+v.z)); }
	//vec3  operator +(double t) const  { return (vec3(x+t, y+t, z+t)); }
	vec3  operator -(const vec3& v) const  { return (vec3(x-v.x, y-v.y, z-v.z)); }
	//vec3  operator -(double t) const  { return (vec3(x-t, y-t, z-t)); }
	vec3  operator *(double t) const  { return (vec3(x*t, y*t, z*t)); }
	vec3  operator /(double t) const  { double f = 1.0 / t;  return (vec3(x*f, y*f, z*f)); }
	double operator &(const vec3& v) const  { return (x*v.x+y*v.y+z*v.z); }
	vec3  operator ^(const vec3& v) const  { return (vec3(y*v.z-z*v.y,  z*v.x-x*v.z,  x*v.y-y*v.x)); }
	vec3  operator *(const vec3& v) const  { return (vec3(x*v.x, y*v.y, z*v.z)); }
	//~ bool operator ==(const vec3& v) const  { return ((x == v.x) && (y == v.y) && (z == v.z)); }
	//~ bool operator !=(const vec3& v) const  { return ((x != v.x) || (y != v.y) || (z != v.z)); }
	bool operator ==(const vec3& v) const  { return ((x-v.x)*(x-v.x) < verySmall) && ((y-v.y)*(y-v.y) < verySmall) && ((z-v.z)*(z-v.z) < verySmall); }
	bool operator !=(const vec3& v) const  { return ((x-v.x)*(x-v.x) >= verySmall) || ((y-v.y)*(y-v.y) >= verySmall) || ((z-v.z)*(z-v.z) >= verySmall); }

	vec3& Normalize(void)  { return (*this /= std::sqrt(x*x+y*y+z*z)); }
	vec3& RotateAboutX(double angle);
	vec3& RotateAboutY(double angle);
	vec3& RotateAboutZ(double angle);
	vec3& RotateAboutAxis(double angle, const vec3& axis);
};



inline vec3 operator *(double t, const vec3& v) { return vec3(t*v.x, t*v.y, t*v.z); }

inline double mag(const vec3& v) { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }

inline double magSqr(const vec3& v) { return (v.x*v.x+v.y*v.y+v.z*v.z); }
inline vec3 norm(const vec3& v)  { return  v/(mag(v)+1.0e-300); }

inline vec3& vec3::RotateAboutAxis(double angle, const vec3& axis)
{
	double s = sinf(angle);
	double c = cosf(angle);
	double k = 1.0 - c;
	double nx = x * (c + k * axis.x * axis.x) + y * (k * axis.x * axis.y - s * axis.z)	+ z * (k * axis.x * axis.z + s * axis.y);
	double ny = x * (k * axis.x * axis.y + s * axis.z) + y * (c + k * axis.y * axis.y)	+ z * (k * axis.y * axis.z - s * axis.x);
	double nz = x * (k * axis.x * axis.z - s * axis.y) + y * (k * axis.y * axis.z + s * axis.x)	+ z * (c + k * axis.z * axis.z);
	x = nx;  y = ny;  z = nz;
	return (*this);
}



#define int3x3  std::array<std::array<int,3>,3>
#define int3  std::array<int,3>

inline vec3  operator*( int3 n, vec3 x) { return vec3(n[0]*x[0], n[1]*x[1], n[2]*x[2]);}
inline int3 operator-( int3 n, int3 x) { return int3{{n[0]-x[0], n[1]-x[1], n[2]-x[2]}};}
inline int3& operator+=( int3& n, int3 x) { n[0]+=x[0]; n[1]+=x[1]; n[2]+=x[2]; return n;}






inline std::ostream& operator<< (std::ostream& out, const vec3& node)
{
	out.flags(std::ios::showpoint);
	out.flags(std::ios::scientific);
	out << std::setprecision(8)  << std::setw(19) << node.x  << std::setw(19) << node.y  << std::setw(19) << node.z;
	return out;
}

inline std::istream& operator>> (std::istream& in, vec3& node)
{
	in >> node.x >> node.y >> node.z;
	return in;
}



#define dblpair  std::pair<double,double>


template<typename T1, typename T2>
inline std::istream& operator>> (std::istream& in, std::pair<T1,T2>& interval)
{
	in >> interval.first >> interval.second;
	return in;
}

template<typename T1, typename T2>
inline std::ostream& operator<< (std::ostream& out, std::pair<T1,T2>& interval)
{
	out  << std::setw(15)<< interval.first <<std::setw(15)<< interval.second;
	return out;
}



#ifndef TOSTR
#define TOSTR
template<typename T> std::string toStr(const T& n){  std::ostringstream stm;  stm<<n;  return stm.str();  }
#endif



#endif


