#include "point.h"
#include <TMath.h>

using TMath::Sqrt;

Point::Point(){
	x = 0;
	y = 0;
	z = 0;
}
Point::Point(double x_, double y_, double z_){
	x = x_;
	y = y_;
	z = z_;
}
Point::Point(const Point& other){
	x = other.x;
	y = other.y;
	z = other.z;
}
Point& Point::operator=(const Point& other){
	x = other.x;
	y = other.y;
	z = other.z;
	return *this;
}
double Point::get_X() const{
	return x;
}
double Point::get_Y() const{
	return y;
}
double Point::get_Z() const{
	return z;
}
double Point::norm() const{
	return Sqrt((x*x)+(y*y)+(z*z));
}
double Point::normSquare() const{
	return (x*x)+(y*y)+(z*z);
}
Point Point::operator-() const{
	return Point(-x,-y,-z);
}
void Point::operator+=(const Point& other){
	x+=other.x;
	y+=other.y;
	z+=other.z;
}
void Point::operator-=(const Point& other){
	x-=other.x;
	y-=other.y;
	z-=other.z;
}
void Point::operator*=(const double& f){
	x*=f;
	y*=f;
	z*=f;
}
void Point::operator/=(const double& f){
	x/=f;
	y/=f;
	z/=f;
}

Point operator+(const Point& P1, const Point& P2){
	Point copie(P1);
	copie+=P2;
	return copie;
}
Point operator-(const Point& P1, const Point& P2){
	Point copie(P1);
	copie-=P2;
	return copie;
}
Point operator*(const Point& P, const double& f){
	Point copie(P);
	copie*=f;
	return copie;
}
Point operator/(const Point& P, const double& f){
	Point copie(P);
	copie/=f;
	return copie;
}
Point operator*(const double& f, const Point& P){
	Point copie(P);
	copie*=f;
	return copie;
}

bool operator<(const Point& P1, const Point& P2){
	if(P1.get_Z()<P2.get_Z()) return true;
	else if(P1.get_Z() == P2.get_Z()){
		if(P1.get_X()<P2.get_X()) return true;
		else if(P1.get_X()==P2.get_X()){
			if(P1.get_Y()<P2.get_Y()) return true;
			else return false;
		}
		else return false;
	}
	else return false;
}

double scalar_product(const Point& P1, const Point& P2){
	return ((P1.x)*(P2.x)) + ((P1.y)*(P2.y)) + ((P1.z)*(P2.z));
}
ostream& operator<<(ostream& os, const Point& P){
	os << "(" << P.x << ", " << P.y << ", " << P.z << ")";
	return os;
}