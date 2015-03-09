#include "point.h"
#include <TMath.h>
#include <limits>

using TMath::Sqrt;
using TMath::Abs;
using std::numeric_limits;

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
bool Point::is_null() const{
	bool epsilon = numeric_limits<double>::epsilon();
	return (Abs(x) <= epsilon && Abs(y) <= epsilon && Abs(z) <= epsilon);
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

Point operator*(const Point& P1, const Point& P2){
	return Point(P1.y*P2.z - P1.z*P2.y,P1.z*P2.x - P1.x*P2.z,P1.x*P2.y - P1.y*P2.x);
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

Line::Line(){
	origin = Point();
	direction = Point();
}
Line::Line(Point first,Point second){
	origin = first;
	direction = second - first;
	direction /= direction.norm();
}
/*
Line::Line(Point origin_,Point direction_){
	origin = origin_;
	direction = direction_;
	direction /= direction.norm();
}
*/
Line::Line(const Line& other){
	origin = other.origin;
	direction = other.direction;
}
bool Line::is_parallel(const Line& other) const{
	return (direction*other.direction).is_null();
}
bool Line::is_coplanar(const Line& other) const{
	return (scalar_product(other.direction,direction*(origin - other.origin)) <= numeric_limits<double>::epsilon());
}
Point Line::PoCA(const Line& other) const{
	if(is_parallel(other)){
		return ((origin + other.origin)*0.5);
	}
	else{
		Point common_term = (other.origin-origin)/(direction.normSquare()*other.direction.normSquare()-scalar_product(direction,other.direction)*scalar_product(direction,other.direction));
		double param_closest_point = scalar_product(common_term,(direction*other.direction.normSquare()-other.direction*scalar_product(direction,other.direction)));
		double param_closest_point_other = -1*scalar_product(common_term,(other.direction*direction.normSquare()-direction*scalar_product(other.direction,direction)));
		Point closest_point = param_closest_point*direction + origin;
		Point closest_point_other = param_closest_point_other*other.direction + other.origin;
		return ((closest_point_other+closest_point)*0.5);
	}
}

Point_2D operator+(const Point_2D& P1, const Point_2D& P2){
	Point_2D copie(P1);
	copie += P2;
	return copie;
}
Point_2D operator-(const Point_2D& P1, const Point_2D& P2){
	Point_2D copie(P1);
	copie -= P2;
	return copie;
}
Point_2D operator*(const Point_2D& P, const double& f){
	Point_2D copie(P);
	copie *= f;
	return copie;
}
Point_2D operator/(const Point_2D& P, const double& f){
	Point_2D copie(P);
	copie /= f;
	return copie;
}
Point_2D operator*(const double& f, const Point_2D& P){
	Point_2D copie(P);
	copie *= f;
	return copie;
}
bool operator<(const Point_2D& P1, const Point_2D& P2){
	if(P1.get_X()<P2.get_X()) return true;
	else if(P1.get_X()==P2.get_X()){
		if(P1.get_Y()<P2.get_Y()) return true;
		else return false;
	}
	else return false;
}
double scalar_product(const Point_2D& P1, const Point_2D& P2){
	return ((P1.x*P2.x) + (P1.y*P2.y));
}
double det(const Point_2D& P1, const Point_2D& P2){
	return ((P1.x*P2.y) - (P2.x*P1.y));
}
ostream& operator<<(ostream& os, const Point_2D& P){
	os << "(" << P.x << ", " << P.y << ")";
	return os;
}

Point_2D::Point_2D(){
	x = 0;
	y = 0;
}
Point_2D::Point_2D(double x_, double y_){
	x = x_;
	y = y_;
}
Point_2D::Point_2D(const Point_2D& other){
	x = other.x;
	y = other.y;
}
double Point_2D::get_X() const{
	return x;
}
double Point_2D::get_Y() const{
	return y;
}
double Point_2D::norm() const{
	return Sqrt((x*x) + (y*y));
}
double Point_2D::normSquare() const{
	return ((x*x) + (y*y));
}
Point_2D Point_2D::operator-() const{
	return Point_2D(-x,-y);
}
void Point_2D::operator+=(const Point_2D& other){
	x += other.x;
	y += other.y;
}
void Point_2D::operator-=(const Point_2D& other){
	x -= other.x;
	y -= other.y;
}
void Point_2D::operator*=(const double& f){
	x *= f;
	y *= f;
}
void Point_2D::operator/=(const double& f){
	x /= f;
	y /= f;
}
Point_2D& Point_2D::operator=(const Point_2D& other){
	x = other.x;
	y = other.y;
	return *this;
}
bool Point_2D::is_null() const{
	bool epsilon = numeric_limits<double>::epsilon();
	return (Abs(x) <= epsilon && Abs(y) <= epsilon);
}

Line_2D::Line_2D(){
	origin = Point_2D();
	direction = Point_2D();
}
Line_2D::Line_2D(Point_2D first,Point_2D second){
	origin = first;
	direction = (second - first);
	direction /= direction.norm();
}
/*
Line_2D::Line_2D(Point_2D origin_,Point_2D direction_){
	origin = origin_;
	direction = direction_/(direction_.norm());
}
*/
Line_2D::Line_2D(const Line_2D& other){
	origin = other.origin;
	direction = other.direction;
}
bool Line_2D::is_parallel(const Line_2D& other) const{
	return (Abs(det(other.direction,direction)) <= numeric_limits<double>::epsilon());
}
Point_2D Line_2D::intersection(const Line_2D& other) const{
	if(is_parallel(other)) return ((origin + other.origin)*0.5);
	else{
		double angle = det(direction,other.direction);
		double common_term = det(origin,direction);
		double common_term_other = det(other.origin,other.direction);
		return ((direction*common_term_other - other.direction*common_term)/angle);
	}
}