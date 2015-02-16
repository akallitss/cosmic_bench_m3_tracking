#ifndef point_h
#define point_h
#include <ostream>

using std::ostream;

class Point;

Point operator+(const Point& P1, const Point& P2);
Point operator-(const Point& P1, const Point& P2);
Point operator*(const Point& P, const double& f);
Point operator/(const Point& P, const double& f);
Point operator*(const double& f, const Point& P);
bool operator<(const Point& P1, const Point& P2);
double scalar_product(const Point& P1, const Point& P2);
ostream& operator<<(ostream& os, const Point& P);


class Point{
	friend double scalar_product(const Point& P1, const Point& P2);
	friend ostream& operator<<(ostream& os, const Point& P);
	public:
		Point();
		Point(double x_, double y_, double z_);
		Point(const Point& other);
		double get_X() const;
		double get_Y() const;
		double get_Z() const;
		double norm() const;
		double normSquare() const;
		Point operator-() const;
		void operator+=(const Point& other);
		void operator-=(const Point& other);
		void operator*=(const double& f);
		void operator/=(const double& f);
		Point& operator=(const Point& other);
	protected:
		double x;
		double y;
		double z;
};

#endif