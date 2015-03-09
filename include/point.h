#ifndef point_h
#define point_h
#include <ostream>

using std::ostream;

class Point;

Point operator*(const Point& P1, const Point& P2);
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
	friend Point operator*(const Point& P1, const Point& P2);
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
		bool is_null() const;
	protected:
		double x;
		double y;
		double z;
};

class Line{
	public:
		Line();
		Line(Point first,Point second);
		//Line(Point origin_,Point direction_);
		Line(const Line& other);
		bool is_parallel(const Line& other) const;
		bool is_coplanar(const Line& other) const;
		Point PoCA(const Line& other) const;
	protected:
		Point origin;
		Point direction;
};


class Point_2D;

Point_2D operator+(const Point_2D& P1, const Point_2D& P2);
Point_2D operator-(const Point_2D& P1, const Point_2D& P2);
Point_2D operator*(const Point_2D& P, const double& f);
Point_2D operator/(const Point_2D& P, const double& f);
Point_2D operator*(const double& f, const Point_2D& P);
bool operator<(const Point_2D& P1, const Point_2D& P2);
double scalar_product(const Point_2D& P1, const Point_2D& P2);
double det(const Point_2D& P1, const Point_2D& P2);
ostream& operator<<(ostream& os, const Point_2D& P);

class Point_2D{
	friend double scalar_product(const Point_2D& P1, const Point_2D& P2);
	friend double det(const Point_2D& P1, const Point_2D& P2);
	friend ostream& operator<<(ostream& os, const Point_2D& P);
	public:
		Point_2D();
		Point_2D(double x_, double y_);
		Point_2D(const Point_2D& other);
		double get_X() const;
		double get_Y() const;
		double norm() const;
		double normSquare() const;
		Point_2D operator-() const;
		void operator+=(const Point_2D& other);
		void operator-=(const Point_2D& other);
		void operator*=(const double& f);
		void operator/=(const double& f);
		Point_2D& operator=(const Point_2D& other);
		bool is_null() const;
	protected:
		double x;
		double y;
};

class Line_2D{
	public:
		Line_2D();
		Line_2D(Point_2D first,Point_2D second);
		//Line_2D(Point_2D origin_,Point_2D direction_);
		Line_2D(const Line_2D& other);
		bool is_parallel(const Line_2D& other) const;
		Point_2D intersection(const Line_2D& other) const;
	//protected:
		Point_2D origin;
		Point_2D direction;
};

#endif