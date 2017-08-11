#ifndef point_h
#define point_h
#include <ostream>

#include <TEllipse.h>

using std::ostream;

//3D point/vector class
class Point;

//vector cross product
Point operator*(const Point& P1, const Point& P2);
Point operator+(const Point& P1, const Point& P2);
Point operator-(const Point& P1, const Point& P2);
Point operator*(const Point& P, const double& f);
Point operator/(const Point& P, const double& f);
Point operator*(const double& f, const Point& P);
//lexicographic order
bool operator<(const Point& P1, const Point& P2);
//vector dot product
double scalar_product(const Point& P1, const Point& P2);
ostream& operator<<(ostream& os, const Point& P);


class Point{
	friend double scalar_product(const Point& P1, const Point& P2);
	friend Point operator*(const Point& P1, const Point& P2);
	friend ostream& operator<<(ostream& os, const Point& P);
	public:
		Point();
		//coordinate constructor
		Point(double x_, double y_, double z_);
		//copy constructor
		Point(const Point& other);
		//get X coord
		double get_X() const;
		//get Y coord
		double get_Y() const;
		//get Z coord
		double get_Z() const;
		//get vector norm (point distance to origin)
		double norm() const;
		//get square of norm()
		double normSquare() const;
		//oposite operator
		Point operator-() const;
		void operator+=(const Point& other);
		void operator-=(const Point& other);
		void operator*=(const double& f);
		void operator/=(const double& f);
		//copy assignment
		Point& operator=(const Point& other);
		//check if vector is null (point is origin)
		bool is_null() const;
	protected:
		double x;
		double y;
		double z;
};

//3D straight line class
class Line{
	public:
		Line();
		//constructor using 2 points by which the line pass
		Line(Point first,Point second);
		//Line(Point origin_,Point direction_);
		Line(const Line& other);
		Line& operator=(const Line& other);
		//check if line is parallel to the given one
		bool is_parallel(const Line& other) const;
		//check if line is defining a proper plane with the other
		bool is_coplanar(const Line& other) const;
		//get point of closest approach with the other line (middle of shortest segment between them)
		Point PoCA(const Line& other) const;
		//get line direction vector
		Point get_direction() const;
		//get a reference point by which the line pass
		Point get_origin() const;
	protected:
		Point origin;
		Point direction;
};

//2D plane in 3D space
class Plane{
	public:
		Plane();
		//constructor with a plane point and a normal vector
		Plane(Point norm, Point origin);
		//constructor with 2 lines defining a plane
		Plane(Line first, Line second);
		//constructor using the ax+by+cz+d=0 equation
		Plane(double a_, double b_, double c_, double d_);
		//copy constructor
		Plane(const Plane& other);
		//copy assignment
		Plane& operator=(const Plane& other);
		//get a factor
		double get_a() const;
		//get b factor
		double get_b() const;
		//get c factor
		double get_c() const;
		//get d factor
		double get_d() const;
		//get plane normal vector
		Point get_norm() const;
		//check if plane is parallel to the given one
		bool is_parallel(Plane other) const;
		//check if plane is parallel to the given line
		bool is_parallel(Line other) const;
		//compute intersection point between plane and given line
		Point intersection(Line other) const;
	protected:
		double a;
		double b;
		double c;
		double d;
};

//2D vector/point class
class Point_2D;

Point_2D operator+(const Point_2D& P1, const Point_2D& P2);
Point_2D operator-(const Point_2D& P1, const Point_2D& P2);
Point_2D operator*(const Point_2D& P, const double& f);
Point_2D operator/(const Point_2D& P, const double& f);
Point_2D operator*(const double& f, const Point_2D& P);
//lexicographic order
bool operator<(const Point_2D& P1, const Point_2D& P2);
double scalar_product(const Point_2D& P1, const Point_2D& P2);
//area of the parallelogram defined by the 2 vectors
double det(const Point_2D& P1, const Point_2D& P2);
ostream& operator<<(ostream& os, const Point_2D& P);

class Point_2D{
	friend double scalar_product(const Point_2D& P1, const Point_2D& P2);
	friend double det(const Point_2D& P1, const Point_2D& P2);
	friend ostream& operator<<(ostream& os, const Point_2D& P);
	public:
		Point_2D();
		//coordinate constructor
		Point_2D(double x_, double y_);
		//copy constructor
		Point_2D(const Point_2D& other);
		//get X coord
		double get_X() const;
		//get Y coord
		double get_Y() const;
		//get vector norm (point distance to origin)
		double norm() const;
		//get square of norm()
		double normSquare() const;
		//oposite operator
		Point_2D operator-() const;
		void operator+=(const Point_2D& other);
		void operator-=(const Point_2D& other);
		void operator*=(const double& f);
		void operator/=(const double& f);
		//copy assignment
		Point_2D& operator=(const Point_2D& other);
		//check if vector is null (point is origin)
		bool is_null() const;
		//check if the point is inside the given ellipse
		bool is_inside(TEllipse shape) const;
	protected:
		double x;
		double y;
};

//2D straight line class
class Line_2D{
	public:
		Line_2D();
		//constructor using 2 points by which the line pass
		Line_2D(Point_2D first,Point_2D second);
		//Line_2D(Point_2D origin_,Point_2D direction_);
		//copy constructor
		Line_2D(const Line_2D& other);
		//copy assignment
		Line_2D& operator=(const Line_2D& other);
		//check if line is parallel to the given one
		bool is_parallel(const Line_2D& other) const;
		//compute the intersetion point between this line and the given one
		Point_2D intersection(const Line_2D& other) const;
	//protected:
		Point_2D origin;
		Point_2D direction;
};

#endif