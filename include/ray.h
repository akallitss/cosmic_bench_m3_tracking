#ifndef ray_h
#define ray_h
#include <vector>
#include <ostream>

#include "point.h"

using std::vector;
using std::ostream;

class Ray;
class Ray_2D;
class Cluster;
class Detector;

ostream& operator<<(ostream& os, const Ray_2D& ray);
ostream& operator<<(ostream& os, const Ray& ray);

class Ray_2D{
	friend class Ray;
	friend ostream& operator<<(ostream& os, const Ray_2D& ray);
	public:
		Ray_2D();
		Ray_2D(char coord_);
		Ray_2D(const Ray_2D& other);
		Ray_2D& operator=(const Ray_2D& other);
		Ray_2D(const Ray& other, char coord_);
		~Ray_2D();
		void add_cluster(Cluster * clus);
		void process();
		double get_chiSquare() const;
		double get_slope() const;
		double get_Z_intercept() const;
		double eval(double z) const;
		double get_residu(Detector * det) const;
		double get_residu_ref(Cluster * clus) const;
		double get_t_mean() const;
		double get_t_sigma() const;
		unsigned int get_clus_n() const;
		void clear();
	protected:
		vector<Cluster*> clusters;
		double chiSquare;
		double slope;
		double Z_intercept;
		char coord;
};

class Ray{
	friend class Ray_2D;
	friend class RayPair;
	friend ostream& operator<<(ostream& os, const Ray& ray);
	public:
		Ray();
		Ray(const Ray& other);
		Ray& operator=(const Ray& other);
		Ray(const Ray_2D& ray1, const Ray_2D& ray2);
		~Ray();
		void add_cluster(Cluster * clus);
		void process();
		void angle_correction();
		double get_chiSquare_X() const;
		double get_chiSquare_Y() const;
		double get_slope_X() const;
		double get_slope_Y() const;
		double get_Z_intercept_X() const;
		double get_Z_intercept_Y() const;
		double eval_X(double z) const;
		double eval_Y(double z) const;
		double eval_X(Detector * det) const;
		double eval_Y(Detector * det) const;
		Point eval_plane(Plane proj) const;
		double get_residu(Detector * det) const;
		double get_residu_ref(Cluster * clus) const;
		double get_t_mean() const;
		double get_t_sigma() const;
		unsigned int get_clus_n() const;
		unsigned int get_clus_x_n() const;
		unsigned int get_clus_y_n() const;
		vector<Cluster*> get_clus() const;
		void clear();
	protected:
		vector<Cluster*> clusters;
		double chiSquare_X;
		double chiSquare_Y;
		double slope_X;
		double slope_Y;
		double Z_intercept_X;
		double Z_intercept_Y;
};

class RayPair{
	public:
		RayPair();
		RayPair(const RayPair& other);
		RayPair& operator=(const RayPair& other);
		RayPair(const Ray& ray1, const Ray& ray2);
		~RayPair();
		void process();
		double get_delta_x() const;
		double get_delta_y() const;
		double get_delta_theta_x() const;
		double get_delta_theta_y() const;
		double get_theta_x_up() const;
		double get_theta_y_up() const;
		double get_theta_x_down() const;
		double get_theta_y_down() const;
		double get_x_up(double z) const;
		double get_y_up(double z) const;
		double get_x_down(double z) const;
		double get_y_down(double z) const;
		double get_doca();
		Point get_PoCA();
		Ray downRay;
		Ray upRay;
	protected:

		//after processing
		double delta_x;
		double delta_y;
		double delta_theta_x;
		double delta_theta_y;
};

#endif