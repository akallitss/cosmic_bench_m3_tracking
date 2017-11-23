#ifndef ray_h
#define ray_h
#include <vector>
#include <ostream>
#include <utility>

#include "point.h"

using std::vector;
using std::ostream;
using std::pair;

class Ray;
class Ray_2D;
class Cluster;
class Detector;
class CosmicBench;

ostream& operator<<(ostream& os, const Ray_2D& ray);
ostream& operator<<(ostream& os, const Ray& ray);

//2D track class
class Ray_2D{
	friend class Ray;
	friend ostream& operator<<(ostream& os, const Ray_2D& ray);
	public:
		Ray_2D();
		//constructor setting the projection plane, X for xz-plane and Y for yz-plane
		Ray_2D(char coord_);
		//copy constructor
		Ray_2D(const Ray_2D& other);
		//copy assignment
		Ray_2D& operator=(const Ray_2D& other);
		//extract the needed projection from the given 3D track
		Ray_2D(const Ray& other, char coord_);
		~Ray_2D();
		//add a cluster to the track
		void add_cluster(const Cluster * const clus);
		//process the track parameter by fitting it to the clusters
		void process();
		//get the track Chi2
		double get_chiSquare() const;
		//get the dz/dcoord
		double get_slope() const;
		//get the coord at which the track cross the Oxy-plane
		double get_Z_intercept() const;
		//get the coord at the given altitude
		double eval(double z) const;
		//get the residual in the given detector (given detector is included into tracking so use it carefully)
		double get_residu(const Detector * const det) const;
		//get the distance between the given cluster and the track
		double get_residu_ref(const Cluster * const clus) const;
		//get mean time of track clusters
		double get_t_mean() const;
		//get time standard deviation of track clusters
		double get_t_sigma() const;
		//get number of cluster used to build the track
		unsigned int get_clus_n() const;
		//get the clusters used to build the track
		vector<Cluster*> get_clus() const;
		//get the index of the extremal detectors (wrt. z axis) which are used in the track (first is lower, second is higher)
		pair<int,int> get_extremal_det(const CosmicBench * const bench) const;
		//delete the clusters
		void clear();
		// check if ray use cluster in given layer
		bool has_layer(int layer) const;
	protected:
		vector<Cluster*> clusters;
		double chiSquare;
		double slope;
		double Z_intercept;
		char coord;
};

//3D track class
class Ray{
	friend class Ray_2D;
	friend class RayPair;
	friend ostream& operator<<(ostream& os, const Ray& ray);
	public:
		Ray();
		//copy constructor
		Ray(const Ray& other);
		//copy assignment
		Ray& operator=(const Ray& other);
		//merge constructor using an X and Y 2D track
		Ray(const Ray_2D& ray1, const Ray_2D& ray2);
		~Ray();
		//add a cluster to the track
		void add_cluster(const Cluster * const clus);
		//process the track parameter by fitting it to the clusters
		void process();
		//make the angular alignment correction of detectors
		void angle_correction();
		//get the xz projection Chi2
		double get_chiSquare_X() const;
		//get the yz projection Chi2
		double get_chiSquare_Y() const;
		//get dz/dx of track
		double get_slope_X() const;
		//get dz/dy of track
		double get_slope_Y() const;
		//get the x coord at which the track cross the Oxy-plane
		double get_Z_intercept_X() const;
		//get the y coord at which the track cross the Oxy-plane
		double get_Z_intercept_Y() const;
		//get the x coord at the given altitude
		double eval_X(double z) const;
		//get the y coord at the given altitude
		double eval_Y(double z) const;
		//get the x coord at the given detector altitude
		double eval_X(const Detector * const det) const;
		//get the y coord at the given detector altitude
		double eval_Y(const Detector * const det) const;
		//get the intersection between this track and the given plane
		Point eval_plane(Plane proj) const;
		//get the residual in the given detector (given detector is included into tracking so use it carefully)
		double get_residu(const Detector * const det) const;
		//get the distance between the given cluster and the track
		double get_residu_ref(const Cluster * const clus) const;
		//get mean time of track clusters
		double get_t_mean() const;
		//get time standard deviation of track clusters
		double get_t_sigma() const;
		//get the total number of cluster used to build the track
		unsigned int get_clus_n() const;
		//get number of cluster measuring the x coord used to build the track
		unsigned int get_clus_x_n() const;
		//get number of cluster measuring the y coord used to build the track
		unsigned int get_clus_y_n() const;
		//get the clusters used to build the track
		vector<Cluster*> get_clus() const;
		//get the index of the extremal detectors (wrt. z axis) which are used in the track (first.first is lower x , first.second is higher x, second.first is lower y , second.second is higher y)
		pair<pair<int,int>,pair<int,int> > get_extremal_det(const CosmicBench * const bench) const;
		//delete the clusters
		void clear();
		// check if ray use cluster in given layer
		bool has_layer(int layer) const;
	protected:
		vector<Cluster*> clusters;
		double chiSquare_X;
		double chiSquare_Y;
		double slope_X;
		double slope_Y;
		double Z_intercept_X;
		double Z_intercept_Y;
};

//Scattered track class
class RayPair{
	public:
		RayPair();
		//copy constructor
		RayPair(const RayPair& other);
		//copy assignment
		RayPair& operator=(const RayPair& other);
		//constructor using before and after scattering tracks
		RayPair(const Ray& ray1, const Ray& ray2);
		~RayPair();
		//compute the scattering infos
		void process();
		//compute delta x as defined by Schultz
		double get_delta_x() const;
		//compute delta y as defined by Schultz
		double get_delta_y() const;
		//compute delta theta x using the x slope difference
		double get_delta_theta_x() const;
		//compute delta theta y using the y slope difference
		double get_delta_theta_y() const;
		//get the x inclination of the incoming track
		double get_theta_x_up() const;
		//get the y inclination of the incoming track
		double get_theta_y_up() const;
		//get the x inclination of the outgoing track
		double get_theta_x_down() const;
		//get the y inclination of the outgoing track
		double get_theta_y_down() const;
		//get the x coord at the given altitude of the incoming track
		double get_x_up(double z) const;
		//get the y coord at the given altitude of the incoming track
		double get_y_up(double z) const;
		//get the x coord at the given altitude of the outgoing track
		double get_x_down(double z) const;
		//get the y coord at the given altitude of the outgoing track
		double get_y_down(double z) const;
		//compute the minimal distance between the two half-tracks
		double get_doca();
		//compute the point of closest approach as the middle of the smaller segment between the two half-trakcs
		Point get_PoCA();
		//incoming track
		Ray downRay;
		//outgoing track
		Ray upRay;
	protected:

		//after processing
		double delta_x;
		double delta_y;
		double delta_theta_x;
		double delta_theta_y;
};

#endif