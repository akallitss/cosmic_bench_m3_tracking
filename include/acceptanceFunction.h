#ifndef acceptanceFunction_h
#define acceptanceFunction_h
#include <TH2D.h>
#include <set>
#include <vector>

using std::set;
using std::vector;

//allow to compute the accpetance and cos^2 factor for rectangular detectors in space coordinates
class acceptanceFunction{
	public:
		//build object with the position of the two extremal detectors hit by the tracks
		//bench_angle is the angle of the z axis with respect to the zenith
		acceptanceFunction(double x_min_Up_,double x_max_Up_,double y_min_Up_,double y_max_Up_,double x_min_Down_,double x_max_Down_,double y_min_Down_,double y_max_Down_,double z_Up_x_,double z_Down_x_,double z_Up_y_,double z_Down_y_,double bench_angle_);
		~acceptanceFunction();
		//make a 3D plot of the factor between hardcoded coordinates
		void plot_3D();
		//compute the factor at a plane crossing the z axis at "z" with an angle around the y axis "y_angle"
		//the 5 first arguments tweak the resulting histogram binning
		TH2D plot_XY(int nbin_x,double x1,double x2,int nbin_y,double y1,double y2, double z, double y_angle = 0);
		//TH2D plot_XY(int nbin_x, int nbin_y, double z, double y_angle = 0);
		//actual operator to compute the factor for tracks passing by the (x,y,z) point
		double operator()(double x,double y,double z);
		double x_min_Up;
		double x_max_Up;
		double y_min_Up;
		double y_max_Up;
		double x_min_Down;
		double x_max_Down;
		double y_min_Down;
		double y_max_Down;
		double z_Up_x;
		double z_Down_x;
		double z_Up_y;
		double z_Down_y;
		double bench_angle;
};
//allow to compute the accpetance and cos^2 factor for rectangular detectors in angular coordinates
class FreeSkyFunction{
	public:
		//build object with the dimensions of the detectors and their position on the z axis
		//bench_angle is the angle of the z axis with respect to the zenith
		FreeSkyFunction(double x_min_,double x_max_,double y_min_,double y_max_,vector<double> z_, double bench_angle_);
		~FreeSkyFunction();
		//unimplemented
		void plot_3D();
		//compute the factor for tracks between phi1 and phi1 and theta1 and theta2 and passing through at least "mult" different detectors
		TH2D plot_PhiTheta(int nbin_phi,double phi1,double phi2,int nbin_theta,double theta1,double theta2,unsigned int mult);
		//same as before with hardcoded limits in phi and theta
		TH2D plot_PhiTheta(int nbin_phi, int nbin_theta, unsigned int mult);
		//actual operator to compute the factor for (phi,theta) tracks passing through detectors distant by delta_z
		double operator()(double phi, double theta, double delta_z);
		double x_min;
		double x_max;
		double y_min;
		double y_max;
		double bench_angle;
		set<double> z;
};

#endif