#ifndef acceptanceFunction_h
#define acceptanceFunction_h
#include <TH2D.h>

class acceptanceFunction{
	public:
		acceptanceFunction(double x_min_,double x_max_,double y_min_,double y_max_,double z_Up_,double z_Down_,double sigma_theta_);
		~acceptanceFunction();
		void plot_3D();
		TH2D plot_XY(int nbin_x,double x1,double x2,int nbin_y,double y1,double y2, double z);
		double operator()(double x,double y,double z);
		double x_min;
		double x_max;
		double y_min;
		double y_max;
		double z_Up;
		double z_Down;
		double sigma_theta;
};

#endif