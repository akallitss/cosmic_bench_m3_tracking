#ifndef acceptanceFunction_h
#define acceptanceFunction_h

class acceptanceFunction{
	public:
		acceptanceFunction(double x_min_,double x_max_,double y_min_,double y_max_,double z_Up_,double z_Down_,double sigma_theta_);
		~acceptanceFunction();
		void plot_3D();
		void plot_XY(double z);
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