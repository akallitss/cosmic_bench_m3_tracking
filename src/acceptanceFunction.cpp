#define acceptanceFunction_cpp

#include "acceptanceFunction.h"
#include "point.h"
#include "tomography.h"

#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRint.h>
#include <TF2.h>
#include <Math/WrappedMultiTF1.h>
#include <Math/AdaptiveIntegratorMultiDim.h>
#include <iostream>

using std::cout;
using std::endl;
using std::flush;

using TMath::ATan;
using TMath::Erf;
using TMath::Sqrt;
using TMath::Max;
using TMath::Min;
using TMath::Pi;
using TMath::Cos;
using TMath::Sin;
using TMath::Abs;

acceptanceFunction::acceptanceFunction(double x_min_,double x_max_,double y_min_,double y_max_,double z_Up_,double z_Down_,double bench_angle_){
	x_min = x_min_;
	x_max = x_max_;
	y_min = y_min_;
	y_max = y_max_;
	z_Up = z_Up_;
	z_Down = z_Down_;
	bench_angle = bench_angle_;
}

acceptanceFunction::~acceptanceFunction(){

}

double acceptanceFunction::operator()(double x,double y, double z){
	double theta_x_min = 0;
	double theta_x_max = 0;
	double theta_y_min = 0;
	double theta_y_max = 0;
	if(z<z_Up && z>z_Down){
		theta_x_min = Max((x-x_min)/(z-z_Up),(x_max-x)/(z_Down-z));
		theta_x_max = Min((x_max-x)/(z_Up-z),(x-x_min)/(z-z_Down));
		if(theta_x_max < theta_x_min) return 0;
		theta_y_min = Max((y-y_min)/(z-z_Up),(y_max-y)/(z_Down-z));
		theta_y_max = Min((y_max-y)/(z_Up-z),(y-y_min)/(z-z_Down));
		if(theta_y_max < theta_y_min) return 0;
	}
	else if(z>=z_Up){
		theta_x_max = Min((x-x_min)/(z-z_Up),(x-x_min)/(z-z_Down));
		theta_x_min = Max((x_max-x)/(z_Up-z),(x_max-x)/(z_Down-z));
		if(theta_x_max < theta_x_min) return 0;
		theta_y_max = Min((y-y_min)/(z-z_Up),(y-y_min)/(z-z_Down));
		theta_y_min = Max((y_max-y)/(z_Up-z),(y_max-y)/(z_Down-z));
		if(theta_y_max < theta_y_min) return 0;
	}
	else{
		theta_x_min = Max((x-x_min)/(z-z_Up),(x-x_min)/(z-z_Down));
		theta_x_max = Min((x_max-x)/(z_Up-z),(x_max-x)/(z_Down-z));
		if(theta_x_max < theta_x_min) return 0;
		theta_y_min = Max((y-y_min)/(z-z_Up),(y-y_min)/(z-z_Down));
		theta_y_max = Min((y_max-y)/(z_Up-z),(y_max-y)/(z_Down-z));
		if(theta_y_max < theta_y_min) return 0;
	}
	/*
	//For gaussian distribution
	theta_x_min = ATan(theta_x_min)/(sigma_theta*Sqrt(2));
	theta_x_max = ATan(theta_x_max)/(sigma_theta*Sqrt(2));
	theta_y_min = ATan(theta_y_min)/(sigma_theta*Sqrt(2));
	theta_y_max = ATan(theta_y_max)/(sigma_theta*Sqrt(2));
	return 0.25*(Erf(theta_x_max)- Erf(theta_x_min))*(Erf(theta_y_max)- Erf(theta_y_min));
	*/
	//For cos square distribution
	theta_x_min = ATan(theta_x_min);
	theta_x_max = ATan(theta_x_max);
	theta_y_min = ATan(theta_y_min);
	theta_y_max = ATan(theta_y_max);
	theta_y_min -= bench_angle;
	theta_y_max -= bench_angle;
	if(theta_x_min<-Pi()/2.) theta_x_min = -Pi()/2.;
	if(theta_x_max<-Pi()/2.) theta_x_max = -Pi()/2.;
	if(theta_x_min>Pi()/2.) theta_x_min = Pi()/2.;
	if(theta_x_max>Pi()/2.) theta_x_max = Pi()/2.;

	TF2 dist("dist","cos(sqrt(x*x+y*y))*cos(sqrt(x*x+y*y))",-Pi()/2.,Pi()/2.,-Pi()/2.,Pi()/2.);
	ROOT::Math::WrappedMultiTF1 wf1(dist);
	ROOT::Math::AdaptiveIntegratorMultiDim ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double thetaMin[] = {theta_x_min,theta_y_min};
	double thetaMax[] = {theta_x_max,theta_y_max};
	return ig.Integral(thetaMin,thetaMax);

	//return ((Sin(theta_x_max-theta_x_min)*Cos(theta_x_max+theta_x_min)) + theta_x_max - theta_x_min)*((Sin(theta_y_max-theta_y_min)*Cos(theta_y_max+theta_y_min)) + theta_y_max - theta_y_min)/(Pi()*Pi());
}

void acceptanceFunction::plot_3D(){
	int step_x = 1000;
	int step_y = 1000;
	int step_z = 1000;
	double x_min_ = -200;
	double x_max_ = 700;
	double y_min_ = -200;
	double y_max_ = 700;
	double z_Up_ = 1500;
	double z_Down_ = -200;
	TCanvas * cDisplay = new TCanvas();
	cDisplay->Divide(3);
	TH2D * proba_XY = new TH2D("proba_XY","proba_XY",step_x,x_min_,x_max_,step_y,y_min_,y_max_);
	TH2D * proba_XZ = new TH2D("proba_XZ","proba_XZ",step_x,x_min_,x_max_,step_z,z_Down_,z_Up_);
	TH2D * proba_YZ = new TH2D("proba_YZ","proba_YZ",step_y,y_min_,y_max_,step_z,z_Down_,z_Up_);
	proba_XY->SetStats(false);
	proba_XZ->SetStats(false);
	proba_YZ->SetStats(false);
	int n = 0;
	for(int i=0;i<step_x && Tomography::can_continue;i++){
		double x = x_min_ + i*(x_max_-x_min_)/step_x;
		for(int j=0;j<step_y && Tomography::can_continue;j++){
			double y = y_min_ + j*(y_max_-y_min_)/step_y;
			for(int k=0;k<step_z && Tomography::can_continue;k++){
				double z = z_Down_ + k*(z_Up_-z_Down_)/step_z;
				double proba = (*this)(x,y,z);
				proba_XY->Fill(x,y,proba);
				proba_XZ->Fill(x,z,proba);
				proba_YZ->Fill(y,z,proba);
				/*
				if((n%10000) == 0){
					cDisplay->cd(1);
					proba_XY->Draw("colz");
					cDisplay->cd(2);
					proba_XZ->Draw("colz");
					cDisplay->cd(3);
					proba_YZ->Draw("colz");
					cDisplay->Modified();
					cDisplay->Update();
				}
				if((n%1000) == 0) cout << "\r" << n << flush;
				n++;
				*/
			}
		}
	}
	proba_XY->Scale(1./step_z);
	proba_XZ->Scale(1./step_y);
	proba_YZ->Scale(1./step_x);
	cout << "\r" << n << endl;
	cDisplay->cd(1);
	proba_XY->Draw("colz");
	cDisplay->cd(2);
	proba_XZ->Draw("colz");
	cDisplay->cd(3);
	proba_YZ->Draw("colz");
	cDisplay->Modified();
	cDisplay->Update();

}

TH2D acceptanceFunction::plot_XY(int nbin_x,double x1,double x2,int nbin_y,double y1,double y2, double z, double y_angle){
	double x_min_plot = x1;
	double x_max_plot = x2;
	double y_min_plot = y1;
	double y_max_plot = y2;
	if(Abs(y_angle)>Pi()/2.){
		return TH2D("error","error",nbin_x,x_min_plot,x_max_plot,nbin_y,y_min_plot,y_max_plot);
	}
	int step_x = 20*nbin_x;
	int step_y = 20*nbin_y;
	TH2D proba_XY("proba_XY","proba_XY",nbin_x,x1,x2,nbin_y,y1,y2);
	proba_XY.SetStats(false);
	double width_x = x_max_plot - x_min_plot;
	double width_y = y_max_plot - y_min_plot;
	x_max_plot += width_x/(nbin_x*100.);
	x_min_plot += width_x/(nbin_x*100.);
	y_max_plot += width_y/(nbin_y*100.);
	y_min_plot += width_y/(nbin_y*100.);
	double width_step_x = (x_max_plot - x_min_plot)/step_x;
	double width_step_y = (y_max_plot - y_min_plot)/step_y;
	for(int i=0;i<step_x && Tomography::can_continue;i++){
		double x = x_min_plot + i*width_step_x;
		if(x<x1 || x>x2) continue;
		for(int j=0;j<step_y && Tomography::can_continue;j++){
			double y = y_min_plot + j*width_step_y;
			if(y<y1 || y>y2) continue;
			double real_x = x;
			double real_y = y*Cos(y_angle);
			double real_z = z - y*Sin(y_angle);
			double proba = (*this)(real_x,real_y,real_z);
			proba_XY.Fill(x,y,proba);
		}
	}
	return proba_XY;
}
TH2D acceptanceFunction::plot_XY(int nbin_x, int nbin_y,double z, double y_angle){
	int step_x = 20*nbin_x;
	int step_y = 20*nbin_y;
	double x_min_plot = x_min;
	double x_max_plot = x_max;
	double y_min_plot = y_min;
	double y_max_plot = y_max;
	if(Abs(y_angle)>Pi()/2.){
		return TH2D("error","error",nbin_x,x_min_plot,x_max_plot,nbin_y,y_min_plot,y_max_plot);
	}
	if(z>z_Up || z<z_Down){
		Point orig(0,0,z);
		Point norm(0,Sin(y_angle),Cos(y_angle));
		Plane proj(norm,orig);
		Line first_line(Point(x_min,y_min,z_Down),Point(x_max,y_max,z_Up));
		Line second_line(Point(x_max,y_max,z_Down),Point(x_min,y_min,z_Up));
		Point corner_a = proj.intersection(first_line);
		Point corner_b = proj.intersection(second_line);
		cout << corner_a << endl;
		cout << corner_b << endl;
		x_max_plot = Max(Abs((corner_a-orig).get_X()),Abs((corner_b-orig).get_X()));
		y_max_plot = Max((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
		x_min_plot = -x_max_plot;
		y_min_plot = Min((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
	}
	y_max_plot = y_max_plot/Cos(y_angle);
	y_min_plot = y_min_plot/Cos(y_angle);
	double width_x = x_max_plot - x_min_plot;
	double width_y = y_max_plot - y_min_plot;
	x_max_plot += 0.05*width_x;
	x_min_plot -= 0.05*width_x;
	y_max_plot += 0.05*width_y;
	y_min_plot -= 0.05*width_y;
	const double x_min_hist = x_min_plot;
	const double x_max_hist = x_max_plot;
	const double y_min_hist = y_min_plot;
	const double y_max_hist = y_max_plot;
	TH2D proba_XY("proba_XY","proba_XY",nbin_x,x_min_hist,x_max_hist,nbin_y,y_min_hist,y_max_hist);
	proba_XY.SetStats(false);
	width_x = x_max_plot - x_min_plot;
	width_y = y_max_plot - y_min_plot;
	x_max_plot += width_x/(nbin_x*100.);
	x_min_plot += width_x/(nbin_x*100.);
	y_max_plot += width_y/(nbin_y*100.);
	y_min_plot += width_y/(nbin_y*100.);
	double width_step_x = (x_max_plot - x_min_plot)/step_x;
	double width_step_y = (y_max_plot - y_min_plot)/step_y;
	for(int i=0;i<step_x && Tomography::can_continue;i++){
		double x = x_min_plot + i*width_step_x;
		if(x<x_min_hist || x>x_max_hist) continue;
		for(int j=0;j<step_y && Tomography::can_continue;j++){
			double y = y_min_plot + j*width_step_y;
			if(y<y_min_hist || y>y_max_hist) continue;
			double real_x = x;
			double real_y = y*Cos(y_angle);
			double real_z = z - y*Sin(y_angle);
			double proba = (*this)(real_x,real_y,real_z);
			proba_XY.Fill(x,y,proba);
		}
	}
	return proba_XY;
}