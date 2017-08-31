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
#include <utility>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::pair;
using std::setprecision;
using std::fixed;
using TMath::ATan;
using TMath::Tan;
using TMath::Erf;
using TMath::Sqrt;
using TMath::Max;
using TMath::Min;
using TMath::Pi;
using TMath::PiOver2;
using TMath::Cos;
using TMath::Sin;
using TMath::Abs;
using TMath::Binomial;

acceptanceFunction::acceptanceFunction(double x_min_Up_,double x_max_Up_,double y_min_Up_,double y_max_Up_,double x_min_Down_,double x_max_Down_,double y_min_Down_,double y_max_Down_,double z_Up_x_,double z_Down_x_,double z_Up_y_,double z_Down_y_,double bench_angle_){
	x_min_Up = x_min_Up_;
	x_max_Up = x_max_Up_;
	y_min_Up = y_min_Up_;
	y_max_Up = y_max_Up_;
	x_min_Down = x_min_Down_;
	x_max_Down = x_max_Down_;
	y_min_Down = y_min_Down_;
	y_max_Down = y_max_Down_;
	z_Up_x = z_Up_x_;
	z_Down_x = z_Down_x_;
	z_Up_y = z_Up_y_;
	z_Down_y = z_Down_y_;
	bench_angle = bench_angle_;
}

acceptanceFunction::~acceptanceFunction(){

}

double acceptanceFunction::operator()(double x,double y, double z){
	double theta_x_min = 0;
	double theta_x_max = 0;
	if(z<z_Up_x && z>z_Down_x){
		theta_x_min = Max((x-x_min_Up)/(z-z_Up_x),(x_max_Down-x)/(z_Down_x-z));
		theta_x_max = Min((x_max_Up-x)/(z_Up_x-z),(x-x_min_Down)/(z-z_Down_x));
		if(theta_x_max < theta_x_min) return 0;
	}
	else if(z>=z_Up_x){
		theta_x_max = Min((x-x_min_Up)/(z-z_Up_x),(x-x_min_Down)/(z-z_Down_x));
		theta_x_min = Max((x_max_Up-x)/(z_Up_x-z),(x_max_Down-x)/(z_Down_x-z));
		if(theta_x_max < theta_x_min) return 0;
	}
	else{
		theta_x_min = Max((x-x_min_Up)/(z-z_Up_x),(x-x_min_Down)/(z-z_Down_x));
		theta_x_max = Min((x_max_Up-x)/(z_Up_x-z),(x_max_Down-x)/(z_Down_x-z));
		if(theta_x_max < theta_x_min) return 0;
	}
	double theta_y_min = 0;
	double theta_y_max = 0;
	if(z<z_Up_y && z>z_Down_y){
		theta_y_min = Max((y-y_min_Up)/(z-z_Up_y),(y_max_Down-y)/(z_Down_y-z));
		theta_y_max = Min((y_max_Up-y)/(z_Up_y-z),(y-y_min_Down)/(z-z_Down_y));
		if(theta_y_max < theta_y_min) return 0;
	}
	else if(z>=z_Up_y){
		theta_y_max = Min((y-y_min_Up)/(z-z_Up_y),(y-y_min_Down)/(z-z_Down_y));
		theta_y_min = Max((y_max_Up-y)/(z_Up_y-z),(y_max_Down-y)/(z_Down_y-z));
		if(theta_y_max < theta_y_min) return 0;
	}
	else{
		theta_y_min = Max((y-y_min_Up)/(z-z_Up_y),(y-y_min_Down)/(z-z_Down_y));
		theta_y_max = Min((y_max_Up-y)/(z_Up_y-z),(y_max_Down-y)/(z_Down_y-z));
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
	for(int i=0;i<step_x && Tomography::get_instance()->get_can_continue();i++){
		double x = x_min_ + i*(x_max_-x_min_)/step_x;
		for(int j=0;j<step_y && Tomography::get_instance()->get_can_continue();j++){
			double y = y_min_ + j*(y_max_-y_min_)/step_y;
			for(int k=0;k<step_z && Tomography::get_instance()->get_can_continue();k++){
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
	x_min_plot -= width_x/(nbin_x*100.);
	y_max_plot += width_y/(nbin_y*100.);
	y_min_plot -= width_y/(nbin_y*100.);
	double width_step_x = (x_max_plot - x_min_plot)/step_x;
	double width_step_y = (y_max_plot - y_min_plot)/step_y;
	cout << "computing acceptance..." << endl;
	for(int i=0;i<step_x && Tomography::get_instance()->get_can_continue();i++){
		double x = x_min_plot + (i+0.5)*width_step_x;
		if(x<x1 || x>x2) continue;
		for(int j=0;j<step_y && Tomography::get_instance()->get_can_continue();j++){
			double y = y_min_plot + (j+0.5)*width_step_y;
			if(y<y1 || y>y2) continue;
			double real_x = x;
			double real_y = y*Cos(y_angle);
			double real_z = z - y*Sin(y_angle);
			double proba = (*this)(real_x,real_y,real_z);
			proba_XY.Fill(x,y,proba);
			cout << "\r" << setprecision(2) << fixed << (i*step_y+j+1)*100./(step_x*step_y) << "%" << flush;
		}
	}
	cout << "done !" << endl;
	return proba_XY;
}
/*
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
	x_min_plot -= width_x/(nbin_x*100.);
	y_max_plot += width_y/(nbin_y*100.);
	y_min_plot -= width_y/(nbin_y*100.);
	double width_step_x = (x_max_plot - x_min_plot)/step_x;
	double width_step_y = (y_max_plot - y_min_plot)/step_y;
	for(int i=0;i<step_x && Tomography::get_instance()->get_can_continue();i++){
		double x = x_min_plot + (i+0.5)*width_step_x;
		if(x<x_min_hist || x>x_max_hist) continue;
		for(int j=0;j<step_y && Tomography::get_instance()->get_can_continue();j++){
			double y = y_min_plot + (j+0.5)*width_step_y;
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
*/

FreeSkyFunction::FreeSkyFunction(double x_min_,double x_max_,double y_min_,double y_max_,vector<double> z_, double bench_angle_){
	x_min = x_min_;
	x_max = x_max_;
	y_min = y_min_;
	y_max = y_max_;
	bench_angle = bench_angle_;
	z = set<double>(z_.begin(),z_.end());
}
FreeSkyFunction::~FreeSkyFunction(){

}
void FreeSkyFunction::plot_3D(){

}
TH2D FreeSkyFunction::plot_PhiTheta(int nbin_phi,double phi1,double phi2,int nbin_theta,double theta1,double theta2,unsigned int mult){
	double phi_min_plot = phi1;
	double phi_max_plot = phi2;
	double theta_min_plot = theta1;
	double theta_max_plot = theta2;
	int step_phi = 20*nbin_phi;
	int step_theta = 20*nbin_theta;
	TH2D proba_PhiTheta("proba_PhiTheta","proba_PhiTheta",nbin_phi,phi1,phi2,nbin_theta,theta1,theta2);
	proba_PhiTheta.SetStats(false);
	double width_phi = phi_max_plot - phi_min_plot;
	double width_theta = theta_max_plot - theta_min_plot;
	phi_max_plot += width_phi/(nbin_phi*100.);
	phi_min_plot -= width_phi/(nbin_phi*100.);
	theta_max_plot += width_theta/(nbin_theta*100.);
	theta_min_plot -= width_theta/(nbin_theta*100.);
	double width_step_phi = (phi_max_plot - phi_min_plot)/step_phi;
	double width_step_theta = (theta_max_plot - theta_min_plot)/step_theta;
	map<double,int> delta_z;
	for(set<double>::iterator z_it=z.begin();z_it!=z.end();++z_it){
		set<double>::iterator z_jt = z.end();
		--z_jt;
		while(z_jt!=z_it){
			if(distance(z_it,z_jt)>(mult-2)){
				delta_z.insert(pair<double,int>(*z_jt - *z_it,Binomial(distance(z_it,z_jt)-1,mult-2)));
			}
			--z_jt;
		}
	}
	for(int i=0;i<step_phi && Tomography::get_instance()->get_can_continue();i++){
		double phi = phi_min_plot + (i+0.5)*width_step_phi;
		if(phi<phi1 || phi>phi2) continue;
		for(int j=0;j<step_theta && Tomography::get_instance()->get_can_continue();j++){
			double theta = theta_min_plot + (j+0.5)*width_step_theta;
			if(theta<theta1 || theta>theta2) continue;
			double proba = 0;
			for(map<double,int>::iterator delta_z_it=delta_z.begin();delta_z_it!=delta_z.end();++delta_z_it){
				proba += ((*this)(phi,theta,delta_z_it->first))*(delta_z_it->second);
			}
			proba_PhiTheta.Fill(phi,theta,proba);
		}
	}
	return proba_PhiTheta;
}
TH2D FreeSkyFunction::plot_PhiTheta(int nbin_phi, int nbin_theta, unsigned int mult){

}
double FreeSkyFunction::operator()(double phi, double theta, double delta_z){
	double eff_x = ((x_max - x_min) - delta_z*Abs(Tan(phi)))/(x_max - x_min);
	double eff_y = ((y_max - y_min) - delta_z*Tan(Abs(theta - bench_angle)))/(y_max - y_min);
	if(eff_x<0 || eff_y<0) return 0;
	double flux_comp = Cos(PiOver2() - theta);
	flux_comp *= flux_comp;
	return eff_x*eff_y*flux_comp;
}
