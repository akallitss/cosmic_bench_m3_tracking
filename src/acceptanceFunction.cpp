#define acceptanceFunction_cpp

#include "acceptanceFunction.h"

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

acceptanceFunction::acceptanceFunction(double x_min_,double x_max_,double y_min_,double y_max_,double z_Up_,double z_Down_,double sigma_theta_){
	x_min = x_min_;
	x_max = x_max_;
	y_min = y_min_;
	y_max = y_max_;
	z_Up = z_Up_;
	z_Down = z_Down_;
	sigma_theta = sigma_theta_;
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
	for(int i=0;i<step_x;i++){
		double x = x_min_ + i*(x_max_-x_min_)/step_x;
		for(int j=0;j<step_y;j++){
			double y = y_min_ + j*(y_max_-y_min_)/step_y;
			for(int k=0;k<step_z;k++){
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

TH2D acceptanceFunction::plot_XY(int nbin_x,double x1,double x2,int nbin_y,double y1,double y2, double z){
	int step_x = 20*nbin_x;
	int step_y = 20*nbin_y;
	TH2D proba_XY("proba_XY","proba_XY",nbin_x,x1,x2,nbin_y,y1,y2);
	proba_XY.SetStats(false);
	for(int i=0;i<step_x;i++){
		double x = x1 + i*(x2-x1)/step_x;
		for(int j=0;j<step_y;j++){
			double y = y1 + j*(y2-y1)/step_y;
			double proba = (*this)(x,y,z);
			proba_XY.Fill(x,y,proba);
		}
	}
	return proba_XY;
}