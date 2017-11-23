#define ray_cpp
#include "ray.h"
#include "cluster.h"
#include "detector.h"
#include "point.h"
#include "tomography.h"

#include <limits>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>


using std::cout;
using std::endl;
using std::numeric_limits;
using std::string;

using TMath::Cos;
using TMath::Sin;
using TMath::Sqrt;
using TMath::ATan;
using TMath::Abs;

Ray_2D::Ray_2D(){
	clusters.clear();
	chiSquare = -1;
	slope = 0;
	Z_intercept = 0;
	coord = 'Z';
}
Ray_2D::Ray_2D(char coord_){
	clusters.clear();
	if(coord_ != 'X' && coord_ != 'Y') return;
	coord = coord_;
	chiSquare = -1;
	slope = 0;
	Z_intercept = 0;
}
Ray_2D::Ray_2D(const Ray_2D& other){
	coord = other.coord;
	chiSquare = other.chiSquare;
	slope = other.slope;
	Z_intercept = other.Z_intercept;
	clusters.clear();
	for(vector<Cluster*>::const_iterator it = other.clusters.begin();it!=other.clusters.end();++it){
		clusters.push_back((*it)->Clone());
	}
}
Ray_2D& Ray_2D::operator=(const Ray_2D& other){
	coord = other.coord;
	chiSquare = other.chiSquare;
	slope = other.slope;
	Z_intercept = other.Z_intercept;
	clear();
	for(vector<Cluster*>::const_iterator it = other.clusters.begin();it!=other.clusters.end();++it){
		clusters.push_back((*it)->Clone());
	}
	return *this;
}
Ray_2D::Ray_2D(const Ray& other, char coord_){
	if(coord_ != 'X' && coord_ != 'Y') return;
	coord = coord_;
	if(coord == 'X'){
		slope = other.slope_X;
		Z_intercept = other.Z_intercept_X;
		chiSquare = other.chiSquare_X;
	}
	else if(coord == 'Y'){
		slope = other.slope_Y;
		Z_intercept = other.Z_intercept_Y;
		chiSquare = other.chiSquare_Y;
	}
	clusters.clear();
	for(vector<Cluster*>::const_iterator it = other.clusters.begin(); it!=other.clusters.end();++it){
		if(((*it)->get_is_X() && coord == 'X')||(coord == 'Y' && !(*it)->get_is_X())){
			clusters.push_back((*it)->Clone());
		}
	}
}
void Ray_2D::clear(){
	for(unsigned int i=0;i<clusters.size();i++){
		delete clusters[i];
	}
	clusters.clear();
}
bool Ray_2D::has_layer(int layer) const{
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		if((*it)->get_layer() == layer) return true;
	}
	return false;
}
Ray_2D::~Ray_2D(){
	clear();
	clusters.clear();
}
void Ray_2D::add_cluster(const Cluster * const clus){
	if(clus->get_type() == Tomography::CM){
		if(dynamic_cast<const CM_Cluster * const>(clus)->get_strip_type() != Tomography::Demux) return;
	}
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		if((clus->get_type() == (*it)->get_type()) && (clus->get_n_in_tree() == (*it)->get_n_in_tree())) return;
	}
	if((clus->get_is_X() && coord == 'X')||(coord == 'Y' && !clus->get_is_X())){
		clusters.push_back(clus->Clone());
	}
}
/*
void Ray_2D::process(){
	if(clusters.size()<2) return;
	TGraph * pos = new TGraph();
	double maxZ = numeric_limits<double>::min();
	double minZ = numeric_limits<double>::max();
	int i =0;
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		if((*it)->get_z()>maxZ) maxZ = (*it)->get_z();
		if((*it)->get_z()<minZ) minZ = (*it)->get_z();
		pos->SetPoint(i,(*it)->get_z(),(*it)->get_pos_mm());
		i++;
	}
	pos->Sort();
	TF1 * line = new TF1("line","pol1(0)",minZ-10,maxZ+10);
	//TF1 * line = new TF1("line","[0] + [1]*x",minZ-10,maxZ+10);
	double maxSlope = Tomography::get_instance()->get_XY_size()/(maxZ-minZ);
	line->SetParameters(0,0);
	
	if(minZ>0) line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.)-(maxSlope*minZ),(Tomography::get_instance()->get_XY_size()/2.)+maxSlope*minZ);
	else if(maxZ<0) line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.)+(maxSlope*maxZ),(Tomography::get_instance()->get_XY_size()/2.)-maxSlope*maxZ);
	else line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.),Tomography::get_instance()->get_XY_size()/2.);
	line->SetParLimits(1,-maxSlope,maxSlope);
	pos->Fit(line,"QN");
	chiSquare = line->GetChisquare();
	Z_intercept = line->GetParameter(0);
	slope = line->GetParameter(1);
	//if(line->Eval(1398.)>600. || line->Eval(1398.)<-100. || line->Eval(27.)>600. || line->Eval(27.)<-100.) chiSquare = numeric_limits<double>::max();
	delete pos; delete line;
}
*/
//MT safe implementation
/*
#include <Fit/BinData.h>
#include <Fit/Fitter.h>
#include <Math/WrappedMultiTF1.h>
#include <Math/WrappedParamFunction.h>
void Ray_2D::process(){
	if(clusters.size()<2) return;
	double maxZ = numeric_limits<double>::min();
	double minZ = numeric_limits<double>::max();
	ROOT::Fit::BinData pos(clusters.size(),1,ROOT::Fit::BinData::kNoError);
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		if((*it)->get_z()>maxZ) maxZ = (*it)->get_z();
		if((*it)->get_z()<minZ) minZ = (*it)->get_z();
		pos.Add((*it)->get_z(),(*it)->get_pos_mm());
	}
	TF1 * line = new TF1("line","pol1(0)",minZ-10,maxZ+10);
	double maxSlope = Tomography::get_instance()->get_XY_size()/(maxZ-minZ);
	line->SetParameters(0,0);
	if(minZ>0) line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.)-(maxSlope*minZ),(Tomography::get_instance()->get_XY_size()/2.)+maxSlope*minZ);
	else if(maxZ<0) line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.)+(maxSlope*maxZ),(Tomography::get_instance()->get_XY_size()/2.)-maxSlope*maxZ);
	else line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.),Tomography::get_instance()->get_XY_size()/2.);
	line->SetParLimits(1,-maxSlope,maxSlope);
	ROOT::Math::WrappedMultiTF1 wline(*line);
	ROOT::Math::IParamMultiFunction & pline = wline;
	ROOT::Fit::Fitter fitter;
	fitter.SetFunction(pline);
	fitter.Fit(pos,pline);
	chiSquare = fitter.Result().Chi2();
	Z_intercept = fitter.Result().Parameter(0);
	slope = fitter.Result().Parameter(1);
	delete line;
}
*/

//This implementation is two order of magnitude faster than the one above
void Ray_2D::process(){
	if(clusters.size()<2) return;
	double mean_x = 0;
	double mean_xx = 0;
	double mean_y = 0;
	double mean_xy = 0;
	int i = clusters.size();
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		mean_x += (*it)->get_z();
		mean_y +=  (*it)->get_pos_mm();
		mean_xy += ((*it)->get_z())*((*it)->get_pos_mm());
		mean_xx += ((*it)->get_z())*((*it)->get_z());
	}
	mean_x /= i;
	mean_y /= i;
	mean_xy /= i;
	mean_xx /= i;
	slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
	Z_intercept = mean_y - slope*mean_x;
	double distance = 0;
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		distance += (slope*((*it)->get_z()) - (*it)->get_pos_mm() + Z_intercept)*(slope*((*it)->get_z()) - (*it)->get_pos_mm() + Z_intercept);
	}
	distance /= (1+slope*slope);
	chiSquare = distance;
}

double Ray_2D::get_chiSquare() const{
	return chiSquare;
}
double Ray_2D::get_slope() const{
	return slope;
}
double Ray_2D::get_Z_intercept() const{
	return Z_intercept;
}
double Ray_2D::eval(double z) const{
	return slope*z+Z_intercept;
}
double Ray_2D::get_residu(const Detector * const det) const{
	if(clusters.size()<2) return numeric_limits<double>::min();
	TGraph * pos = new TGraph();
	double maxZ = numeric_limits<double>::min();
	double minZ = numeric_limits<double>::max();
	int i =0;
	double det_coord = 0;
	bool is_in_ray = false;
	for(vector<Cluster*>::const_iterator it = clusters.begin(); it!=clusters.end();++it){
		if((*it)->get_z()>maxZ) maxZ = (*it)->get_z();
		if((*it)->get_z()<minZ) minZ = (*it)->get_z();
		if(!((*it)->is_in_det(det))){
			pos->SetPoint(i,(*it)->get_z(),(*it)->get_pos_mm());
			i++;
		}
		else{
			det_coord = (*it)->get_pos_mm();
			is_in_ray = true;
		}
	}
	if(!is_in_ray){
		delete pos;
		return numeric_limits<double>::min();
	}
	TF1 * line = new TF1("line","[0]*x+[1]",minZ-10,maxZ+10);
	double maxSlope = Tomography::get_instance()->get_XY_size()/(maxZ-minZ);
	line->SetParameters(0,0);
	line->SetParLimits(0,-(Tomography::get_instance()->get_XY_size()/2.)-(maxSlope*minZ),(Tomography::get_instance()->get_XY_size()/2.)+maxSlope*minZ);
	line->SetParLimits(1,-maxSlope,maxSlope);
	pos->Fit(line,"QN");
	double residu = det_coord - line->Eval(det->get_z());
	delete pos; delete line;
	return residu;
}
double Ray_2D::get_residu_ref(const Cluster * const clus) const{
	bool is_in_ray = false;
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		if((*it)->get_type() == clus->get_type() && (*it)->get_n_in_tree() == clus->get_n_in_tree()){
			is_in_ray = true;
			break;
		}
	}
	if(is_in_ray) return numeric_limits<double>::max();
	return clus->get_pos_mm() - (Z_intercept + slope*clus->get_z());
}
double Ray_2D::get_t_mean() const{
	double t = 0;
	double size = clusters.size();
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		t+=(*it)->get_t();
	}
	t/=size;
	return t;
}
double Ray_2D::get_t_sigma() const{
	double sigma = 0;
	double size = clusters.size();
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		sigma+=((*it)->get_t())*((*it)->get_t());
	}
	sigma/=size;
	double mean = get_t_mean();
	sigma-=mean*mean;
	return Sqrt(sigma);
}
unsigned int Ray_2D::get_clus_n() const{
	return clusters.size();
}
vector<Cluster*> Ray_2D::get_clus() const{
	vector<Cluster*> return_vector;
	for(vector<Cluster*>::const_iterator it = clusters.begin();it != clusters.end(); ++it){
		return_vector.push_back((*it)->Clone());
	}
	return return_vector;
}
pair<int,int> Ray_2D::get_extremal_det(const CosmicBench * const bench) const{
	pair<int,int> return_pair(-1,-1);
	pair<double,double> extremal_z(numeric_limits<double>::max(), numeric_limits<double>::min());
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		double current_z = (*it)->get_z();
		if(current_z < extremal_z.first){
			extremal_z.first = current_z;
			return_pair.first = bench->find_det(*it);
		}
		if(current_z > extremal_z.second){
			extremal_z.second = current_z;
			return_pair.second = bench->find_det(*it);
		}
	}
	return return_pair;
}
ostream& operator<<(ostream& os, const Ray_2D& ray){
	os << "RAY_" << ray.coord << "(slope : " << ray.slope << " ; intercept : " << ray.Z_intercept << " ; ChiSquare : " << ray.chiSquare << " ; NClus : " << ray.clusters.size() << ")";
	return os;
}

Ray::Ray(){
	clusters.clear();
	chiSquare_X = -1;
	slope_X = 0;
	Z_intercept_X = 0;
	chiSquare_X = -1;
	slope_X = 0;
	Z_intercept_X = 0;
}
Ray::Ray(const Ray& other){
	chiSquare_X = other.chiSquare_X;
	chiSquare_Y = other.chiSquare_Y;
	slope_X = other.slope_X;
	slope_Y = other.slope_Y;
	Z_intercept_X = other.Z_intercept_X;
	Z_intercept_Y = other.Z_intercept_Y;
	clusters.clear();
	for(vector<Cluster*>::const_iterator it = other.clusters.begin(); it!=other.clusters.end();++it){
		clusters.push_back((*it)->Clone());
	}
}
Ray& Ray::operator=(const Ray& other){
	chiSquare_X = other.chiSquare_X;
	chiSquare_Y = other.chiSquare_Y;
	slope_X = other.slope_X;
	slope_Y = other.slope_Y;
	Z_intercept_X = other.Z_intercept_X;
	Z_intercept_Y = other.Z_intercept_Y;
	clear();
	for(vector<Cluster*>::const_iterator it = other.clusters.begin(); it!=other.clusters.end();++it){
		clusters.push_back((*it)->Clone());
	}
	return *this;
}
Ray::Ray(const Ray_2D& ray1, const Ray_2D& ray2){
	if(ray1.coord == ray2.coord) return;
	Ray_2D * rayX = (ray1.coord == 'X') ? new Ray_2D(ray1) : new Ray_2D(ray2);
	Ray_2D * rayY = (ray1.coord == 'Y') ? new Ray_2D(ray1) : new Ray_2D(ray2);
	chiSquare_X = rayX->get_chiSquare();
	chiSquare_Y = rayY->get_chiSquare();
	slope_X = rayX->get_slope();
	slope_Y = rayY->get_slope();
	Z_intercept_X = rayX->get_Z_intercept();
	Z_intercept_Y = rayY->get_Z_intercept();
	clusters.clear();
	for(vector<Cluster*>::const_iterator it = ray1.clusters.begin(); it!=ray1.clusters.end();++it){
		clusters.push_back((*it)->Clone());
    }
    for(vector<Cluster*>::const_iterator it = ray2.clusters.begin(); it!=ray2.clusters.end();++it){
    	clusters.push_back((*it)->Clone());
    }
	delete rayX; delete rayY;
}
void Ray::clear(){
	for(unsigned int i=0;i<clusters.size();i++){
		delete clusters[i];
	}
	clusters.clear();
}
Ray::~Ray(){
	clear();
	clusters.clear();
}
void Ray::add_cluster(const Cluster * const clus){
	if(clus->get_type() == Tomography::CM){
		if(dynamic_cast<const CM_Cluster * const>(clus)->get_strip_type() != Tomography::Demux) return;
	}
	for(vector<Cluster*>::iterator it = clusters.begin(); it!=clusters.end();++it){
		if((clus->get_type() == (*it)->get_type()) && (clus->get_n_in_tree() == (*it)->get_n_in_tree())) return;
	}
	clusters.push_back(clus->Clone());
}
void Ray::process(){
	if(clusters.size()<4) return;
	Ray_2D * rayX = new Ray_2D(*this,'X');
	Ray_2D * rayY = new Ray_2D(*this,'Y');
	rayX->process();
	rayY->process();
	clear();
	*this = Ray(*rayX,*rayY);
	delete rayX; delete rayY;
}
void Ray::angle_correction(){
	double delta_param = 100;
	double old_delta_param = 100;
	while(delta_param>0.1 && delta_param<old_delta_param){
		double current_Z_intercept_X = Z_intercept_X;
		double current_Z_intercept_Y = Z_intercept_Y;
		double current_slope_X = slope_X;
		double current_slope_Y = slope_Y;

		for(vector<Cluster*>::iterator it = clusters.begin();it!=clusters.end();++it){
			(*it)->set_perp_pos_mm(*this);
		}
		this->process();
		old_delta_param = delta_param;
		delta_param = Abs(Z_intercept_X - current_Z_intercept_X)+Abs(Z_intercept_Y - current_Z_intercept_Y)+Abs(slope_X - current_slope_X)+Abs(slope_Y - current_slope_Y);
	}
	for(vector<Cluster*>::iterator it = clusters.begin();it!=clusters.end();++it){
		if((*it)->get_is_X()){
			(*it)->set_perp_pos_mm(this->eval_Y((*it)->get_z()));
		}
		else{
			(*it)->set_perp_pos_mm(this->eval_X((*it)->get_z()));
		}
	}
	this->process();
}
double Ray::get_chiSquare_X() const{
	return chiSquare_X;
}
double Ray::get_chiSquare_Y() const{
	return chiSquare_Y;
}
double Ray::get_slope_X() const{
	return slope_X;
}
double Ray::get_slope_Y() const{
	return slope_Y;
}
double Ray::get_Z_intercept_X() const{
	return Z_intercept_X;
}
double Ray::get_Z_intercept_Y() const{
	return Z_intercept_Y;
}
double Ray::eval_X(double z) const{
	return slope_X*z+Z_intercept_X;
}
double Ray::eval_Y(double z) const{
	return slope_Y*z+Z_intercept_Y;
}
double Ray::eval_X(const Detector * const det) const{
	double n_x = Cos(det->get_angle_z())*Sin(det->get_angle_y()) + Sin(det->get_angle_z())*Sin(det->get_angle_x())*Cos(det->get_angle_y());
	double n_y = Cos(det->get_angle_z())*Sin(det->get_angle_x())*Cos(det->get_angle_y()) - Sin(det->get_angle_z())*Sin(det->get_angle_y());
	double n_z = Cos(det->get_angle_x())*Cos(det->get_angle_y());
	double z = (n_z*det->get_z() - Z_intercept_X*n_x - Z_intercept_Y*n_y)/(n_z + slope_X*n_x + slope_Y*n_y);
	return slope_X*z + Z_intercept_X;
}
double Ray::eval_Y(const Detector * const det) const{
	double n_x = Cos(det->get_angle_z())*Sin(det->get_angle_y()) + Sin(det->get_angle_z())*Sin(det->get_angle_x())*Cos(det->get_angle_y());
	double n_y = Cos(det->get_angle_z())*Sin(det->get_angle_x())*Cos(det->get_angle_y()) - Sin(det->get_angle_z())*Sin(det->get_angle_y());
	double n_z = Cos(det->get_angle_x())*Cos(det->get_angle_y());
	double z = (n_z*det->get_z() - Z_intercept_X*n_x - Z_intercept_Y*n_y)/(n_z + slope_X*n_x + slope_Y*n_y);
	return slope_Y*z + Z_intercept_Y;
}
Point Ray::eval_plane(Plane proj) const{
	Point first(Z_intercept_X,Z_intercept_Y,0);
	Point second(Z_intercept_X + slope_X,Z_intercept_Y + slope_Y,1); // <-- this is not a direction vector, this is a second point defining the line
	Line ray_line(first,second);
	return proj.intersection(ray_line);
}
double Ray::get_residu(const Detector * const det) const{
	bool is_in_ray = false;
	for(vector<Cluster*>::const_iterator it = clusters.begin(); it!=clusters.end();++it){
		if((*it)->is_in_det(det)){
			is_in_ray = true;
		}
	}
	if(!is_in_ray){
		return numeric_limits<double>::min();
	}
	else{
		bool det_is_x = det->get_is_X();
		if(det_is_x){
			Ray_2D sub_ray = Ray_2D(*this,'X');
			return sub_ray.get_residu(det);
		}
		else{
			Ray_2D sub_ray = Ray_2D(*this,'Y');
			return sub_ray.get_residu(det);
		}
	}
}
double Ray::get_residu_ref(const Cluster * const clus) const{
	if(clus->get_is_X()){
		Ray_2D sub_ray = Ray_2D(*this,'X');
		return sub_ray.get_residu_ref(clus);
	}
	else{
		Ray_2D sub_ray = Ray_2D(*this,'Y');
		return sub_ray.get_residu_ref(clus);
	}
}
double Ray::get_t_mean() const{
	double t = 0;
	double size = clusters.size();
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		t+=(*it)->get_t();
	}
	t/=size;
	return t;
}
double Ray::get_t_sigma() const{
	double sigma = 0;
	double size = clusters.size();
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		sigma+=((*it)->get_t())*((*it)->get_t());
	}
	sigma/=size;
	double mean = get_t_mean();
	sigma-=mean*mean;
	return Sqrt(sigma);
}
unsigned int Ray::get_clus_n() const{
	return clusters.size();
}
unsigned int Ray::get_clus_x_n() const{
	int n = 0;
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		if((*it)->get_is_X()) n++;
	}
	return n;
}
unsigned int Ray::get_clus_y_n() const{
	return (get_clus_n() - get_clus_x_n());
}
vector<Cluster*> Ray::get_clus() const{
	vector<Cluster*> return_vector;
	for(vector<Cluster*>::const_iterator it = clusters.begin();it != clusters.end(); ++it){
		return_vector.push_back((*it)->Clone());
	}
	return return_vector;
}
pair<pair<int,int>,pair<int,int> > Ray::get_extremal_det(const CosmicBench * const bench) const{
	pair<pair<int,int>,pair<int,int> > return_pair(pair<int,int>(-1,-1),pair<int,int>(-1,-1));
	pair<double,double> extremal_z(numeric_limits<double>::max(), numeric_limits<double>::min());
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		double current_z = (*it)->get_z();
		if((*it)->get_is_X()){
			if(current_z < extremal_z.first){
				extremal_z.first = current_z;
				return_pair.first.first = bench->find_det(*it);
				return_pair.second.first = bench->find_det(bench->get_detector(return_pair.first.first)->get_perp_type(),bench->get_detector(return_pair.first.first)->get_perp_n());
			}
			if(current_z > extremal_z.second){
				extremal_z.second = current_z;
				return_pair.first.second = bench->find_det(*it);
				return_pair.second.second = bench->find_det(bench->get_detector(return_pair.first.second)->get_perp_type(),bench->get_detector(return_pair.first.second)->get_perp_n());
			}
		}
		else{
			if(current_z < extremal_z.first){
				extremal_z.first = current_z;
				return_pair.second.first = bench->find_det(*it);
				return_pair.first.first = bench->find_det(bench->get_detector(return_pair.second.first)->get_perp_type(),bench->get_detector(return_pair.second.first)->get_perp_n());
			}
			if(current_z > extremal_z.second){
				extremal_z.second = current_z;
				return_pair.second.second = bench->find_det(*it);
				return_pair.first.second = bench->find_det(bench->get_detector(return_pair.second.second)->get_perp_type(),bench->get_detector(return_pair.second.second)->get_perp_n());
			}
		}
	}
	return return_pair;
}
bool Ray::has_layer(int layer) const{
	for(vector<Cluster*>::const_iterator it = clusters.begin();it!=clusters.end();++it){
		if((*it)->get_layer() == layer) return true;
	}
	return false;
}

ostream& operator<<(ostream& os, const Ray& ray){
	os << "[" << Ray_2D(ray,'X') << " ; " << Ray_2D(ray,'Y') << "]";
	return os;
}

RayPair::RayPair(){
	downRay = Ray();
	upRay = Ray();
	delta_x = 0;
	delta_y = 0;
	delta_theta_x = 0;
	delta_theta_y = 0;
}
RayPair::RayPair(const RayPair& other){
	delta_x = other.delta_x;
	delta_y = other.delta_y;
	delta_theta_x = other.delta_theta_y;
	delta_theta_y = other.delta_theta_y;
	downRay = other.downRay;
	upRay = other.upRay;
}
RayPair& RayPair::operator=(const RayPair& other){
	delta_x = other.delta_x;
	delta_y = other.delta_y;
	delta_theta_x = other.delta_theta_y;
	delta_theta_y = other.delta_theta_y;
	downRay = other.downRay;
	upRay = other.upRay;
	return *this;
}
RayPair::RayPair(const Ray& ray1, const Ray& ray2){
	downRay = Ray();
	upRay = Ray();
	delta_x = 0;
	delta_y = 0;
	delta_theta_x = 0;
	delta_theta_y = 0;
	bool is_1_up = Tomography::get_instance()->get_is_up((*(ray1.clusters.begin()))->get_layer());
	for(vector<Cluster*>::const_iterator it = ray1.clusters.begin()+1;it!=ray1.clusters.end();++it){
		if(Tomography::get_instance()->get_is_up((*it)->get_layer()) != is_1_up) return;
	}
	bool is_2_up = Tomography::get_instance()->get_is_up((*(ray2.clusters.begin()))->get_layer());
	for(vector<Cluster*>::const_iterator it = ray2.clusters.begin()+1;it!=ray2.clusters.end();++it){
		if(Tomography::get_instance()->get_is_up((*it)->get_layer()) != is_2_up) return;
	}
	if(!(is_1_up != is_2_up)) return;
	upRay = (is_1_up) ? Ray(ray1) : Ray(ray2);
	downRay = (is_2_up) ? Ray(ray1) : Ray(ray2);
}
RayPair::~RayPair(){
	
}
void RayPair::process(){
	//upRay.process();
	//downRay.process();
	double theta_up_x = ATan(upRay.slope_X);
	double theta_up_y = ATan(upRay.slope_Y);
	double theta_down_x = ATan(downRay.slope_X);
	double theta_down_y = ATan(downRay.slope_Y);
	delta_theta_x = theta_down_x - theta_up_x;
	delta_theta_y = theta_down_y - theta_up_y;
	double z_up = numeric_limits<double>::max();
	double z_down = numeric_limits<double>::min();
	for(vector<Cluster*>::const_iterator it = upRay.clusters.begin();it!=upRay.clusters.end();++it){
		if((*it)->get_z()<z_up) z_up = (*it)->get_z();
	}
	for(vector<Cluster*>::const_iterator it = downRay.clusters.begin();it!=downRay.clusters.end();++it){
		if((*it)->get_z()>z_down) z_down = (*it)->get_z();
	}
	double L_xy = (z_up-z_down)*Sqrt(1+(upRay.slope_X*upRay.slope_X)+(upRay.slope_Y*upRay.slope_Y));
	delta_x = (downRay.eval_X(z_down)-upRay.eval_X(z_down))*Cos(theta_up_x)*L_xy*Cos(delta_theta_x+theta_up_x)/Cos(delta_theta_x);
	delta_y = (downRay.eval_Y(z_down)-upRay.eval_Y(z_down))*Cos(theta_up_y)*L_xy*Cos(delta_theta_y+theta_up_y)/Cos(delta_theta_y);
}
double RayPair::get_delta_x() const{
	return delta_x;
}
double RayPair::get_delta_y() const{
	return delta_y;
}
double RayPair::get_delta_theta_x() const{
	return delta_theta_x;
}
double RayPair::get_delta_theta_y() const{
	return delta_theta_y;
}
double RayPair::get_theta_x_up() const{
	return ATan(upRay.slope_X);
}
double RayPair::get_theta_y_up() const{
	return ATan(upRay.slope_Y);
}
double RayPair::get_theta_x_down() const{
	return ATan(downRay.slope_X);
}
double RayPair::get_theta_y_down() const{
	return ATan(downRay.slope_Y);
}
double RayPair::get_x_up(double z) const{
	return upRay.eval_X(z);
}
double RayPair::get_y_up(double z) const{
	return upRay.eval_Y(z);
}
double RayPair::get_x_down(double z) const{
	return downRay.eval_X(z);
}
double RayPair::get_y_down(double z) const{
	return downRay.eval_Y(z);
}
double RayPair::get_doca(){
	//this->process();
	double perp_x = upRay.slope_Y-downRay.slope_X;
	double perp_y = downRay.slope_X-upRay.slope_Y;
	double perp_z = upRay.slope_X*downRay.slope_Y - upRay.slope_Y*downRay.slope_X;
	double norm = Sqrt(perp_x*perp_x+perp_y*perp_y+perp_z*perp_z);
	return Abs(perp_x*(upRay.Z_intercept_X-downRay.Z_intercept_X)+perp_y*(upRay.Z_intercept_Y-downRay.Z_intercept_Y))/norm;
}
Point RayPair::get_PoCA(){
	Point n_UpRay(upRay.slope_X,upRay.slope_Y,1);
	Point a_UpRay(upRay.Z_intercept_X,upRay.Z_intercept_Y,0);
	Point n_DownRay(downRay.slope_X,downRay.slope_Y,1);
	Point a_DownRay(downRay.Z_intercept_X,downRay.Z_intercept_Y,0);

	Point common_term = (a_DownRay-a_UpRay)/(n_UpRay.normSquare()*n_DownRay.normSquare()-scalar_product(n_UpRay,n_DownRay)*scalar_product(n_UpRay,n_DownRay));

	double param_closest_point_UpRay = scalar_product(common_term,(n_UpRay*n_DownRay.normSquare()-n_DownRay*scalar_product(n_UpRay,n_DownRay)));
	double param_closest_point_DownRay = -1*scalar_product(common_term,(n_DownRay*n_UpRay.normSquare()-n_UpRay*scalar_product(n_DownRay,n_UpRay)));

	Point closest_point_UpRay = param_closest_point_UpRay*n_UpRay + a_UpRay;
	Point closest_point_DownRay = param_closest_point_DownRay*n_DownRay + a_DownRay;
	return (closest_point_DownRay+closest_point_UpRay)*0.5;
}
