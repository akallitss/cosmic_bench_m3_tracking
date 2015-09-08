#define cluster_cpp
#include "cluster.h"
#include "Tanalyse_R.h"
#include "detector.h"

#include <TMath.h>
#include <iostream>

using std::cout;
using std::endl;

using TMath::Min;
using TMath::Tan;
using TMath::Cos;
using TMath::Sin;

Cluster::Cluster(){
	evn = -1;
	number = -1;
	ampl = -1;
	size = -1;
	pos = -1;
	maxStripAmpl = -1;
	maxSample = -1;
	maxStrip = -1;
	TOT = -1;
	t = -1;
	z = -1;
	is_X = -1;
	is_up = -1;
	offset = 0;
	direction = false;
	perp_pos_mm = -1;
	angle_x = 0;
	angle_y = 0;
	angle_z = 0;
	n_in_tree = -1;
}
Cluster::Cluster(const Cluster& other){
	type = other.type;
	evn = other.evn;
	number = other.number;
	ampl = other.ampl;
	size = other.size;
	pos = other.pos;
	maxStripAmpl = other.maxStripAmpl;
	maxSample= other.maxSample;
	maxStrip = other.maxStrip;
	TOT = other.TOT;
	t = other.t;
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	offset = other.offset;
	direction = other.direction;
	perp_pos_mm = other.perp_pos_mm;
	angle_x = other.angle_x;
	angle_y = other.angle_y;
	angle_z = other.angle_z;
	n_in_tree = other.n_in_tree;
}
Cluster& Cluster::operator=(const Cluster& other){
	type = other.type;
	evn = other.evn;
	number = other.number;
	ampl = other.ampl;
	size = other.size;
	pos = other.pos;
	maxStripAmpl = other.maxStripAmpl;
	maxSample= other.maxSample;
	maxStrip = other.maxStrip;
	TOT = other.TOT;
	t = other.t;
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	offset = other.offset;
	direction = other.direction;
	perp_pos_mm = other.perp_pos_mm;
	angle_x = other.angle_x;
	angle_y = other.angle_y;
	angle_z = other.angle_z;
	n_in_tree = other.n_in_tree;
	return *this;
}
Cluster::Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	number = number_;
	ampl = -1;
	size = -1;
	pos = -1;
	maxStripAmpl = -1;
	maxSample = -1;
	maxStrip = -1;
	TOT = -1;
	t = -1;
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	perp_pos_mm = -1;
	angle_x = det->get_angle_x();
	angle_y = det->get_angle_y();
	angle_z = det->get_angle_z();
	n_in_tree = det->get_n_in_tree();
}
Cluster::Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_){
	evn = -1;
	number = number_;
	ampl = ampl_;
	size = size_;
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	perp_pos_mm = -1;
	angle_x = det->get_angle_x();
	angle_y = det->get_angle_y();
	angle_z = det->get_angle_z();
	pos = pos_;
	maxSample = maxSample_;
	maxStripAmpl = maxStripAmpl_;
	maxStrip = maxStrip_;
	TOT = TOT_;
	t = t_;
	n_in_tree = det->get_n_in_tree();
}
Tomography::det_type Cluster::get_type() const{
	return type;
}
bool Cluster::get_is_X() const{
	return is_X;
}
double Cluster::get_ampl() const{
	return ampl;
}
double Cluster::get_size() const{
	return size;
}
double Cluster::get_pos() const{
	return pos;
}
double Cluster::get_TOT() const{
	return TOT;
}
double Cluster::get_t() const{
	return t;
}
double Cluster::get_maxSample() const{
	return maxSample;
}
double Cluster::get_maxStripAmpl() const{
	return maxStripAmpl;
}
int Cluster::get_maxStrip() const{
	return maxStrip;
}
bool Cluster::get_is_up() const{
	return is_up;
}
int Cluster::get_n_in_tree() const{
	return n_in_tree;
}
void Cluster::set_perp_pos_mm(double coord){
	perp_pos_mm = coord;
}
void Cluster::set_perp_pos_mm(Ray ray){
	double n_x = Cos(angle_z)*Sin(angle_y) + Sin(angle_z)*Sin(angle_x)*Cos(angle_y);
	double n_y = Cos(angle_z)*Sin(angle_x)*Cos(angle_y) - Sin(angle_z)*Sin(angle_y);
	double n_z = Cos(angle_x)*Cos(angle_y);
	double true_z = (n_z*z - ray.get_Z_intercept_X()*n_x - ray.get_Z_intercept_Y()*n_y)/(n_z + ray.get_slope_X()*n_x + ray.get_slope_Y()*n_y);
	if(is_X){
		perp_pos_mm = ray.get_slope_Y()*true_z + ray.get_Z_intercept_Y();
	}
	else{
		perp_pos_mm = ray.get_slope_X()*true_z + ray.get_Z_intercept_X();
	}

}
double Cluster::get_perp_pos_mm() const{
	return perp_pos_mm;
}
int Cluster::find_det(const vector<Detector*> det_array) const{
	int n_pos = -1;
	for(unsigned int i = 0; i<det_array.size();i++){
		if(is_in_det(det_array[i])) n_pos = i;
	}
	return n_pos;
}
bool Cluster::is_in_det(const Detector * det) const{
	return ((type == det->get_type()) && (n_in_tree == det->get_n_in_tree()));
}
Cluster::~Cluster(){

}

CM_Cluster::CM_Cluster(): Cluster(){
	type = Tomography::CM;
}
CM_Cluster::CM_Cluster(const CM_Cluster& other): Cluster(other){
	type = Tomography::CM;
	strip_type = other.strip_type;
}
CM_Cluster& CM_Cluster::operator=(const CM_Cluster& other){
	Cluster::operator=(other);
	type = other.type;
	strip_type = other.strip_type;
	return *this;
}
CM_Cluster::CM_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry): Cluster(treeObject,number_,det,entry){
	if(det->get_type() != Tomography::CM){
		*this = CM_Cluster();
		return;
	}
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	ampl = treeObject->CM_ClusAmpl[n_in_tree][number];
	size = treeObject->CM_ClusSize[n_in_tree][number];
	pos = treeObject->CM_ClusPos[n_in_tree][number];
	maxStripAmpl = treeObject->CM_ClusMaxStripAmpl[n_in_tree][number];
	maxSample = treeObject->CM_ClusMaxSample[n_in_tree][number];
	TOT = treeObject->CM_ClusTOT[n_in_tree][number];
	t = treeObject->CM_ClusT[n_in_tree][number];
	maxStrip = treeObject->CM_ClusMaxStrip[n_in_tree][number];
	type = Tomography::CM;
	(pos>31) ? strip_type = Tomography::Wide : strip_type = Tomography::Thin;
}
CM_Cluster::CM_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_): Cluster(det, number_, pos_, size_, ampl_, maxSample_, maxStripAmpl_, TOT_, t_, maxStrip_){
	if(det->get_type() != Tomography::CM){
		*this = CM_Cluster();
		return;
	}
	type = Tomography::CM;
	(pos>31) ? strip_type = Tomography::Wide : strip_type = Tomography::Thin;
}
CM_Cluster::~CM_Cluster(){
	
}
Cluster * CM_Cluster::Clone() const{
	return new CM_Cluster(*this);
}
Tomography::strip_type CM_Cluster::get_strip_type() const{
	return strip_type;
}
double CM_Cluster::get_pos_mm() const{
	if(strip_type == Tomography::Wide){
		double pos_mm = 0;
		if(direction){
			pos_mm = ((31.-maxStrip)*CM_Detector::wideStripPitch);
		}
		else{
			pos_mm = (maxStrip*CM_Detector::wideStripPitch);
		}
		pos_mm -= CM_Detector::size/2.;
		if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
		else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
		pos_mm += offset;
		return pos_mm;
	}
	else return -1;
}
double CM_Cluster::correct_strip_nb(int strip_nb) const{
	if(strip_type == Tomography::Wide){
		double pos_mm = 0;
		if(direction){
			pos_mm = ((31.-strip_nb)*CM_Detector::wideStripPitch);
		}
		else{
			pos_mm = (strip_nb*CM_Detector::wideStripPitch);
		}
		pos_mm -= CM_Detector::size/2.;
		if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
		else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
		pos_mm += offset;
		return pos_mm;
	}
	else return -1;
}
double CM_Cluster::get_z() const{
	if(strip_type == Tomography::Wide){
		double real_z = z;
		double pos_mm = 0;
		if(direction){
			pos_mm = ((31.-maxStrip)*CM_Detector::wideStripPitch);
		}
		else{
			pos_mm = (maxStrip*CM_Detector::wideStripPitch);
		}
		pos_mm -= CM_Detector::size/2.;
		if(is_X) real_z += (pos_mm)*((Sin(angle_y)/Cos(angle_x)) - Tan(angle_x)*Tan(angle_z)*Cos(angle_y)) + (perp_pos_mm)*(Tan(angle_x)/Cos(angle_z));
		else real_z += (pos_mm*(Cos(angle_z)*Sin(angle_x)*Cos(angle_y) + Sin(angle_z)*Sin(angle_y)) + perp_pos_mm*Cos(angle_x)*Sin(angle_y))/(Cos(angle_z)*Cos(angle_y) + Sin(angle_z)*Sin(angle_x)*Sin(angle_y));
		return real_z;
	}
	else return z;
}

CM_Demux_Cluster::CM_Demux_Cluster(): CM_Cluster(){
	type = Tomography::CM;
	strip_type = Tomography::Demux;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Demux_Cluster& other): CM_Cluster(other){
	type = Tomography::CM;
	strip_type = Tomography::Demux;
}
CM_Demux_Cluster& CM_Demux_Cluster::operator=(const CM_Demux_Cluster& other){
	CM_Cluster::operator=(other);
	type = Tomography::CM;
	strip_type = Tomography::Demux;
	return *this;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Cluster& thinStrip_clus, const CM_Cluster& wideStrip_clus){
	type = Tomography::CM;
	strip_type = Tomography::Demux;
	if(thinStrip_clus.pos>31 || wideStrip_clus.pos<32) return;
	if(thinStrip_clus.evn != wideStrip_clus.evn) return;
	if(thinStrip_clus.n_in_tree != wideStrip_clus.n_in_tree) return;
	n_in_tree = thinStrip_clus.n_in_tree;
	z = thinStrip_clus.z;
	is_X = thinStrip_clus.is_X;
	is_up = thinStrip_clus.is_up;
	number = thinStrip_clus.number;
	offset = thinStrip_clus.offset;
	direction = thinStrip_clus.direction;
	angle_x = thinStrip_clus.angle_x;
	angle_y = thinStrip_clus.angle_y;
	angle_z = thinStrip_clus.angle_z;
	ampl = thinStrip_clus.ampl + wideStrip_clus.ampl;
	size = thinStrip_clus.size*wideStrip_clus.size;
	pos = (31.-thinStrip_clus.pos)+32.*(63.-wideStrip_clus.maxStrip);
	//pos = thinStrip_clus.pos+32.*(63.-wideStrip_clus.maxStrip);
	maxStripAmpl = Min(thinStrip_clus.maxStripAmpl,wideStrip_clus.maxStripAmpl);
	maxSample = (thinStrip_clus.maxSample+wideStrip_clus.maxSample)/2;
	TOT = Min(thinStrip_clus.TOT,wideStrip_clus.TOT);
	t = (thinStrip_clus.t+wideStrip_clus.t)/2;
	maxStrip = (thinStrip_clus.maxStripAmpl>wideStrip_clus.maxStripAmpl) ? thinStrip_clus.maxStrip : wideStrip_clus.maxStrip;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Cluster& wideStrip_clus){
	type = Tomography::CM;
	strip_type = Tomography::Demux;
	if(wideStrip_clus.pos<32) return;
	n_in_tree = wideStrip_clus.n_in_tree;
	z = wideStrip_clus.z;
	is_X = wideStrip_clus.is_X;
	is_up = wideStrip_clus.is_up;
	number = wideStrip_clus.number;
	offset = wideStrip_clus.offset;
	direction = wideStrip_clus.direction;
	angle_x = wideStrip_clus.angle_x;
	angle_y = wideStrip_clus.angle_y;
	angle_z = wideStrip_clus.angle_z;
	ampl = wideStrip_clus.ampl;
	size = 32*wideStrip_clus.size;
	pos = 15.5+32.*(63.-wideStrip_clus.maxStrip);
	maxStripAmpl = wideStrip_clus.maxStripAmpl;
	maxSample = wideStrip_clus.maxSample;
	TOT = wideStrip_clus.TOT;
	t = wideStrip_clus.t;
	maxStrip = wideStrip_clus.maxStrip;
}
double CM_Demux_Cluster::get_pos_mm() const{
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*CM_Detector::thinStripPitch;
	}
	else{
		pos_mm = (1023-pos)*CM_Detector::thinStripPitch;
	}
	pos_mm -= CM_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double CM_Demux_Cluster::correct_strip_nb(int strip_nb) const{
	double pos_mm = 0;
	if(direction){
		pos_mm = strip_nb*CM_Detector::thinStripPitch;
	}
	else{
		pos_mm = (1023-strip_nb)*CM_Detector::thinStripPitch;
	}
	pos_mm -= CM_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double CM_Demux_Cluster::get_z() const{
	double real_z = z;
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*CM_Detector::thinStripPitch;
	}
	else{
		pos_mm = (1023-pos)*CM_Detector::thinStripPitch;
	}
	pos_mm -= CM_Detector::size/2.;
	if(is_X) real_z += (pos_mm)*((Sin(angle_y)/Cos(angle_x)) - Tan(angle_x)*Tan(angle_z)*Cos(angle_y)) + (perp_pos_mm)*(Tan(angle_x)/Cos(angle_z));
	else real_z += (pos_mm*(Cos(angle_z)*Sin(angle_x)*Cos(angle_y) + Sin(angle_z)*Sin(angle_y)) + perp_pos_mm*Cos(angle_x)*Sin(angle_y))/(Cos(angle_z)*Cos(angle_y) + Sin(angle_z)*Sin(angle_x)*Sin(angle_y));
	return real_z;
}
Cluster * CM_Demux_Cluster::Clone() const{
	return new CM_Demux_Cluster(*this);
}
CM_Demux_Cluster::~CM_Demux_Cluster(){
	
}

MG_Cluster::MG_Cluster(): Cluster(){
	type = Tomography::MG;
}
MG_Cluster::MG_Cluster(const MG_Cluster& other): Cluster(other){
	type = Tomography::MG;
}
MG_Cluster& MG_Cluster::operator=(const MG_Cluster& other){
	Cluster::operator=(other);
	type = Tomography::MG;
	return *this;
}
MG_Cluster::MG_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry): Cluster(treeObject,number_,det,entry){
	if(det->get_type() != Tomography::MG){
		*this = MG_Cluster();
		return;
	}
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	ampl = treeObject->MG_ClusAmpl[n_in_tree][number];
	size = treeObject->MG_ClusSize[n_in_tree][number];
	pos = treeObject->MG_ClusPos[n_in_tree][number];
	maxStripAmpl = treeObject->MG_ClusMaxStripAmpl[n_in_tree][number];
	maxSample = treeObject->MG_ClusMaxSample[n_in_tree][number];
	TOT = treeObject->MG_ClusTOT[n_in_tree][number];
	t = treeObject->MG_ClusT[n_in_tree][number];
	maxStrip = treeObject->MG_ClusMaxStrip[n_in_tree][number];
	type = Tomography::MG;
}
MG_Cluster::MG_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_): Cluster(det, number_, pos_, size_, ampl_, maxSample_, maxStripAmpl_, TOT_, t_, maxStrip_){
	if(det->get_type() != Tomography::MG){
		*this = MG_Cluster();
		return;
	}
	type = Tomography::MG;
}
Cluster * MG_Cluster::Clone() const{
	return new MG_Cluster(*this);
}
MG_Cluster::~MG_Cluster(){
	
}
double MG_Cluster::get_pos_mm() const{
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*MG_Detector::StripPitch;
	}
	else{
		pos_mm = (1023-pos)*MG_Detector::StripPitch;
	}
	pos_mm -= MG_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double MG_Cluster::correct_strip_nb(int strip_nb) const{
	double pos_mm = 0;
	if(direction){
		pos_mm = strip_nb*MG_Detector::StripPitch;
	}
	else{
		pos_mm = (1023-strip_nb)*MG_Detector::StripPitch;
	}
	pos_mm -= MG_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double MG_Cluster::get_z() const{
	double real_z = z;
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*MG_Detector::StripPitch;
	}
	else{
		pos_mm = (1023-pos)*MG_Detector::StripPitch;
	}
	pos_mm -= MG_Detector::size/2.;
	if(is_X) real_z += (pos_mm)*((Sin(angle_y)/Cos(angle_x)) - Tan(angle_x)*Tan(angle_z)*Cos(angle_y)) + (perp_pos_mm)*(Tan(angle_x)/Cos(angle_z));
	else real_z += (pos_mm*(Cos(angle_z)*Sin(angle_x)*Cos(angle_y) + Sin(angle_z)*Sin(angle_y)) + perp_pos_mm*Cos(angle_x)*Sin(angle_y))/(Cos(angle_z)*Cos(angle_y) + Sin(angle_z)*Sin(angle_x)*Sin(angle_y));
	return real_z;
}

MGv2_Cluster::MGv2_Cluster(): Cluster(){
	type = Tomography::MGv2;
}
MGv2_Cluster::MGv2_Cluster(const MGv2_Cluster& other): Cluster(other){
	type = Tomography::MGv2;
}
MGv2_Cluster& MGv2_Cluster::operator=(const MGv2_Cluster& other){
	Cluster::operator=(other);
	type = Tomography::MGv2;
	return *this;
}
MGv2_Cluster::MGv2_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry): Cluster(treeObject,number_,det,entry){
	if(det->get_type() != Tomography::MGv2){
		*this = MGv2_Cluster();
		return;
	}
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	ampl = treeObject->MGv2_ClusAmpl[n_in_tree][number];
	size = treeObject->MGv2_ClusSize[n_in_tree][number];
	pos = treeObject->MGv2_ClusPos[n_in_tree][number];
	maxStripAmpl = treeObject->MGv2_ClusMaxStripAmpl[n_in_tree][number];
	maxSample = treeObject->MGv2_ClusMaxSample[n_in_tree][number];
	TOT = treeObject->MGv2_ClusTOT[n_in_tree][number];
	t = treeObject->MGv2_ClusT[n_in_tree][number];
	maxStrip = treeObject->MGv2_ClusMaxStrip[n_in_tree][number];
	type = Tomography::MGv2;
}
MGv2_Cluster::MGv2_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_): Cluster(det, number_, pos_, size_, ampl_, maxSample_, maxStripAmpl_, TOT_, t_, maxStrip_){
	if(det->get_type() != Tomography::MGv2){
		*this = MGv2_Cluster();
		return;
	}
	type = Tomography::MGv2;
}
Cluster * MGv2_Cluster::Clone() const{
	return new MGv2_Cluster(*this);
}
MGv2_Cluster::~MGv2_Cluster(){
	
}
double MGv2_Cluster::get_pos_mm() const{
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*MGv2_Detector::StripPitch;
	}
	else{
		pos_mm = (1036-pos)*MGv2_Detector::StripPitch;
	}
	pos_mm -= MGv2_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double MGv2_Cluster::correct_strip_nb(int strip_nb) const{
	double pos_mm = 0;
	if(direction){
		pos_mm = strip_nb*MGv2_Detector::StripPitch;
	}
	else{
		pos_mm = (1036-strip_nb)*MGv2_Detector::StripPitch;
	}
	pos_mm -= MGv2_Detector::size/2.;
	if(is_X) pos_mm = (pos_mm*Cos(angle_y)/Cos(angle_z)) - perp_pos_mm*Tan(angle_z);
	else pos_mm = (pos_mm*Cos(angle_x) + perp_pos_mm*(Sin(angle_z) - Cos(angle_z)*Sin(angle_x)*Tan(angle_y)))/(Cos(angle_z) + Sin(angle_z)*Sin(angle_x)*Tan(angle_y));
	pos_mm += offset;
	return pos_mm;
}
double MGv2_Cluster::get_z() const{
	double real_z = z;
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*MGv2_Detector::StripPitch;
	}
	else{
		pos_mm = (1036-pos)*MGv2_Detector::StripPitch;
	}
	pos_mm -= MGv2_Detector::size/2.;
	if(is_X) real_z += (pos_mm)*((Sin(angle_y)/Cos(angle_x)) - Tan(angle_x)*Tan(angle_z)*Cos(angle_y)) + (perp_pos_mm)*(Tan(angle_x)/Cos(angle_z));
	else real_z += (pos_mm*(Cos(angle_z)*Sin(angle_x)*Cos(angle_y) + Sin(angle_z)*Sin(angle_y)) + perp_pos_mm*Cos(angle_x)*Sin(angle_y))/(Cos(angle_z)*Cos(angle_y) + Sin(angle_z)*Sin(angle_x)*Sin(angle_y));
	return real_z;
}