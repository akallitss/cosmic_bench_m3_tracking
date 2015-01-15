#define cluster_cpp
#include "cluster.h"
#include <string>
#include "T.h"
#include "detector.h"
#include "ray.h"
#include <TMath.h>
#include <iostream>

using std::cout;
using std::endl;
using std::string;
using TMath::Min;
using TMath::Tan;

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
	angle = 0;
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
	angle = other.angle;
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
	angle = other.angle;
	return *this;
}
Cluster::Cluster(T * treeObject, int entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
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
	angle = 0;
}
Cluster::Cluster(double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_){
	pos = pos_;
	size = size_;
	ampl = ampl_;
	maxSample = maxSample_;
	maxStripAmpl = maxStripAmpl_;
	maxStrip = maxStrip_;
	TOT = TOT_;
	t = t_;
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
double Cluster::get_z() const{
	return z;
}
void Cluster::set_perp_pos_mm(double coord){
	perp_pos_mm = coord;
}
double Cluster::get_perp_pos_mm() const{
	return perp_pos_mm;
}
int Cluster::find_det(const vector<Detector*> det_array) const{
	int pos = -1;
	for(unsigned int i = 0; i<det_array.size();i++){
		if(is_in_det(det_array[i])) pos = i;
	}
	return pos;
}
Cluster::~Cluster(){

}

CM_Cluster::CM_Cluster(): Cluster(){
	type = Tomography::CM;
	cm_n_in_tree = -1;
}
CM_Cluster::CM_Cluster(const CM_Cluster& other): Cluster(other){
	type = Tomography::CM;
	strip_type = other.strip_type;
	cm_n_in_tree = other.cm_n_in_tree;
}
CM_Cluster& CM_Cluster::operator=(const CM_Cluster& other){
	Cluster::operator=(other);
	type = other.type;
	strip_type = other.strip_type;
	cm_n_in_tree = other.cm_n_in_tree;
	return *this;
}
CM_Cluster::CM_Cluster(T * treeObject,int number_,CM_Detector * det, int entry): Cluster(treeObject,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	number = number_;
	cm_n_in_tree = det->get_cm_n_in_tree();
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	angle = det->get_angle();
	ampl = treeObject->CM_ClusAmpl[cm_n_in_tree][number];
	size = treeObject->CM_ClusSize[cm_n_in_tree][number];
	pos = treeObject->CM_ClusPos[cm_n_in_tree][number];
	maxStripAmpl = treeObject->CM_ClusMaxStripAmpl[cm_n_in_tree][number];
	maxSample = treeObject->CM_ClusMaxSample[cm_n_in_tree][number];
	TOT = treeObject->CM_ClusTOT[cm_n_in_tree][number];
	t = treeObject->CM_ClusT[cm_n_in_tree][number];
	maxStrip = treeObject->CM_ClusMaxStrip[cm_n_in_tree][number];
	type = Tomography::CM;
	(pos>31) ? strip_type = Tomography::Wide : strip_type = Tomography::Thin;
}
CM_Cluster::CM_Cluster(CM_Detector * det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_): Cluster(pos_,size_,ampl_,maxSample_, maxStripAmpl_,TOT_,t_, maxStrip_){
	number = number_;
	cm_n_in_tree = det->get_cm_n_in_tree();
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	angle = det->get_angle();
	type = Tomography::CM;
	(pos>31) ? strip_type = Tomography::Wide : strip_type = Tomography::Thin;
}
CM_Cluster::~CM_Cluster(){
	
}
bool CM_Cluster::is_suitable(T * treeObject,int number_,CM_Detector * detector, int entry){
	//cout << "blah" << endl;
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	int n_in_tree = detector->get_cm_n_in_tree();
	if(treeObject->CM_ClusPos[n_in_tree][number_]>63 || treeObject->CM_ClusPos[n_in_tree][number_]<0) return false;
	if(treeObject->CM_ClusMaxStrip[n_in_tree][number_]>63 || treeObject->CM_ClusMaxStrip[n_in_tree][number_]<0) return false;
	if(!(detector->test_ClusTOT(treeObject->CM_ClusTOT[n_in_tree][number_]))) return false;
	if(!(detector->test_ClusMaxSample(treeObject->CM_ClusMaxSample[n_in_tree][number_]))) return false;
	if(treeObject->CM_ClusPos[n_in_tree][number_]>31){
		if(!(detector->test_ClusMaxStripAmpl_Wide(treeObject->CM_ClusMaxStripAmpl[n_in_tree][number_]))) return false;
		if(!(detector->test_ClusSize_Wide(treeObject->CM_ClusSize[n_in_tree][number_]))) return false;
	}
	return true;
}
bool CM_Cluster::is_suitable(CM_Detector * detector){
	if(!is_in_det(detector)) return false;
	if(pos>1023 || pos<0) return false;
	if(!(detector->test_ClusTOT(TOT))) return false;
	if(!(detector->test_ClusMaxSample(maxSample))) return false;
	if(pos>31){
		if(!(detector->test_ClusMaxStripAmpl_Wide(maxStripAmpl))) return false;
		if(!(detector->test_ClusSize_Wide(size))) return false;
	}
	return true;
}
bool CM_Cluster::is_in_det(Detector * det) const{
	if(det->get_type() != Tomography::CM) return false;
	return ((dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree() == cm_n_in_tree);
}
Tomography::strip_type CM_Cluster::get_strip_type() const{
	return strip_type;
}
double CM_Cluster::get_pos_mm() const{
	if(strip_type == Tomography::Wide){
		double pos_mm = 0;
		if(direction){
			pos_mm = ((63.-maxStrip)*500./32.);
		}
		else{
			pos_mm = (maxStrip*500./32.);
		}
		if(perp_pos_mm>-1){
			pos_mm += (perp_pos_mm - 250)*Tan(angle);
		}
		pos_mm += offset;
		return pos_mm;
	}
	else return -1;
}
double CM_Cluster::correct_strip_nb(int strip_nb) const{
	if(strip_type == Tomography::Wide){
		double pos_mm = 0;
		if(direction){
			pos_mm = ((63.-strip_nb)*500./32.);
		}
		else{
			pos_mm = (strip_nb*500./32.);
		}
		if(perp_pos_mm>-1){
			pos_mm += (perp_pos_mm - 250)*Tan(angle);
		}
		pos_mm += offset;
		return pos_mm;
	}
	else return -1;
}
int CM_Cluster::get_n_in_tree() const{
	return cm_n_in_tree;
}

CM_Demux_Cluster::CM_Demux_Cluster(): CM_Cluster(){
	type = Tomography::CM_Demux;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Demux_Cluster& other): CM_Cluster(other){
	type = Tomography::CM_Demux;
	strip_type = Tomography::Demux;
}
CM_Demux_Cluster& CM_Demux_Cluster::operator=(const CM_Demux_Cluster& other){
	CM_Cluster::operator=(other);
	type = Tomography::CM_Demux;
	strip_type = Tomography::Demux;
	return *this;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Cluster& thinStrip_clus, const CM_Cluster& wideStrip_clus){
	type = Tomography::CM_Demux;
	strip_type = Tomography::Demux;
	if(thinStrip_clus.pos>31 || wideStrip_clus.pos<32) return;
	if(thinStrip_clus.evn != wideStrip_clus.evn) return;
	if(thinStrip_clus.cm_n_in_tree != wideStrip_clus.cm_n_in_tree) return;
	cm_n_in_tree = thinStrip_clus.cm_n_in_tree;
	z = thinStrip_clus.z;
	is_X = thinStrip_clus.is_X;
	is_up = thinStrip_clus.is_up;
	number = thinStrip_clus.number;
	offset = thinStrip_clus.offset;
	direction = thinStrip_clus.direction;
	angle = thinStrip_clus.angle;
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
	type = Tomography::CM_Demux;
	strip_type = Tomography::Demux;
	if(wideStrip_clus.pos<32) return;
	cm_n_in_tree = wideStrip_clus.cm_n_in_tree;
	z = wideStrip_clus.z;
	is_X = wideStrip_clus.is_X;
	is_up = wideStrip_clus.is_up;
	number = wideStrip_clus.number;
	offset = wideStrip_clus.offset;
	direction = wideStrip_clus.direction;
	angle = wideStrip_clus.angle;
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
		pos_mm = pos*500./1024.;
	}
	else{
		pos_mm = (1024-pos)*500./1024.;
	}
	if(perp_pos_mm>-1){
		pos_mm += (perp_pos_mm - 250)*Tan(angle);
	}
	pos_mm += offset;
	return pos_mm;
}
double CM_Demux_Cluster::correct_strip_nb(int strip_nb) const{
	double pos_mm = 0;
	if(direction){
		pos_mm = strip_nb*500./1024.;
	}
	else{
		pos_mm = (1024-strip_nb)*500./1024.;
	}
	if(perp_pos_mm>-1){
		pos_mm += (perp_pos_mm - 250)*Tan(angle);
	}
	pos_mm += offset;
	return pos_mm;
}
CM_Demux_Cluster::~CM_Demux_Cluster(){
	
}

MG_Cluster::MG_Cluster(): Cluster(){
	type = Tomography::MG;
	mg_n_in_tree = -1;
}
MG_Cluster::MG_Cluster(const MG_Cluster& other): Cluster(other){
	type = Tomography::MG;
	mg_n_in_tree = other.mg_n_in_tree;
}
MG_Cluster& MG_Cluster::operator=(const MG_Cluster& other){
	Cluster::operator=(other);
	type = Tomography::MG;
	mg_n_in_tree = other.mg_n_in_tree;
	return *this;
}
MG_Cluster::MG_Cluster(T * treeObject,int number_,MG_Detector * det, int entry): Cluster(treeObject,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	number = number_;
	mg_n_in_tree = det->get_mg_n_in_tree();
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	angle = det->get_angle();
	ampl = treeObject->MG_ClusAmpl[mg_n_in_tree][number];
	size = treeObject->MG_ClusSize[mg_n_in_tree][number];
	pos = treeObject->MG_ClusPos[mg_n_in_tree][number];
	maxStripAmpl = treeObject->MG_ClusMaxStripAmpl[mg_n_in_tree][number];
	maxSample = treeObject->MG_ClusMaxSample[mg_n_in_tree][number];
	TOT = treeObject->MG_ClusTOT[mg_n_in_tree][number];
	t = treeObject->MG_ClusT[mg_n_in_tree][number];
	maxStrip = treeObject->MG_ClusMaxStrip[mg_n_in_tree][number];
	type = Tomography::MG;
}
MG_Cluster::MG_Cluster(MG_Detector * det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_): Cluster(pos_,size_,ampl_,maxSample_, maxStripAmpl_,TOT_,t_,maxStrip_){
	number = number_;
	mg_n_in_tree = det->get_mg_n_in_tree();
	z = det->get_z();
	is_X = det->get_is_X();
	is_up = det->get_is_up();
	offset = det->get_offset();
	direction = det->get_direction();
	angle = det->get_angle();
	type = Tomography::MG;
}

MG_Cluster::~MG_Cluster(){
	
}
bool MG_Cluster::is_suitable(T * treeObject,int number_,MG_Detector * detector, int entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	int n_in_tree = detector->get_mg_n_in_tree();
	if(treeObject->MG_ClusPos[n_in_tree][number_]>1023 || treeObject->MG_ClusPos[n_in_tree][number_]<0) return false;
	if(!(detector->test_ClusTOT(treeObject->MG_ClusTOT[n_in_tree][number_]))) return false;
	if(!(detector->test_ClusMaxSample(treeObject->MG_ClusMaxSample[n_in_tree][number_]))) return false;
	if(!(detector->test_ClusSize(treeObject->MG_ClusSize[n_in_tree][number_]))) return false;
	return true;
}
bool MG_Cluster::is_suitable(MG_Detector * detector){
	if(!is_in_det(detector)) return false;
	if(pos>1023 || pos<0) return false;
	if(!(detector->test_ClusTOT(TOT))) return false;
	if(!(detector->test_ClusMaxSample(maxSample))) return false;
	if(!(detector->test_ClusSize(size))) return false;
	return true;
}
bool MG_Cluster::is_in_det(Detector * det) const{
	if(det->get_type() != Tomography::MG) return false;
	return ((dynamic_cast<MG_Detector*>(det))->get_mg_n_in_tree() == mg_n_in_tree);
}
double MG_Cluster::get_pos_mm() const{
	double pos_mm = 0;
	if(direction){
		pos_mm = pos*500./1024.;
	}
	else{
		pos_mm = (1024-pos)*500./1024.;
	}
	if(perp_pos_mm>-1){
		pos_mm += (perp_pos_mm - 250)*Tan(angle);
	}
	pos_mm += offset;
	return pos_mm;
}
double MG_Cluster::correct_strip_nb(int strip_nb) const{
	double pos_mm = 0;
	if(direction){
		pos_mm = strip_nb*500./1024.;
	}
	else{
		pos_mm = (1024-strip_nb)*500./1024.;
	}
	if(perp_pos_mm>-1){
		pos_mm += (perp_pos_mm - 250)*Tan(angle);
	}
	pos_mm += offset;
	return pos_mm;
}
int MG_Cluster::get_n_in_tree() const{
	return mg_n_in_tree;
}