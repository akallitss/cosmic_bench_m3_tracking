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

Cluster::Cluster(){
	type.clear();
	evn = -1;
	number = -1;
	ampl = -1;
	size = -1;
	pos = -1;
	maxStripAmpl = -1;
	maxSample = -1;
	TOT = -1;
	t = -1;
	type = "";
	z = -1;
	is_X = -1;
	is_up = -1;
	offset = 0;
	direction = false;
}
Cluster::Cluster(const Cluster& other){
	type.clear();
	type = "";
	type = other.type;
	evn = other.evn;
	number = other.number;
	ampl = other.ampl;
	size = other.size;
	pos = other.pos;
	maxStripAmpl = other.maxStripAmpl;
	maxSample= other.maxSample;
	TOT = other.TOT;
	t = other.t;
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	offset = other.offset;
	direction = other.direction;
}
Cluster& Cluster::operator=(const Cluster& other){
	type.clear();
	type = "";
	type = other.type;
	evn = other.evn;
	number = other.number;
	ampl = other.ampl;
	size = other.size;
	pos = other.pos;
	maxStripAmpl = other.maxStripAmpl;
	maxSample= other.maxSample;
	TOT = other.TOT;
	t = other.t;
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	offset = other.offset;
	direction = other.direction;
	return *this;
}
Cluster::Cluster(T * treeObject, int entry){
	type.clear();
	type = "";
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
	TOT = -1;
	t = -1;
	z = -1;
	is_X = -1;
	is_up = -1;
	offset = 0;
	direction = false;
}
string Cluster::get_type() const{
	return type;
}
bool Cluster::get_is_X() const{
	return is_X;
}
Cluster::~Cluster(){

}

CM_Cluster::CM_Cluster(): Cluster(){
	type.clear();
	type = "CM";
	maxStrip = -1;
	strip_type = "";
	cm_n_in_tree = -1;
}
CM_Cluster::CM_Cluster(const CM_Cluster& other): Cluster(other){
	maxStrip = other.maxStrip;
	type.clear();
	type = "CM";
	strip_type = other.strip_type;
	cm_n_in_tree = other.cm_n_in_tree;
}
CM_Cluster& CM_Cluster::operator=(const CM_Cluster& other){
	Cluster::operator=(other);
	maxStrip = other.maxStrip;
	type.clear();
	type = "CM";
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
	ampl = treeObject->CM_ClusAmpl[cm_n_in_tree][number];
	size = treeObject->CM_ClusSize[cm_n_in_tree][number];
	pos = treeObject->CM_ClusPos[cm_n_in_tree][number];
	maxStripAmpl = treeObject->CM_ClusMaxStripAmpl[cm_n_in_tree][number];
	maxSample = treeObject->CM_ClusMaxSample[cm_n_in_tree][number];
	TOT = treeObject->CM_ClusTOT[cm_n_in_tree][number];
	t = treeObject->CM_ClusT[cm_n_in_tree][number];
	maxStrip = treeObject->CM_ClusMaxStrip[cm_n_in_tree][number];
	type.clear();
	type = "CM";
	(pos>31) ? strip_type = "Wide" : strip_type = "Thin";
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
	if(treeObject->CM_ClusTOT[n_in_tree][number_]<detector->ClusTOTCut_Min){
		//cout << "TOT" << endl;
		return false;
	}
	if(treeObject->CM_ClusMaxSample[n_in_tree][number_]<detector->ClusMaxSampleCut_Min){
		//cout << "ClusMaxSample Min" << endl;
		return false;
	}
	if(treeObject->CM_ClusMaxSample[n_in_tree][number_]>detector->ClusMaxSampleCut_Max){
		//cout << "ClusMaxSample Max" << endl;
		return false;
	}
	if(treeObject->CM_ClusPos[n_in_tree][number_]>31 && treeObject->CM_ClusMaxStripAmpl[n_in_tree][number_]<detector->ClusMaxStripAmplCut_Min_Wide){
		//cout << "ClusMaxStripAmpl Min" << endl;
		return false;
	}
	if(treeObject->CM_ClusPos[n_in_tree][number_]>31 && treeObject->CM_ClusSize[n_in_tree][number_]>detector->ClusSizeCut_Max_Wide){
		//cout << "ClusMaxStripAmpl Max" << endl;
		return false;
	}
	//cout << "gnah" << endl;
	return true;
}
bool CM_Cluster::is_in_det(Detector * det) const{
	if(det->get_type() != "CM") return false;
	return ((dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree() == cm_n_in_tree);
}
string CM_Cluster::get_strip_type() const{
	return strip_type;
}
double CM_Cluster::get_pos_mm() const{
	if(strip_type == "Wide"){
		if(direction){
			return ((63.-maxStrip)*500./32.) + offset;
		}
		else{
			return (maxStrip*500./32.) + offset;
		}
	}
	else return -1;
}

CM_Demux_Cluster::CM_Demux_Cluster(): CM_Cluster(){
	type.clear();
	type = "CM_Demux";
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Demux_Cluster& other): CM_Cluster(other){
	type.clear();
	type = "CM_Demux";
	strip_type = "Demux";
}
CM_Demux_Cluster& CM_Demux_Cluster::operator=(const CM_Demux_Cluster& other){
	CM_Cluster::operator=(other);
	type.clear();
	type = "CM_Demux";
	strip_type = "Demux";
	return *this;
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Cluster& thinStrip_clus, const CM_Cluster& wideStrip_clus){
	type.clear();
	type = "CM_Demux";
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
	ampl = thinStrip_clus.ampl + wideStrip_clus.ampl;
	size = thinStrip_clus.size*wideStrip_clus.size;
	pos = (31.-thinStrip_clus.pos)+32.*(63.-wideStrip_clus.maxStrip);
	//pos = thinStrip_clus.pos+32.*(63.-wideStrip_clus.maxStrip);
	maxStripAmpl = Min(thinStrip_clus.maxStripAmpl,wideStrip_clus.maxStripAmpl);
	maxSample = (thinStrip_clus.maxSample+wideStrip_clus.maxSample)/2;
	TOT = Min(thinStrip_clus.TOT,wideStrip_clus.TOT);
	t = (thinStrip_clus.t+wideStrip_clus.t)/2;
	maxStrip = (thinStrip_clus.maxStripAmpl>wideStrip_clus.maxStripAmpl) ? thinStrip_clus.maxStrip : wideStrip_clus.maxStrip;
	strip_type = "Demux";
}
CM_Demux_Cluster::CM_Demux_Cluster(const CM_Cluster& wideStrip_clus){
	type.clear();
	type = "CM_Demux";
	if(wideStrip_clus.pos<32) return;
	cm_n_in_tree = wideStrip_clus.cm_n_in_tree;
	z = wideStrip_clus.z;
	is_X = wideStrip_clus.is_X;
	is_up = wideStrip_clus.is_up;
	number = wideStrip_clus.number;
	offset = wideStrip_clus.offset;
	direction = wideStrip_clus.direction;
	ampl = wideStrip_clus.ampl;
	size = 32*wideStrip_clus.size;
	pos = 15.5+32.*(63.-wideStrip_clus.maxStrip);
	maxStripAmpl = wideStrip_clus.maxStripAmpl;
	maxSample = wideStrip_clus.maxSample;
	TOT = wideStrip_clus.TOT;
	t = wideStrip_clus.t;
	maxStrip = wideStrip_clus.maxStrip;
	strip_type = "Demux";
}
double CM_Demux_Cluster::get_pos_mm() const{
	if(direction){
		return (pos*500./1024.)+offset;
	}
	else{
		return ((1024-pos)*500./1024.)+offset;
	}
}
CM_Demux_Cluster::~CM_Demux_Cluster(){
	
}

MG_Cluster::MG_Cluster(): Cluster(){
	type.clear();
	type = "MG";
	mg_n_in_tree = -1;
}
MG_Cluster::MG_Cluster(const MG_Cluster& other): Cluster(other){
	type.clear();
	type = "MG";
	mg_n_in_tree = other.mg_n_in_tree;
}
MG_Cluster& MG_Cluster::operator=(const MG_Cluster& other){
	Cluster::operator=(other);
	type.clear();
	type = "MG";
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
	ampl = treeObject->MG_ClusAmpl[mg_n_in_tree][number];
	size = treeObject->MG_ClusSize[mg_n_in_tree][number];
	pos = treeObject->MG_ClusPos[mg_n_in_tree][number];
	maxStripAmpl = treeObject->MG_ClusMaxStripAmpl[mg_n_in_tree][number];
	maxSample = treeObject->MG_ClusMaxSample[mg_n_in_tree][number];
	TOT = treeObject->MG_ClusTOT[mg_n_in_tree][number];
	t = treeObject->MG_ClusT[mg_n_in_tree][number];
	type.clear();
	type = "MG";
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
	if(treeObject->MG_ClusTOT[n_in_tree][number_]<detector->ClusTOTCut_Min) return false;
	if(treeObject->MG_ClusMaxSample[n_in_tree][number_]<detector->ClusMaxSampleCut_Min) return false;
	if(treeObject->MG_ClusMaxSample[n_in_tree][number_]>detector->ClusMaxSampleCut_Max) return false;
	if(treeObject->MG_ClusSize[n_in_tree][number_]<detector->ClusSizeCut_Min) return false;
	return true;
}
bool MG_Cluster::is_in_det(Detector * det) const{
	if(det->get_type() != "MG") return false;
	return ((dynamic_cast<MG_Detector*>(det))->get_mg_n_in_tree() == mg_n_in_tree);
}
double MG_Cluster::get_pos_mm() const{
	if(direction){
		return (pos*500./1024.)+offset;
	}
	else{
		return ((1024-pos)*500./1024.)+offset;
	}
}