#define detector_cpp
#include "detector.h"
#include <string>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

bool operator==(Detector const &det1, Detector const &det2){
	if(det1.get_type()!= det2.get_type()) return false;
	else{
		if(det1.get_type() == "CM"){
				if(dynamic_cast<const CM_Detector*>(&det1)->get_cm_n_in_tree() != dynamic_cast<const CM_Detector*>(&det2)->get_cm_n_in_tree()) return false;
				else return true;
		}
		else if(det1.get_type() == "MG"){
				if(dynamic_cast<const MG_Detector*>(&det1)->get_mg_n_in_tree() != dynamic_cast<const MG_Detector*>(&det2)->get_mg_n_in_tree()) return false;
				else return true;
		}
		else return false;
	}
}
bool operator!=(Detector const &det1, Detector const &det2){
	return !(det1==det2);
}

double Detector::get_z() const{
	return z;
}
bool Detector::get_is_X() const{
	return is_X;
}
bool Detector::get_is_up() const{
	return is_up;
}
bool Detector::get_is_ref() const{
	return is_ref;
}
double Detector::get_offset() const{
	return offset;
}
bool Detector::get_direction() const{
	return direction;
}
void Detector::set_ClusTOTCut_Min(double cut){
	ClusTOTCut_Min = cut;
}
void Detector::set_ClusMaxSampleCut_Min(double cut){
	ClusMaxSampleCut_Min = cut;
}
void Detector::set_ClusMaxSampleCut_Max(double cut){
	ClusMaxSampleCut_Max = cut;
}
Detector::Detector(){
	z = -1;
	is_X = false;
	is_up = false;
	is_ref = false;
	offset = 0;
	direction = false;
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
}
Detector::Detector(const Detector& other){
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	is_ref = other.is_ref;
	offset = other.offset;
	direction = other.direction;
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
}
Detector& Detector::operator=(const Detector& other){
	z = other.z;
	is_X = other.is_X;
	is_up = other.is_up;
	is_ref = other.is_ref;
	offset = other.offset;
	direction = other.direction;
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	return *this;
}
Detector::Detector(double z_, bool is_X_, bool is_up_, bool is_ref_, double offset_, bool direction_){
	z = z_;
	is_X = is_X_;
	is_up = is_up_;
	is_ref = is_ref_;
	offset = offset_;
	direction = direction_;
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
}
Detector::~Detector(){

}

CM_Detector::CM_Detector() : Detector(){
	cm_n_in_tree = -1;
	use_thin_strip = false;
	ClusMaxStripAmplCut_Min_Wide = -1;
	ClusSizeCut_Max_Wide = -1;
}
CM_Detector::CM_Detector(const CM_Detector& other) : Detector(other){
	cm_n_in_tree = other.cm_n_in_tree;
	use_thin_strip = other.use_thin_strip;
	ClusMaxStripAmplCut_Min_Wide = other.ClusMaxStripAmplCut_Min_Wide;
	ClusSizeCut_Max_Wide = other.ClusSizeCut_Max_Wide;
}
CM_Detector& CM_Detector::operator=(const CM_Detector& other){
	Detector::operator=(other);
	cm_n_in_tree = other.cm_n_in_tree;
	use_thin_strip = other.use_thin_strip;
	ClusMaxStripAmplCut_Min_Wide = other.ClusMaxStripAmplCut_Min_Wide;
	ClusSizeCut_Max_Wide = other.ClusSizeCut_Max_Wide;
	return *this;
}
CM_Detector::CM_Detector(double z_, bool is_X_, bool is_up_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_) :Detector(z_,is_X_,is_up_, is_ref_, offset_, direction_){
	cm_n_in_tree = cm_n;
	use_thin_strip = use_thin_strip_;
	ClusMaxStripAmplCut_Min_Wide = -1;
	ClusSizeCut_Max_Wide = -1;
}
CM_Detector::~CM_Detector(){

}
const double CM_Detector::thinStripPitch = 500./1024.;
const double CM_Detector::wideStripPitch = 500./32.;
void CM_Detector::set_ClusMaxStripAmplCut_Min_Wide(double cut){
	ClusMaxStripAmplCut_Min_Wide = cut;
}
void CM_Detector::set_ClusSizeCut_Max_Wide(double cut){
	ClusSizeCut_Max_Wide = cut;
}
int CM_Detector::get_cm_n_in_tree() const{
	return cm_n_in_tree;
}
bool CM_Detector::get_use_thin_strip() const{
	return use_thin_strip;
}
string CM_Detector::get_type() const{
	return "CM";
}

MG_Detector::MG_Detector(): Detector(){
	mg_n_in_tree = -1;
	ClusSizeCut_Min = -1;
}
MG_Detector::MG_Detector(const MG_Detector& other): Detector(other){
	mg_n_in_tree = other.mg_n_in_tree;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
}
MG_Detector& MG_Detector::operator=(const MG_Detector& other){
	Detector::operator=(other);
	mg_n_in_tree = other.mg_n_in_tree;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
	return *this;
}
MG_Detector::MG_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_): Detector(z_,is_X_,is_up_, is_ref_, offset_,direction_){
	mg_n_in_tree = mg_n;
	ClusSizeCut_Min = -1;
}
MG_Detector::~MG_Detector(){

}
const double MG_Detector::StripPitch = 500./1024.; // distance between the middle of 2 adjacent strips
void MG_Detector::set_ClusSizeCut_Min(double cut){
	ClusSizeCut_Min = cut;
}
int MG_Detector::get_mg_n_in_tree() const{
	return mg_n_in_tree;
}
string MG_Detector::get_type() const{
	return "MG";
}

CosmicBench::CosmicBench(){
	detectors.clear();
	CM_N = 0;
	MG_N = 0;
}
CosmicBench::~CosmicBench(){
	for(unsigned int i=0;i<detectors.size();i++){
		delete detectors[i];
	}
	detectors.clear();
	CM_N = 0;
	MG_N = 0;
}
/*
void CosmicBench::add_MM(Detector * det){
	if(det->get_type() == "MG") MG_N++;
	else if(det->get_type() == "CM") CM_N++;
	detectors.push_back(det);
}
*/
int CosmicBench::get_CM_N() const{
	return CM_N;
}
int CosmicBench::get_MG_N() const{
	return MG_N;
}
Detector * CosmicBench::get_detector(unsigned int i) const{
	return detectors[i];
}