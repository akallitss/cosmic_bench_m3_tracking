#define detector_cpp
#include "detector.h"

#include <string>
#include <vector>
#include <iostream>

#include <TMath.h>

using std::cout;
using std::endl;

using TMath::Exp;
using TMath::Log;

bool operator==(Detector const &det1, Detector const &det2){
	if(det1.get_type()!= det2.get_type()) return false;
	else{
		if(det1.get_type() == Tomography::CM){
				if(dynamic_cast<const CM_Detector*>(&det1)->get_cm_n_in_tree() != dynamic_cast<const CM_Detector*>(&det2)->get_cm_n_in_tree()) return false;
				else return true;
		}
		else if(det1.get_type() == Tomography::MG){
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
double Detector::get_angle() const{
	return angle;
}
int Detector::get_perp_n() const{
	return perp_n;
}
int Detector::get_clustering_holes() const{
	return clustering_holes;
}
double Detector::get_RMS(int i) const{
	return RMS[i];
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
bool Detector::test_ClusTOT(double TOT) const{
	if(TOT<ClusTOTCut_Min) return false;
	return true;
}
bool Detector::test_ClusMaxSample(double maxSample) const{
	if(maxSample<ClusMaxSampleCut_Min) return false;
	if(maxSample>ClusMaxSampleCut_Max) return false;
	return true;
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
	angle = 0;
	perp_n = -1;
	clustering_holes = 0;
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
	angle = other.angle;
	RMS = other.RMS;
	perp_n = other.perp_n;
	clustering_holes = other.clustering_holes;
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
	angle = other.angle;
	RMS = other.RMS;
	perp_n = other.perp_n;
	clustering_holes = other.clustering_holes;
	return *this;
}
Detector::Detector(double z_, bool is_X_, bool is_up_, bool is_ref_, double offset_, bool direction_, double angle_, int perp_n_, int clustering_holes_){
	z = z_;
	is_X = is_X_;
	is_up = is_up_;
	is_ref = is_ref_;
	offset = offset_;
	direction = direction_;
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	angle = angle_;
	perp_n = perp_n_;
	clustering_holes = clustering_holes_;
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
CM_Detector::CM_Detector(double z_, bool is_X_, bool is_up_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_, double angle_) :Detector(z_,is_X_,is_up_, is_ref_, offset_, direction_, angle_,-1,0){
	cm_n_in_tree = cm_n;
	use_thin_strip = use_thin_strip_;
	ClusMaxStripAmplCut_Min_Wide = -1;
	ClusSizeCut_Max_Wide = -1;
}
CM_Detector::~CM_Detector(){

}
void CM_Detector::set_ClusMaxStripAmplCut_Min_Wide(double cut){
	ClusMaxStripAmplCut_Min_Wide = cut;
}
void CM_Detector::set_ClusSizeCut_Max_Wide(double cut){
	ClusSizeCut_Max_Wide = cut;
}
bool CM_Detector::test_ClusMaxStripAmpl_Wide(double maxStripAmpl) const{
	if(maxStripAmpl<ClusMaxStripAmplCut_Min_Wide) return false;
	return true;
}
bool CM_Detector::test_ClusSize_Wide(double size) const{
	if(size>ClusSizeCut_Max_Wide) return false;
	return true;
}
int CM_Detector::get_cm_n_in_tree() const{
	return cm_n_in_tree;
}
bool CM_Detector::get_use_thin_strip() const{
	return use_thin_strip;
}
Tomography::det_type CM_Detector::get_type() const{
	return Tomography::CM;
}
void CM_Detector::set_RMS(vector<double> RMS_){
	if(RMS_.size()!=64) return;
	RMS = RMS_;
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
MG_Detector::MG_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_, int perp_n_, int clustering_holes_): Detector(z_,is_X_,is_up_, is_ref_, offset_,direction_, angle_, perp_n_, clustering_holes_){
	mg_n_in_tree = mg_n;
	ClusSizeCut_Min = -1;
}
MG_Detector::~MG_Detector(){

}
static vector<unsigned int> generate_StripToChannel(){
	int p=61; int n=1024;
	int MultiplexSeries[]={30,10,15,19,5,20,27,22,11,24,18,12,9,6,3,4,1,16,8,2,23,21,7,25,14,28,17,26,13,29};
	vector<unsigned int> Detector(n,0); // strip to channel correspondance
	for(int i=0;i<(p-1)/2;i++){
		for(int j=0;j<p;j++){
			if(i*p+j<n){
				Detector[i*p+j]=(0+MultiplexSeries[i]*j)%p;
			}
		}
	}
	return Detector;
}
const vector<unsigned int> MG_Detector::StripToChannel = generate_StripToChannel();
unsigned int MG_Detector::StripToChannel_f(unsigned int strip_nb){
	if(strip_nb>=1024) return -1;
	int p=61; int n=1024;
	int MultiplexSeries[]={30,10,15,19,5,20,27,22,11,24,18,12,9,6,3,4,1,16,8,2,23,21,7,25,14,28,17,26,13,29};
	int Detector[n]; // strip to channel correspondance
	for(int i=0;i<(p-1)/2;i++){
		for(int j=0;j<p;j++){
			if(i*p+j<n){
				Detector[i*p+j]=(0+MultiplexSeries[i]*j)%p;
			}
		}
	}
	return Detector[strip_nb];
}
vector<unsigned int> MG_Detector::ChannelToStrip(unsigned int channel_nb){
	vector<unsigned int> channel_list;
	if(channel_nb>=61) return channel_list;
	int p=61; int n=1024;
	int MultiplexSeries[]={30,10,15,19,5,20,27,22,11,24,18,12,9,6,3,4,1,16,8,2,23,21,7,25,14,28,17,26,13,29};
	unsigned int Detector[n]; // strip to channel correspondance
	for(int i=0;i<(p-1)/2;i++){
		for(int j=0;j<p;j++){
			if(i*p+j<n){
				Detector[i*p+j]=(0+MultiplexSeries[i]*j)%p;
			}
		}
	}
	for(unsigned int i=0;i<61;i++){
		if(Detector[i] == channel_nb) channel_list.push_back(i);
	}
	return channel_list;
}
void MG_Detector::set_ClusSizeCut_Min(double cut){
	ClusSizeCut_Min = cut;
}
bool MG_Detector::test_ClusSize(double size){
	if(size<ClusSizeCut_Min) return false;
	return true;
}
int MG_Detector::get_mg_n_in_tree() const{
	return mg_n_in_tree;
}
Tomography::det_type MG_Detector::get_type() const{
	return Tomography::MG;
}
void MG_Detector::set_RMS(vector<double> RMS_){
	if(RMS_.size()!=61) return;
	RMS = RMS_;
}
void MG_Detector::set_SRF(double offset, double gauss, double lorentz, double ratio){
	srf_offset = offset;
	srf_gauss_width = gauss;
	srf_lorentz_width = lorentz;
	srf_ratio = ratio;	
}
double MG_Detector::SRF_fit(double * x, double * p){
	double position = (x[0]-p[0])/p[1];
	double alpha = srf_lorentz_width/srf_gauss_width;
	double return_value = Exp(-4*Log(2)*(1-srf_ratio)*position*position);
	return_value /= 1 + (4*srf_ratio*position*position/(alpha*alpha));
	return_value *= 1 - srf_offset;
	return_value += srf_offset;
	return return_value;
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