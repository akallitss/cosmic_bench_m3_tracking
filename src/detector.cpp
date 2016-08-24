#define detector_cpp
#include "detector.h"
#include "Tanalyse_R.h"
#include "event.h"
#include "cluster.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <limits>

//Boost
#include <boost/foreach.hpp>

#include <TMath.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::pair;
using std::numeric_limits;

using TMath::Exp;
using TMath::Log;
/*
static map<const Tomography::det_type,const Detector* const> Static_Detector_build(){
	map<const Tomography::det_type,const Detector* const> return_map;
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::CM,new CM_Detector()));
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::MG,new MG_Detector()));
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::MGv2,new MGv2_Detector()));
	return return_map;
}

map<const Tomography::det_type,const Detector* const> Tomography::Static_Detector = Static_Detector_build();
*/
bool operator==(Detector const &det1, Detector const &det2){
	return ((det1.get_type() == det2.get_type()) && det1.get_n_in_tree() == det1.get_n_in_tree());
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
int Detector::get_layer() const{
	return layer;
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
double Detector::get_angle_x() const{
	return angle_x;
}
double Detector::get_angle_y() const{
	return angle_y;
}
double Detector::get_angle_z() const{
	return angle_z;
}
int Detector::get_perp_n() const{
	return perp_n;
}
Tomography::det_type Detector::get_perp_type() const{
	return this->get_type();
}
int Detector::get_clustering_holes() const{
	return clustering_holes;
}
double Detector::get_RMS(int i) const{
	return RMS[i];
}
int Detector::get_n_in_tree() const{
	return n_in_tree;
}
vector<pair<int,bool> > Detector::get_asic_n() const{
	return asic_n;
}
Detector::Detector(){
	z = -1;
	is_X = false;
	layer = -1;
	is_ref = false;
	offset = 0;
	direction = false;
	angle_x = 0;
	angle_y = 0;
	angle_z = 0;
	perp_n = -1;
	clustering_holes = 0;
	n_in_tree = -1;
	asic_n.clear();
}
Detector::Detector(const Detector& other){
	z = other.z;
	is_X = other.is_X;
	layer = other.layer;
	is_ref = other.is_ref;
	offset = other.offset;
	direction = other.direction;
	angle_x = other.angle_x;
	angle_y = other.angle_y;
	angle_z = other.angle_z;
	RMS = other.RMS;
	perp_n = other.perp_n;
	clustering_holes = other.clustering_holes;
	n_in_tree = other.n_in_tree;
	asic_n = other.asic_n;
}
Detector& Detector::operator=(const Detector& other){
	z = other.z;
	is_X = other.is_X;
	layer = other.layer;
	is_ref = other.is_ref;
	offset = other.offset;
	direction = other.direction;
	angle_x = other.angle_x;
	angle_y = other.angle_y;
	angle_z = other.angle_z;
	RMS = other.RMS;
	perp_n = other.perp_n;
	clustering_holes = other.clustering_holes;
	n_in_tree = other.n_in_tree;
	asic_n = other.asic_n;
	return *this;
}
Detector::Detector(double z_, bool is_X_, int layer_, int det_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, vector<pair<int,bool> > asic_n_){
	z = z_;
	is_X = is_X_;
	layer = layer_;
	is_ref = is_ref_;
	offset = offset_;
	direction = direction_;
	angle_x = angle_x_;
	angle_y = angle_y_;
	angle_z = angle_z_;
	perp_n = perp_n_;
	clustering_holes = clustering_holes_;
	n_in_tree = det_n;
	asic_n = asic_n_;
}
Detector::~Detector(){

}

CM_Detector::CM_Detector() : Detector(){
	use_thin_strip = false;
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusMaxStripAmplCut_Min_Wide = -1;
	ClusSizeCut_Max_Wide = -1;
}
CM_Detector::CM_Detector(const CM_Detector& other) : Detector(other){
	use_thin_strip = other.use_thin_strip;
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusMaxStripAmplCut_Min_Wide = other.ClusMaxStripAmplCut_Min_Wide;
	ClusSizeCut_Max_Wide = other.ClusSizeCut_Max_Wide;
}
CM_Detector& CM_Detector::operator=(const CM_Detector& other){
	Detector::operator=(other);
	use_thin_strip = other.use_thin_strip;
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusMaxStripAmplCut_Min_Wide = other.ClusMaxStripAmplCut_Min_Wide;
	ClusSizeCut_Max_Wide = other.ClusSizeCut_Max_Wide;
	return *this;
}
CM_Detector::CM_Detector(double z_, bool is_X_, int layer_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int asic_n_, bool connector_direction) :Detector(z_,is_X_,layer_,cm_n, is_ref_, offset_, direction_, angle_x_, angle_y_, angle_z_,-1,0,vector<pair<int,bool> >(1,pair<int,bool>(asic_n_,connector_direction))){
	use_thin_strip = use_thin_strip_;
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusMaxStripAmplCut_Min_Wide = -1;
	ClusSizeCut_Max_Wide = -1;
}
Detector * CM_Detector::Clone() const{
	return new CM_Detector(*this);
}
Event * CM_Detector::build_event(Tanalyse_R * treeObject, int entry) const{
	return new CM_Event(treeObject, this, entry);
}
Event * CM_Detector::build_event(const Tanalyse_R * const treeObject) const{
	return new CM_Event(treeObject, this);
}
Event * CM_Detector::build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const{
	return new CM_Event(this, strip_ampl_,evn_,evttime_);
}
CM_Detector::~CM_Detector(){

}
Detector * CM_Detector::build_det(const ptree::value_type& child) const{
	CM_Detector * current_det  = new CM_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("layer"),child.second.get<int>("cm_n"),child.second.get<bool>("use_thin_strip"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle_x"),child.second.get<double>("angle_y"),child.second.get<double>("angle_z"),child.second.get<int>("asic_n"), child.second.get<bool>("connector_direction"));
	current_det->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
	current_det->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
	current_det->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
	current_det->set_ClusMaxStripAmplCut_Min_Wide(child.second.get<double>("ClusMaxStripAmplCut_Min_Wide"));
	current_det->set_ClusSizeCut_Max_Wide(child.second.get<double>("ClusSizeCut_Max_Wide"));
	return current_det;
}
void CM_Detector::set_ClusMaxStripAmplCut_Min_Wide(double cut){
	ClusMaxStripAmplCut_Min_Wide = cut;
}
void CM_Detector::set_ClusSizeCut_Max_Wide(double cut){
	ClusSizeCut_Max_Wide = cut;
}
void CM_Detector::set_ClusTOTCut_Min(double cut){
	ClusTOTCut_Min = cut;
}
void CM_Detector::set_ClusMaxSampleCut_Min(double cut){
	ClusMaxSampleCut_Min = cut;
}
void CM_Detector::set_ClusMaxSampleCut_Max(double cut){
	ClusMaxSampleCut_Max = cut;
}
bool CM_Detector::get_use_thin_strip() const{
	return use_thin_strip;
}
Tomography::det_type CM_Detector::get_type() const{
	return Tomography::CM;
}
void CM_Detector::set_RMS(vector<double> RMS_){
	if(RMS_.size()!=Nchannel) return;
	RMS = RMS_;
}
double CM_Detector::get_size() const{
	return size;
}
int CM_Detector::get_Nchannel() const{
	return Nchannel;
}
int CM_Detector::get_Nstrip() const{
	return Nstrip;
}
double CM_Detector::get_StripPitch() const{
	return -1;
}
int CM_Detector::get_CMN_div() const{
	return CMN_div;
}
unsigned int CM_Detector::StripToChannel(unsigned int i) const{
	return i;
}
bool CM_Detector::is_suitable(const Cluster * const clus) const{
	if(!(clus->is_in_det(this))) return false;
	if(clus->get_pos() < 0) return false;
	if(clus->get_pos() > 63) return false;
	if(clus->get_TOT() < ClusTOTCut_Min) return false;
	if(clus->get_maxSample() < ClusMaxSampleCut_Min) return false;
	if(clus->get_maxSample() > ClusMaxSampleCut_Max) return false;
	if(dynamic_cast<const CM_Cluster * const>(clus)->get_strip_type() == Tomography::Wide){
		if(clus->get_size() > ClusSizeCut_Max_Wide) return false;
		if(clus->get_maxStripAmpl() < ClusMaxStripAmplCut_Min_Wide) return false;
	}
	return true;
}
int CM_Detector::feminos_mapping(int channel, bool connector_direction) const{
	if(connector_direction) return channel;
	else return Tomography::Nchannel - 1 - channel;
}
int CM_Detector::dream_mapping(int channel, bool connector_direction) const{
	if(connector_direction) return channel;
	else return Tomography::Nchannel - 1 - channel;
}
string CM_Detector::Name() const{
	return "CM";
}
int CM_Detector::get_MaxNClus() const{
	return MaxNClus;
}

MG_Detector::MG_Detector(): Detector(){
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusSizeCut_Min = -1;
}
MG_Detector::MG_Detector(const MG_Detector& other): Detector(other){
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
}
MG_Detector& MG_Detector::operator=(const MG_Detector& other){
	Detector::operator=(other);
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
	return *this;
}
MG_Detector::MG_Detector(double z_, bool is_X_, int layer_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_, bool connector_direction): Detector(z_,is_X_,layer_,mg_n, is_ref_, offset_,direction_, angle_x_, angle_y_, angle_z_, perp_n_, clustering_holes_,vector<pair<int,bool> >(1,pair<int,bool>(asic_n_,connector_direction))){
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusSizeCut_Min = -1;
}
Detector * MG_Detector::Clone() const{
	return new MG_Detector(*this);
}
Event * MG_Detector::build_event(Tanalyse_R * treeObject, int entry) const{
	return new MG_Event(treeObject, this, entry);
}
Event * MG_Detector::build_event(const Tanalyse_R * const treeObject) const{
	return new MG_Event(treeObject, this);
}
Event * MG_Detector::build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const{
	return new MG_Event(this, strip_ampl_,evn_,evttime_);
}
MG_Detector::~MG_Detector(){

}
Detector * MG_Detector::build_det(const ptree::value_type& child) const{
	MG_Detector * current_det = new MG_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("layer"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle_x"),child.second.get<double>("angle_y"),child.second.get<double>("angle_z"),child.second.get<int>("2D_perp_n"),child.second.get<int>("clustering_holes"),child.second.get<int>("asic_n"), child.second.get<bool>("connector_direction"));
	current_det->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
	current_det->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
	current_det->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
	current_det->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
	return current_det;
}
static vector<unsigned int> generate_StripToChannel_MG(){
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
const vector<unsigned int> MG_Detector::StripToChannel_a = generate_StripToChannel_MG();
unsigned int MG_Detector::StripToChannel(unsigned int i) const{
	return StripToChannel_a[i];
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
void MG_Detector::set_ClusTOTCut_Min(double cut){
	ClusTOTCut_Min = cut;
}
void MG_Detector::set_ClusMaxSampleCut_Min(double cut){
	ClusMaxSampleCut_Min = cut;
}
void MG_Detector::set_ClusMaxSampleCut_Max(double cut){
	ClusMaxSampleCut_Max = cut;
}
double MG_Detector::get_ClusSizeCut_Min() const{
	return ClusSizeCut_Min;
}
Tomography::det_type MG_Detector::get_type() const{
	return Tomography::MG;
}
void MG_Detector::set_RMS(vector<double> RMS_){
	if(RMS_.size()!=Nchannel) return;
	RMS = RMS_;
}
void MG_Detector::set_SRF(double offset_, double gauss, double lorentz, double ratio){
	srf_offset = offset_;
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
double MG_Detector::get_size() const{
	return size;
}
int MG_Detector::get_Nchannel() const{
	return Nchannel;
}
int MG_Detector::get_Nstrip() const{
	return Nstrip;
}
double MG_Detector::get_StripPitch() const{
	return StripPitch;
}
int MG_Detector::get_CMN_div() const{
	return CMN_div;
}
bool MG_Detector::is_suitable(const Cluster * const clus) const{
	if(!(clus->is_in_det(this))) return false;
	if(clus->get_pos() < 0) return false;
	if(clus->get_pos() > 1023) return false;
	if(clus->get_TOT() < ClusTOTCut_Min) return false;
	if(clus->get_maxSample() < ClusMaxSampleCut_Min) return false;
	if(clus->get_maxSample() > ClusMaxSampleCut_Max) return false;
	if(clus->get_size() < ClusSizeCut_Min) return false;
	return true;
}
int MG_Detector::feminos_mapping(int channel, bool connector_direction) const{
	int tmpchan = (connector_direction) ? channel : (Tomography::Nchannel - 1 - channel);
	tmpchan = tmpchan - 2 - (tmpchan>13) - (tmpchan>24) - (tmpchan>47) - (tmpchan>58);
	if(tmpchan>15 && tmpchan<48) return tmpchan;
	else return (tmpchan + 1 - (2*(tmpchan%2)));
}
int MG_Detector::dream_mapping(int channel, bool connector_direction) const{
	int tmpchan = (connector_direction) ? channel : (Tomography::Nchannel - 1 - channel);
	return (tmpchan + 1 - (2*(tmpchan%2)));
}
string MG_Detector::Name() const{
	return "MG";
}
int MG_Detector::get_MaxNClus() const{
	return MaxNClus;
}

MGv2_Detector::MGv2_Detector(): Detector(){
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusSizeCut_Min = -1;
}
MGv2_Detector::MGv2_Detector(const MGv2_Detector& other): Detector(other){
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
}
MGv2_Detector& MGv2_Detector::operator=(const MGv2_Detector& other){
	Detector::operator=(other);
	ClusTOTCut_Min = other.ClusTOTCut_Min;
	ClusMaxSampleCut_Min = other.ClusMaxSampleCut_Min;
	ClusMaxSampleCut_Max = other.ClusMaxSampleCut_Max;
	ClusSizeCut_Min = other.ClusSizeCut_Min;
	return *this;
}
MGv2_Detector::MGv2_Detector(double z_, bool is_X_, int layer_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_, bool connector_direction): Detector(z_,is_X_,layer_,mg_n, is_ref_, offset_,direction_, angle_x_, angle_y_, angle_z_, perp_n_, clustering_holes_,vector<pair<int,bool> >(1,pair<int,bool>(asic_n_,connector_direction))){
	ClusTOTCut_Min = -1;
	ClusMaxSampleCut_Min = -1;
	ClusMaxSampleCut_Max = -1;
	ClusSizeCut_Min = -1;
}
Detector * MGv2_Detector::Clone() const{
	return new MGv2_Detector(*this);
}
Event * MGv2_Detector::build_event(Tanalyse_R * treeObject, int entry) const{
	return new MGv2_Event(treeObject, this, entry);
}
Event * MGv2_Detector::build_event(const Tanalyse_R * const treeObject) const{
	return new MGv2_Event(treeObject, this);
}
Event * MGv2_Detector::build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const{
	return new MGv2_Event(this, strip_ampl_,evn_,evttime_);
}
MGv2_Detector::~MGv2_Detector(){

}
Detector * MGv2_Detector::build_det(const ptree::value_type& child) const{
	MGv2_Detector * current_det = new MGv2_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("layer"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle_x"),child.second.get<double>("angle_y"),child.second.get<double>("angle_z"),child.second.get<int>("2D_perp_n"),child.second.get<int>("clustering_holes"),child.second.get<int>("asic_n"),child.second.get<bool>("connector_direction"));
	current_det->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
	current_det->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
	current_det->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
	current_det->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
	return current_det;
}
static vector<unsigned int> generate_StripToChannel_MGv2(){
	int p=61; int n=1037;
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
const vector<unsigned int> MGv2_Detector::StripToChannel_a = generate_StripToChannel_MGv2();
unsigned int MGv2_Detector::StripToChannel(unsigned int i) const{
	return StripToChannel_a[i];
}

vector<unsigned int> MGv2_Detector::ChannelToStrip(unsigned int channel_nb){
	vector<unsigned int> channel_list;
	if(channel_nb>=61) return channel_list;
	int p=61; int n=1037;
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
void MGv2_Detector::set_ClusSizeCut_Min(double cut){
	ClusSizeCut_Min = cut;
}
void MGv2_Detector::set_ClusTOTCut_Min(double cut){
	ClusTOTCut_Min = cut;
}
void MGv2_Detector::set_ClusMaxSampleCut_Min(double cut){
	ClusMaxSampleCut_Min = cut;
}
void MGv2_Detector::set_ClusMaxSampleCut_Max(double cut){
	ClusMaxSampleCut_Max = cut;
}
double MGv2_Detector::get_ClusSizeCut_Min() const{
	return ClusSizeCut_Min;
}
Tomography::det_type MGv2_Detector::get_type() const{
	return Tomography::MGv2;
}
void MGv2_Detector::set_RMS(vector<double> RMS_){
	if(RMS_.size()!=Nchannel) return;
	RMS = RMS_;
}
void MGv2_Detector::set_SRF(double offset_, double gauss, double lorentz, double ratio){
	srf_offset = offset_;
	srf_gauss_width = gauss;
	srf_lorentz_width = lorentz;
	srf_ratio = ratio;	
}
double MGv2_Detector::SRF_fit(double * x, double * p){
	double position = (x[0]-p[0])/p[1];
	double alpha = srf_lorentz_width/srf_gauss_width;
	double return_value = Exp(-4*Log(2)*(1-srf_ratio)*position*position);
	return_value /= 1 + (4*srf_ratio*position*position/(alpha*alpha));
	return_value *= 1 - srf_offset;
	return_value += srf_offset;
	return return_value;
}
double MGv2_Detector::get_size() const{
	return size;
}
int MGv2_Detector::get_Nchannel() const{
	return Nchannel;
}
int MGv2_Detector::get_Nstrip() const{
	return Nstrip;
}
double MGv2_Detector::get_StripPitch() const{
	return StripPitch;
}
int MGv2_Detector::get_CMN_div() const{
	return CMN_div;
}
bool MGv2_Detector::is_suitable(const Cluster * const clus) const{
	if(!(clus->is_in_det(this))) return false;
	if(clus->get_pos() < 0) return false;
	if(clus->get_pos() > 1036) return false;
	if(clus->get_TOT() < ClusTOTCut_Min) return false;
	if(clus->get_maxSample() < ClusMaxSampleCut_Min) return false;
	if(clus->get_maxSample() > ClusMaxSampleCut_Max) return false;
	if(clus->get_size() < ClusSizeCut_Min) return false;
	return true;
}
int MGv2_Detector::feminos_mapping(int channel, bool connector_direction) const{
	// Warning !!!
	//This part was never tested
	int tmpchan = (connector_direction) ? channel : (Tomography::Nchannel - 1 - channel);
	tmpchan = tmpchan - 2 - (tmpchan>13) - (tmpchan>24) - (tmpchan>47) - (tmpchan>58);
	if(tmpchan>15 && tmpchan<48) return tmpchan;
	else return (tmpchan + 1 - (2*(tmpchan%2)));
}
int MGv2_Detector::dream_mapping(int channel, bool connector_direction) const{
	int tmpchan = (connector_direction) ? channel : (Tomography::Nchannel - 1 - channel);
	if(tmpchan == 6) return 8;
	else if(tmpchan == 8) return 6;
	else return tmpchan;
	/*
	if(channel == 57) return 8;
	else if(channel == 55) return 6;
	else return (63 - channel);
	*/
}
string MGv2_Detector::Name() const{
	return "MGv2";
}
int MGv2_Detector::get_MaxNClus() const{
	return MaxNClus;
}


CosmicBench::CosmicBench(){
	for(unsigned int i=0;i<detectors.size();i++){
		delete detectors[i];
	}
	detectors.clear();
	det_n.clear();
	non_ref_n = 0;
}
CosmicBench::CosmicBench(const CosmicBench& other){
	for(unsigned int i=0;i<detectors.size();i++){
		delete detectors[i];
	}
	detectors.clear();
	det_n = other.det_n;
	for(vector<Detector*>::const_iterator it = other.detectors.begin();it!=other.detectors.end();++it){
		detectors.push_back((*it)->Clone());
	}
	non_ref_n = other.non_ref_n;
}
CosmicBench& CosmicBench::operator=(const CosmicBench& other){
	for(unsigned int i=0;i<detectors.size();i++){
		delete detectors[i];
	}
	detectors.clear();
	det_n = other.det_n;
	for(vector<Detector*>::const_iterator it = other.detectors.begin();it!=other.detectors.end();++it){
		detectors.push_back((*it)->Clone());
	}
	non_ref_n = other.non_ref_n;
	return *this;
}
CosmicBench::~CosmicBench(){
	for(unsigned int i=0;i<detectors.size();i++){
		delete detectors[i];
	}
	detectors.clear();
	det_n.clear();
}
CosmicBench::CosmicBench(ptree config_tree){
	Init(config_tree);
}
map<Tomography::det_type,vector<vector<double> > > CosmicBench::read_pedfile(string filename, map<Tomography::det_type,unsigned short> det_n_){
	ifstream in;
	in.open(filename.c_str());
	int current_det = 0;
	int current_strip = 0;
	int last_det = -1;
	double current_rms = 0;
	map<Tomography::det_type,vector<vector<double> > > RMS;
	for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_n_.begin();type_it!=det_n_.end();++type_it){
		if(type_it->second == 0) continue;
		RMS[type_it->first] = vector<vector<double> >(type_it->second,vector<double>(Tomography::Static_Detector[type_it->first]->get_Nchannel()));
	}
	map<Tomography::det_type,unsigned short>::iterator det_it = det_n_.begin();
	while(det_it->second == 0 && det_it!=det_n_.end()) ++det_it;
	while(in.good() && det_it!=det_n_.end()){
		in >> current_det >> current_strip >> current_rms;
		if(current_det<0 || current_strip<0){
			cout << "problem reading pedfile" << endl;
		}
		if(current_det < last_det && current_strip == 0){
			++det_it;
			while(det_it->second == 0 && det_it!=det_n_.end()) ++det_it;
			if(det_it==det_n_.end()) break;
		}
		if(current_det >= det_it->second){
			cout << "too much " << det_it->first << " detectors" << endl;
		}
		if(current_strip > Tomography::Static_Detector[det_it->first]->get_Nchannel()){
			cout << "too much strip for " << det_it->first << "_" << current_det << endl;
		}
		RMS[det_it->first][current_det][current_strip] = current_rms;
		last_det = current_det;
	}
	in.close();
	return RMS;
}
void CosmicBench::Init(ptree config_tree){
	//int total_CM_N = config_tree.get<int>("total_CM_N");
	//int total_MG_N = config_tree.get<int>("total_MG_N");
	string RMSName = config_tree.get<string>("RMSPed");
	/*
	ifstream in;
	in.open(RMSName.c_str());
	int rms_strip, det;
	vector<vector<double> > RMS;
	for(int i=0;i<total_CM_N;i++){
		RMS.push_back(vector<double>(CM_Detector::Nchannel,0));
	}
	for(int i=0;i<total_MG_N;i++){
		RMS.push_back(vector<double>(MG_Detector::Nchannel,0));
	}
	int n_lines = 0;
	while(in.good() && n_lines<((total_CM_N*CM_Detector::Nchannel)+(total_MG_N*MG_Detector::Nchannel))){
		double current_rms;
		in >> det >> rms_strip >> current_rms;
		if(det<0 || det>(total_MG_N+total_CM_N-1)){
			cout << "problem reading RMS file" << endl;
		}
		else if(det<total_CM_N && rms_strip>63){
			cout << "problem reading RMS file" << endl;
		}
		else if(rms_strip>60){
			cout << "problem reading RMS file" << endl;
		}
		RMS[det][rms_strip] = current_rms;
		n_lines++;
	}
	in.close();
	*/
	detectors.clear();
	non_ref_n = 0;
	for(map<const Tomography::det_type,const Detector* const>::iterator type_it=Tomography::Static_Detector.begin();type_it!=Tomography::Static_Detector.end();++type_it){
		ostringstream childname;
		childname << "CosmicBench." << type_it->first;
		if(config_tree.get_child_optional(childname.str())){
			BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child(childname.str())){
				detectors.push_back(type_it->second->build_det(child));
				if(!((detectors.back())->get_is_ref())) non_ref_n++;
				det_n[type_it->first]++;
			}
		}
	}
	map<Tomography::det_type,vector<vector<double> > > RMS = read_pedfile(RMSName,det_n);
	for(vector<Detector*>::iterator det_it=detectors.begin();det_it!=detectors.end();++det_it){
		(*det_it)->set_RMS(RMS[(*det_it)->get_type()][(*det_it)->get_n_in_tree()]);
	}
	/*
	detectors.clear();
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		detectors.push_back(new CM_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("layer"),child.second.get<int>("cm_n"),child.second.get<bool>("use_thin_strip"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle_x"),child.second.get<double>("angle_y"),child.second.get<double>("angle_z"),child.second.get<double>("asic_n")));
		CM_Detector * current_det = dynamic_cast<CM_Detector*>(detectors.back());
		current_det->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		current_det->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		current_det->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		current_det->set_ClusMaxStripAmplCut_Min_Wide(child.second.get<double>("ClusMaxStripAmplCut_Min_Wide"));
		current_det->set_ClusSizeCut_Max_Wide(child.second.get<double>("ClusSizeCut_Max_Wide"));
		detectors.back()->set_RMS(RMS[child.second.get<int>("cm_n")]);
		det_n[Tomography::CM]++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		detectors.push_back(new MG_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("layer"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle_x"),child.second.get<double>("angle_y"),child.second.get<double>("angle_z"),child.second.get<int>("2D_perp_n"),child.second.get<int>("clustering_holes"),child.second.get<double>("angle_z")));
		MG_Detector * current_det = dynamic_cast<MG_Detector*>(detectors.back());
		current_det->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		current_det->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		current_det->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		current_det->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
		detectors.back()->set_RMS(RMS[total_CM_N+child.second.get<int>("mg_n")]);
		det_n[Tomography::MG]++;
	}
	if((total_CM_N!=det_n[Tomography::CM]) || (total_MG_N!=det_n[Tomography::MG])){
		cout << "problem in detectors number" << endl;
		return;
	}
	*/
}
int CosmicBench::get_det_N(Tomography::det_type det_t) const{
	return (det_n.find(det_t))->second;
}
map<Tomography::det_type,unsigned short> CosmicBench::get_det_N() const{
	return det_n;
}
int CosmicBench::get_det_N_tot() const{
	int tot = 0;
	for(map<Tomography::det_type,unsigned short>::const_iterator it=det_n.begin();it!=det_n.end();++it){
		tot += it->second;
	}
	return tot;
}
Detector * CosmicBench::get_detector(unsigned int i) const{
	return detectors[i];
}
int CosmicBench::get_non_ref_N() const{
	return non_ref_n;
}
double CosmicBench::get_z_Up() const{
	double z_Up = numeric_limits<double>::max();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if(Tomography::get_instance()->get_is_up((*it)->get_layer()) && (*it)->get_z()<z_Up){
			z_Up = (*it)->get_z();
		}
	}
	return z_Up;
}
double CosmicBench::get_z_Down() const{
	double z_Down = numeric_limits<double>::min();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if(!(Tomography::get_instance()->get_is_up((*it)->get_layer())) && (*it)->get_z()>z_Down){
			z_Down = (*it)->get_z();
		}
	}
	return z_Down;
}
