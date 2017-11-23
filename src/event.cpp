#define event_cpp
#include "event.h"

#include "Tanalyse_R.h"
#include "detector.h"
#include "cluster.h"

#include <limits>
#include <iostream>
#include <utility>
#include <sstream>
#include <algorithm>
#include <set>

#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>

using std::pair;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::flush;
using std::ostringstream;
using std::max_element;
using std::set;
using std::multimap;

using TMath::Min;
using TMath::Max;
using TMath::Abs;
using TMath::FloorNint;
using TMath::CeilNint;
using TMath::Sqrt;
/*
vector<map<double,int> > CosmicBenchEvent::combinaisons(map<double,int> sizes, bool allow_drop){
	map<double,int> partial_product;
	int current_product = 1;
	for(map<double,int>::iterator it = sizes.begin();it!=sizes.end();++it){
		if(it->second == 0) continue;
		partial_product[it->first] = current_product;
		current_product*=it->second;
	}
	vector<map<double,int> > result(current_product,map<double,int>());
	for(int i=0;i<current_product;i++){
		for(map<double,int>::iterator it = partial_product.begin();it!=partial_product.end();++it){
			result[i][it->first] = (i/(it->second)) % sizes[it->first];
		}
	}
	if(allow_drop){
		for(map<double,int>::iterator it = sizes.begin();it!=sizes.end();++it){
			if(it->second == 0) continue;
			map<double,int> current_sizes(sizes);
			current_sizes[it->first] = 0;
			vector<map<double,int> > current_result = combinaisons(current_sizes,false);
			result.insert(result.end(),current_result.begin(),current_result.end());
		}
	}
	return result;
}
*/
template<typename T>
vector<map<T,int> > CosmicBenchEvent::combinaisons(map<T,int> sizes, bool allow_drop){
	map<T,int> partial_product;
	int current_product = 1;
	for(typename map<T,int>::iterator it = sizes.begin();it!=sizes.end();++it){
		if(it->second == 0) continue;
		partial_product[it->first] = current_product;
		current_product*=it->second;
	}
	vector<map<T,int> > result(current_product,map<T,int>());
	for(int i=0;i<current_product;i++){
		for(typename map<T,int>::iterator it = partial_product.begin();it!=partial_product.end();++it){
			result[i][it->first] = (i/(it->second)) % sizes[it->first];
		}
	}
	if(allow_drop){
		for(typename map<T,int>::iterator it = sizes.begin();it!=sizes.end();++it){
			if(it->second == 0) continue;
			map<T,int> current_sizes(sizes);
			current_sizes[it->first] = 0;
			vector<map<T,int> > current_result = combinaisons(current_sizes,false);
			result.insert(result.end(),current_result.begin(),current_result.end());
		}
	}
	return result;
}
//template<> vector<map<double,int> > CosmicBenchEvent::combinaisons(map<double,int> sizes, bool allow_drop);
//template<> vector<map<int,int> > CosmicBenchEvent::combinaisons(map<int,int> sizes, bool allow_drop);

Event::Event(int evn_){
	evn = evn_;
	evttime = 0;
	has_spark = true;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	//if(detector) delete detector;
	detector = NULL;
}
Event::Event(const Event& other){
	evn = other.evn;
	evttime = other.evttime;
	type = other.type;
	has_spark = other.has_spark;
	strip_ampl = other.strip_ampl;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	for(vector<Cluster*>::const_iterator clus_it = other.clusters.begin();clus_it != other.clusters.end();++clus_it){
		clusters.push_back((*clus_it)->Clone());
	}
	//if(detector != NULL) delete detector;
	detector = (other.detector)->Clone();
}
Event& Event::operator=(const Event& other){
	evn = other.evn;
	evttime = other.evttime;
	type = other.type;
	has_spark = other.has_spark;
	strip_ampl = other.strip_ampl;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	for(vector<Cluster*>::const_iterator clus_it = other.clusters.begin();clus_it != other.clusters.end();++clus_it){
		clusters.push_back((*clus_it)->Clone());
	}
	if(detector != NULL) delete detector;
	detector = (other.detector)->Clone();
	return *this;
}
Event::Event(Tanalyse_R * treeObject, const Detector * const det,long entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	evttime = treeObject->evttime;
	has_spark = true;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//if(detector != NULL) delete detector;
	detector = det->Clone();
}
Event::Event(const Tanalyse_R * const treeObject, const Detector * const det){
	evn = treeObject->evn;
	evttime = treeObject->evttime;
	has_spark = true;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//if(detector != NULL) delete detector;
	detector = det->Clone();
}
Event::Event(const Detector * const detector_,int evn_, double evttime_){
	evn = evn_;
	evttime = evttime_;
	has_spark = true;
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//if(detector != NULL) delete detector;
	detector = detector_->Clone();
}
Event::~Event(){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	if(detector != NULL) delete detector;
}
int Event::get_evn() const{
	return evn;
}
double Event::get_evttime() const{
	return evttime;
}
Tomography::det_type Event::get_type() const{
	return type;
}
int Event::get_n_in_tree() const{
	return detector->get_n_in_tree();
}
bool Event::get_is_ref() const{
	return detector->get_is_ref();
}
bool Event::get_is_X() const{
	return detector->get_is_X();
}
double Event::get_StripPitch() const{
	return detector->get_StripPitch();
}
double Event::get_z() const{
	return detector->get_z();
}
int Event::get_NClus() const{
	return clusters.size();
}
void Event::do_cuts(){
	vector<Cluster*>::iterator clus_it = clusters.begin();
	while(clus_it!= clusters.end()){
		if(detector->is_suitable(*clus_it)) ++clus_it;
		else{
			delete *clus_it;
			clus_it = clusters.erase(clus_it);
		}
	}
}
vector<Cluster*> Event::get_clusters() const{
	vector<Cluster*> return_vector;
	for(vector<Cluster*>::const_iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		return_vector.push_back((*clus_it)->Clone());
	}
	return return_vector;
}
Detector * Event::get_det() const{
	return detector->Clone();
}

dummy_Event::dummy_Event(): Event(){
	type = Tomography::dummy;
	strip_ampl = vector<vector<double> >(Tomography::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
dummy_Event::dummy_Event(const dummy_Event& other): Event(other){
	type = Tomography::dummy;
}
dummy_Event& dummy_Event::operator=(const dummy_Event& other){	
	Event::operator=(other);
	type = Tomography::dummy;
	return *this;
}
dummy_Event::~dummy_Event(){

}
void dummy_Event::MultiCluster(){

}
void dummy_Event::ConvCluster(){

}
void dummy_Event::HoughCluster(int hole_nb){

}
void dummy_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=Tomography::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
TH1D * dummy_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_dummy_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1,0,0);
}
TH1D * dummy_Event::get_TOT_hist() const{
	ostringstream name;
	name << "ampl_dummy_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1,0,0);
}
Event * dummy_Event::Clone() const{
	return new dummy_Event(*this);
}

CM_Event::CM_Event(): Event(){
	type = Tomography::CM;
	detector = new CM_Detector();
	strip_ampl = vector<vector<double> >(CM_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
CM_Event::CM_Event(const CM_Event& other): Event(other){
	type = Tomography::CM;
}
CM_Event& CM_Event::operator=(const CM_Event& other){
	Event::operator=(other);
	type = Tomography::CM;
	return *this;
}
CM_Event::CM_Event(Tanalyse_R * treeObject, const CM_Detector * const det,long entry): Event(treeObject,det,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	strip_ampl = vector<vector<double> >(CM_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
	for(int i=0;i<treeObject->NClus[Tomography::CM][detector->get_n_in_tree()];i++){
		clusters.push_back(new CM_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	has_spark = (treeObject->Spark[Tomography::CM][detector->get_n_in_tree()]==1) ? true : false;
	type = Tomography::CM;
}
CM_Event::CM_Event(const Tanalyse_R * const treeObject, const CM_Detector * const det): Event(treeObject,det){
	for(int i=0;i<treeObject->NClus.find(Tomography::CM)->second[detector->get_n_in_tree()];i++){
		clusters.push_back(new CM_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	strip_ampl = vector<vector<double> >(CM_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
	has_spark = (treeObject->Spark.find(Tomography::CM)->second[detector->get_n_in_tree()]==1) ? true : false;
	type = Tomography::CM;
}
CM_Event::CM_Event(const CM_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_): Event(detector_,evn_, evttime_){
	if(strip_ampl_.size()!=CM_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
	type = Tomography::CM;
}
void CM_Event::MultiCluster(){
	// TODO : implement multicluster for CM
}
void CM_Event::ConvCluster(){
	// TODO : implement convcluster for CM
}
void CM_Event::HoughCluster(int hole_nb){
	// TODO : implement houghcluster for CM
}
void CM_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=CM_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
TH1D * CM_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_CM_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1024,-CM_Detector::size/2.,CM_Detector::size/2.);
}
TH1D * CM_Event::get_TOT_hist() const{
	ostringstream name;
	name << "TOT_CM_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1024,-CM_Detector::size/2.,CM_Detector::size/2.);
}
Event * CM_Event::Clone() const{
	return new CM_Event(*this);
}
CM_Event::~CM_Event(){

}

CM_Demux_Event::CM_Demux_Event(): Event(){
	type = Tomography::CM;
	strip_ampl = vector<vector<double> >(CM_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
CM_Demux_Event::CM_Demux_Event(const CM_Demux_Event& other): Event(other){
	type = Tomography::CM;
}
CM_Demux_Event& CM_Demux_Event::operator=(const CM_Demux_Event& other){
	Event::operator=(other);
	type = Tomography::CM;
	return *this;
}
CM_Demux_Event::CM_Demux_Event(const CM_Event& rawEvent): Event(rawEvent.detector,rawEvent.evn, rawEvent.evttime){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	if(dynamic_cast<CM_Detector*>(rawEvent.detector)->get_use_thin_strip()){
		for(vector<Cluster*>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(dynamic_cast<CM_Cluster*>(*it)->get_strip_type() == Tomography::Wide){
				for(vector<Cluster*>::const_iterator jt = rawEvent.clusters.begin();jt!=rawEvent.clusters.end();++jt){
					if(dynamic_cast<CM_Cluster*>(*jt)->get_strip_type() == Tomography::Thin){
						clusters.push_back(new CM_Demux_Cluster(*dynamic_cast<CM_Cluster*>(*jt),*dynamic_cast<CM_Cluster*>(*it)));
					}
				}
			}
		}
	}
	else{
		for(vector<Cluster*>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(dynamic_cast<CM_Cluster*>(*it)->get_strip_type() == Tomography::Wide){
				clusters.push_back(new CM_Demux_Cluster(*dynamic_cast<CM_Cluster*>(*it)));
			}
		}
	}
	has_spark = rawEvent.has_spark;
	strip_ampl = rawEvent.strip_ampl;
	type = Tomography::CM;
}
void CM_Demux_Event::MultiCluster(){
	// TODO : implement multicluster for CM
}
void CM_Demux_Event::ConvCluster(){
	// TODO : implement convcluster for CM
}
void CM_Demux_Event::HoughCluster(int hole_nb){
	// TODO : implement houghcluster for CM
}
void CM_Demux_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=64){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
Event * CM_Demux_Event::Clone() const{
	return new CM_Demux_Event(*this);
}
TH1D * CM_Demux_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_CM_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1024,-CM_Detector::size/2.,CM_Detector::size/2.);
}
TH1D * CM_Demux_Event::get_TOT_hist() const{
	ostringstream name;
	name << "TOT_CM_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	return new TH1D(name.str().c_str(),name.str().c_str(),1024,-CM_Detector::size/2.,CM_Detector::size/2.);
}
CM_Demux_Event::~CM_Demux_Event(){

}

MG_Event::MG_Event(): Event(){
	type = Tomography::MG;
	detector = new MG_Detector();
	strip_ampl = vector<vector<double> >(MG_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
MG_Event::MG_Event(const MG_Event& other): Event(other){
	type = Tomography::MG;
}
MG_Event& MG_Event::operator=(const MG_Event& other){
	Event::operator=(other);
	type = Tomography::MG;
	return *this;
}
MG_Event::MG_Event(Tanalyse_R * treeObject, const MG_Detector * const det,long entry): Event(treeObject,det,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	for(int i=0;i<treeObject->NClus[Tomography::MG][detector->get_n_in_tree()];i++){
		clusters.push_back(new MG_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	has_spark = (treeObject->Spark[Tomography::MG][detector->get_n_in_tree()]==1) ? true : false;
	type = Tomography::MG;
	strip_ampl = vector<vector<double> >(MG_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
MG_Event::MG_Event(const Tanalyse_R * const treeObject, const MG_Detector * const det): Event(treeObject,det){
	for(int i=0;i<treeObject->NClus.find(Tomography::MG)->second[detector->get_n_in_tree()];i++){
		clusters.push_back(new MG_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	has_spark = (treeObject->Spark.find(Tomography::MG)->second[detector->get_n_in_tree()]==1) ? true : false;
	type = Tomography::MG;
	strip_ampl = vector<vector<double> >(MG_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
MG_Event::MG_Event(const MG_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_): Event(detector_,evn_, evttime_){
	if(strip_ampl_.size()!=MG_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
	type = Tomography::MG;
}
void MG_Event::MultiCluster(){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	int p = 61;
	int n = 1024;
	map<int,bool> channelOverThreshold;
	map<int,StripInfo> allChannels;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=0;j<SampleMin;j++){
			current_strip.signal_sample[j] = false;
		}
		for(int j=SampleMax;j<Tomography::get_instance()->get_Nsample();j++){
			current_strip.signal_sample[j] = false;
		}
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				// time calculation with maximum
				/*
				if(j>0 && j<31){
					double a = (0.5*strip_ampl[i][j+1]) - strip_ampl[i][j] + (0.5*strip_ampl[i][j-1]);
					//double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					double b = (0.5*(strip_ampl[i][j+1] - strip_ampl[i][j-1])) - (2*a*j);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
				// --
			}
			if(strip_ampl[i][j]>(sigma*(detector->get_RMS(i)))){
				current_strip.TOT++;
				current_strip.signal_sample[j] = true;
			}
			else current_strip.signal_sample[j] = false;
		}
		// time calculation with rising edge
		
		int k=current_strip.MaxSample;
		double mean_xx = 0;
		double mean_xy = 0;
		double mean_x = 0;
		double mean_y = 0;
		int point_n = 0;
		/*
		while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k--;
		}
		*/
		/*
		while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
			k--;
		}
		k++;
		while(point_n<4 && strip_ampl[i][k]<current_strip.MaxAmpl){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k++;
			point_n++;
		}
		*/
		while(k>=SampleMin && strip_ampl[i][k]>0.9*current_strip.MaxAmpl){
			k--;
		}
		while(k>=SampleMin && strip_ampl[i][k]>0.1*current_strip.MaxAmpl){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k--;
			point_n++;
		}
		if(point_n>0){
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = k;
		
		// --
		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	vector<int> global_used_channel;
	unsigned int max_hole_size = detector->get_clustering_holes();
	while(global_used_channel.size()<channelOverThreshold.size()){
		pair<int,int> biggest_current_cluster(0,-1);
		vector<int> current_used_channel;
		for(int i=0;i<n;i++){
			if(channelOverThreshold.count(MG_Detector::StripToChannel_a[i])>0 && find(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel_a[i]) == global_used_channel.end()){
				pair<int,int> current_cluster(i,i);
				vector<int> used_channel;
				used_channel.push_back(MG_Detector::StripToChannel_a[i]);
				map<int,int> hole_channel;
				for(int j=i+1;j<n;j++){
					if(/*find(used_channel.begin(),used_channel.end(),MG_Detector::StripToChannel_a[j]) == used_channel.end() &&*/ find(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel_a[j]) == global_used_channel.end() /*&& !(hole_channel.count(MG_Detector::StripToChannel_a[j])>0)*/){
						if(channelOverThreshold.count(MG_Detector::StripToChannel_a[j])>0){
							current_cluster.second = j;
							used_channel.push_back(MG_Detector::StripToChannel_a[j]);
						}
						else if(hole_channel.size()<max_hole_size){
							hole_channel.insert(pair<int,int>(MG_Detector::StripToChannel_a[j],j));
						}
						else break;
					}		
					else break;
				}
				map<int,int>::iterator hole_it = hole_channel.begin();
				while(hole_it!=hole_channel.end()){
					if(hole_it->second > current_cluster.second){
						hole_channel.erase(hole_it);
						hole_it = hole_channel.begin();
					}
					else{
						++hole_it;
					}
				}
				for(hole_it = hole_channel.begin();hole_it!=hole_channel.end();++hole_it){
					used_channel.push_back(hole_it->first);
				}
				if((current_cluster.second - current_cluster.first)>(biggest_current_cluster.second - biggest_current_cluster.first)){
					biggest_current_cluster = current_cluster;
					current_used_channel = used_channel;
				}
			}
		}
		if(biggest_current_cluster.second != -1){
			cluster_list.push_back(biggest_current_cluster);
			global_used_channel.insert(global_used_channel.end(),current_used_channel.begin(),current_used_channel.end());
		}
	}
	//third loop : store the clusters and their caracteristics
	//NClus = cluster_list.size();
	for(unsigned int i=0;i<cluster_list.size();i++){
		/*
		TF1 * SRFfit;
		TGraphErrors * SRFgraph;
		int graph_point_n = 0;
		if(use_srf){
			SRFfit = new TF1("SRFfit",dynamic_cast<MG_Detector*>(detector),&MG_Detector::SRF_fit,0,1024,2,"MG_Detector","SRF_fit");
			SRFgraph = new TGraphErrors();
		}
		*/
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		//Micro TPC
		/*
		vector<double> pos_TPC(Tomography::get_instance()->get_Nsample(),0);
		vector<double> ampl_TPC(Tomography::get_instance()->get_Nsample(),0);
		vector<bool> used_sample_TPC(Tomography::get_instance()->get_Nsample(),false);
		double first_time = numeric_limits<double>::max();
		*/
		//--
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MG_Detector::StripToChannel_a[j]];
			double effective_ampl = current_strip.MaxAmpl/count(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel_a[j]);
			ClusPos = (ClusPos*ClusAmpl + j*effective_ampl)/(ClusAmpl + effective_ampl);
			ClusAmpl += effective_ampl;

			//Micro TPC
			/*
			if(current_strip.Time < first_time && current_strip.TOT>0) first_time = current_strip.Time;

			for(int k=SampleMin;k<SampleMax;k++){
				if(!(current_strip.signal_sample[k])) continue;
				double current_ampl = strip_ampl[MG_Detector::StripToChannel_a[j]][k];
				pos_TPC[k] = (pos_TPC[k]*ampl_TPC[k] + j*current_ampl)/(ampl_TPC[k]+current_ampl);
				ampl_TPC[k] += current_ampl;
				used_sample_TPC[k] = true;
			}
			*/
			// --

			if(effective_ampl>ClusMaxStripAmpl){
				ClusMaxStripAmpl = effective_ampl;
				ClusMaxSample = current_strip.MaxSample;
				ClusT = current_strip.Time;
				ClusMaxStrip = j;
				//ClusTOT[i] = current_strip.TOT;
			}
			if(current_strip.TOT>ClusTOT) ClusTOT = current_strip.TOT;
			/*
			if(use_srf){
				SRFgraph->SetPoint(graph_point_n,j,effective_ampl);
				SRFgraph->SetPointError(graph_point_n,0.5*MG_Detector::StripPitch,detector->get_RMS(MG_Detector::StripToChannel_a[j]));
				graph_point_n++;
			}
			*/
		}

		//Micro TPC
		/*
		double mean_xx = 0;
		double mean_xy = 0;
		double mean_x = 0;
		double mean_y = 0;
		unsigned short used_sample = 0;
		for(int k=SampleMin;k<SampleMax;k++){
			if(used_sample_TPC[k]){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*pos_TPC[k];
				mean_y += pos_TPC[k];
				used_sample++;
			}
		}
		if(used_sample>0){
			mean_xx /= used_sample;
			mean_xy /= used_sample;
			mean_y /= used_sample;
			mean_x /= used_sample;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			ClusPos = mean_y + (first_time - mean_x)*(mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			ClusT = first_time;
		}
		*/
		//--
		/*
		if(graph_point_n>2 && use_srf){
			for(int iPoint=0;iPoint<graph_point_n;iPoint++){
				double x_current = -1;
				double y_current = -1;
				SRFgraph->GetPoint(iPoint,x_current,y_current);
				y_current /= ClusMaxStripAmpl;
				//y_current = strip_ampl[MG_Detector::StripToChannel_a[static_cast<int>(x_current)]][static_cast<int>(ClusMaxSample)]/ClusMaxStripAmpl;
				SRFgraph->SetPoint(iPoint,x_current,y_current);
				SRFgraph->SetPointError(iPoint,MG_Detector::StripPitch,(SRFgraph->GetEY())[iPoint]/ClusMaxStripAmpl);
			}
			SRFfit->SetParameter(0,ClusPos);
			SRFfit->SetParLimits(0,ClusPos-0.5*ClusSize,ClusPos+0.5*ClusSize);
			SRFfit->SetParameter(1,0.25*ClusSize);
			SRFfit->SetParLimits(1,0,ClusSize);
			SRFgraph->Fit(SRFfit,"QN");
			
			if(detector.get_n_in_tree() == 0 && i == 0){
				TCanvas * cSRF = new TCanvas();
				TLine * posLine = new TLine(ClusPos[i],0,ClusPos[i],1);
				posLine->SetLineStyle(2);
				posLine->SetLineColor(4);
				cSRF->cd();
				SRFgraph->Draw("AP");
				SRFfit->Draw("same");
				posLine->Draw();
				cSRF->Modified();
				cSRF->Update();
				cout << graph_point_n << " " << SRFfit->GetChisquare() << endl;
				gROOT->GetApplication()->Run(true);
				delete posLine;
				delete cSRF;
			}
			
			ClusPos = SRFfit->GetParameter(0);
			//ClusSize = SRFfit->GetParameter(1);
		}
		if(use_srf){
			delete SRFfit; delete SRFgraph;
		}
		*/
		clusters.push_back(new MG_Cluster(detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
		//clusters.push_back(MG_Cluster(&detector,i,pos_TPC,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MG_Event::ConvCluster(){
	const double sigma_gaus = 2;
	TH1D * convHist = new TH1D("convHist","convHist",MG_Detector::Nstrip,0,MG_Detector::Nstrip);
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	map<int,StripInfo> allChannels;
	for(int i=0;i<MG_Detector::Nchannel;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				/*
				if(j>0 && j<31){
					double a = 0.5*strip_ampl[i][j+1] - 2*strip_ampl[i][j] + strip_ampl[i][j-1];
					double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
			}
			if(strip_ampl[i][j]>sigma*(detector->get_RMS(i))) current_strip.TOT++;
		}

		if(current_strip.TOT>2){
			int k=current_strip.MaxSample;
			double mean_xx = 0;
			double mean_xy = 0;
			double mean_x = 0;
			double mean_y = 0;
			while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*strip_ampl[i][k];
				mean_y += strip_ampl[i][k];
				k--;
			}
			int point_n = current_strip.MaxSample - k;
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = 0;

		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}


	for(int i=0;i<MG_Detector::Nstrip;i++){
		for(int j=0;j<MG_Detector::Nstrip;j++){
			convHist->Fill(i,(allChannels[MG_Detector::StripToChannel_a[j]].MaxAmpl)*TMath::Gaus(j,i,sigma_gaus));
		}
	}
	double current_derivative = 0;
	double old_derivative = 0;
	int clus_n = 0;
	for(int i=1;i<MG_Detector::Nstrip;i++){
		current_derivative = convHist->GetBinContent(i+1) - convHist->GetBinContent(i);
		if(old_derivative>0 && current_derivative<0){
			clusters.push_back(new MG_Cluster(detector,clus_n,i-1,sigma_gaus,convHist->GetBinContent(i),allChannels[MG_Detector::StripToChannel_a[i-1]].MaxSample,allChannels[MG_Detector::StripToChannel_a[i-1]].MaxAmpl,allChannels[MG_Detector::StripToChannel_a[i-1]].TOT,allChannels[MG_Detector::StripToChannel_a[i-1]].Time,i-1));
			clus_n++;
		}
		old_derivative = current_derivative;
	}
	delete convHist;

}
void MG_Event::HoughCluster(int hole_nb){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	int p = 61;
	int n = 1024;
	map<int,bool> channelOverThreshold;
	map<int,StripInfo> allChannels;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				/*
				if(j>0 && j<31){
					double a = 0.5*strip_ampl[i][j+1] - 2*strip_ampl[i][j] + strip_ampl[i][j-1];
					double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
			}
			if(strip_ampl[i][j]>sigma*(detector->get_RMS(i))) current_strip.TOT++;
		}

		if(current_strip.TOT>2){
			int k=current_strip.MaxSample;
			double mean_xx = 0;
			double mean_xy = 0;
			double mean_x = 0;
			double mean_y = 0;
			while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*strip_ampl[i][k];
				mean_y += strip_ampl[i][k];
				k--;
			}
			int point_n = current_strip.MaxSample - k;
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = 0;

		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	int k = 0;
	while(k<n){
		if(channelOverThreshold.count(MG_Detector::StripToChannel_a[k])>0){
			pair<int,int> current_cluster(k,k);
			unsigned int current_hole_nb = 0;
			for(int j=k+1;j<n;j++){
				if(channelOverThreshold.count(MG_Detector::StripToChannel_a[j])>0){
					current_cluster.second = j;
				}
				else if(current_hole_nb<hole_nb){
					current_hole_nb++;
				}
				else break;
			}
			if((1 + current_cluster.second - current_cluster.first) > (dynamic_cast<MG_Detector*>(detector)->get_ClusSizeCut_Min())) cluster_list.push_back(current_cluster);
			k = current_cluster.second+1;
		}
		else k++;
	}
	k=0;
	if(cluster_list.empty()){
		while(k<n){
			if(channelOverThreshold.count(MG_Detector::StripToChannel_a[k])>0){
				pair<int,int> current_cluster(k,k);
				unsigned int current_hole_nb = 0;
				for(int j=k+1;j<n;j++){
					if(channelOverThreshold.count(MG_Detector::StripToChannel_a[j])>0){
						current_cluster.second = j;
					}
					else if(current_hole_nb<hole_nb){
						current_hole_nb++;
					}
					else break;
				}
				cluster_list.push_back(current_cluster);
				k = current_cluster.second+1;
			}
			else k++;
		}
	}
	//third loop : store the clusters and their caracteristics
	//NClus = cluster_list.size();
	for(unsigned int i=0;i<cluster_list.size();i++){
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MG_Detector::StripToChannel_a[j]];
			double effective_ampl = current_strip.MaxAmpl;
			ClusPos = (ClusPos*ClusAmpl + j*effective_ampl)/(ClusAmpl + effective_ampl);
			ClusAmpl += effective_ampl;
			if(effective_ampl>ClusMaxStripAmpl){
				ClusMaxStripAmpl = effective_ampl;
				ClusMaxSample = current_strip.MaxSample;
				ClusT = current_strip.Time;
				ClusMaxStrip = j;
				//ClusTOT[i] = current_strip.TOT;
			}
			if(current_strip.TOT>ClusTOT) ClusTOT = current_strip.TOT;
		}

		clusters.push_back(new MG_Cluster(detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MG_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=MG_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
TH1D * MG_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_MG_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	TH1D * histo = new TH1D(name.str().c_str(),name.str().c_str(),1024,detector->get_offset() - MG_Detector::size/2.,detector->get_offset() + MG_Detector::size/2.);
	vector<pair<int,int> > cluster_edges;
	for(vector<Cluster*>::const_iterator it=clusters.begin();it!=clusters.end();++it){
		cluster_edges.push_back(pair<int,int>(FloorNint((*it)->get_pos()-(*it)->get_size()),CeilNint((*it)->get_pos()+(*it)->get_size())));
	}
	vector<bool> is_used(1024,false);
	for(vector<pair<int,int> >::iterator it = cluster_edges.begin();it!=cluster_edges.end();++it){
		if(it->first<0) it->first = 0;
		if(it->second<0) it->second = 0;
		if(it->first>1023) it->first = 1023;
		if(it->second>1023) it->second = 1023;
		for(int strip=it->first;strip<=it->second;strip++){
			if(is_used[strip]) continue;
			int channel = MG_Detector::StripToChannel_a[strip];
			histo->Fill(strip*MG_Detector::StripPitch + detector->get_offset() - MG_Detector::size/2.,*max_element(strip_ampl[channel].begin(),strip_ampl[channel].end()));
			is_used[strip] = true;
		}
	}
	return histo;
}
TH1D * MG_Event::get_TOT_hist() const{
	ostringstream name;
	name << "TOT_MG_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	TH1D * histo = new TH1D(name.str().c_str(),name.str().c_str(),1024,0,1024);
	vector<pair<int,int> > cluster_edges;
	
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int p = 61;
	int n = 1024;
	map<int,int> Channels_TOT;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>(sigma*(detector->get_RMS(i)))){
				current_strip.TOT++;
			}
		}
		Channels_TOT.insert(pair<int,int>(i,current_strip.TOT));
	}
	for(int i=0;i<n;i++){
		histo->Fill(i,Channels_TOT[MG_Detector::StripToChannel_a[i]]);
	}
	return histo;
}
Event * MG_Event::Clone() const{
	return new MG_Event(*this);
}
MG_Event::~MG_Event(){

}

MGv2_Event::MGv2_Event(): Event(){
	type = Tomography::MGv2;
	detector = new MGv2_Detector();
	strip_ampl = vector<vector<double> >(MGv2_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
}
MGv2_Event::MGv2_Event(const MGv2_Event& other): Event(other){
	type = Tomography::MGv2;
}
MGv2_Event& MGv2_Event::operator=(const MGv2_Event& other){
	Event::operator=(other);
	type = Tomography::MGv2;
	return *this;
}
MGv2_Event::MGv2_Event(Tanalyse_R * treeObject, const MGv2_Detector * const det,long entry): Event(treeObject,det,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	strip_ampl = vector<vector<double> >(MGv2_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
	for(int i=0;i<treeObject->NClus[Tomography::MGv2][detector->get_n_in_tree()];i++){
		clusters.push_back(new MGv2_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	has_spark = (treeObject->Spark[Tomography::MGv2][detector->get_n_in_tree()]==1) ? true : false;
	type = Tomography::MGv2;
}
MGv2_Event::MGv2_Event(const Tanalyse_R * const treeObject, const MGv2_Detector * const det): Event(treeObject,det){
	for(int i=0;i<treeObject->NClus.find(Tomography::MGv2)->second[detector->get_n_in_tree()];i++){
		clusters.push_back(new MGv2_Cluster(treeObject,i,det));
		if(!(det->is_suitable(clusters.back()))){
			delete clusters.back();
			clusters.pop_back();
		}
	}
	has_spark = (treeObject->Spark.find(Tomography::MGv2)->second[detector->get_n_in_tree()]==1) ? true : false;
	strip_ampl = vector<vector<double> >(MGv2_Detector::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
	type = Tomography::MGv2;
}
MGv2_Event::MGv2_Event(const MGv2_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_): Event(detector_,evn_,evttime_){
	if(strip_ampl_.size()!=MGv2_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
	type = Tomography::MGv2;
}
void MGv2_Event::MultiCluster(){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	double noise_RMS = Tomography::get_instance()->get_noise_RMS();
	int p = 61;
	int n = 1037;
	map<int,bool> channelOverThreshold;
	map<int,bool> noisyChannels;
	map<int,StripInfo> allChannels;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		if((detector->get_RMS(i))>noise_RMS){
			noisyChannels.insert(pair<int,bool>(i,true));
			//cout << "noisy channel : " << i << " in det : " << get_n_in_tree() << " with RMS : " << detector->get_RMS(i) << endl;
		}
		for(int j=0;j<SampleMin;j++){
			current_strip.signal_sample[j] = false;
		}
		for(int j=SampleMax;j<Tomography::get_instance()->get_Nsample();j++){
			current_strip.signal_sample[j] = false;
		}
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				// time calculation with maximum
				/*
				if(j>0 && j<31){
					double a = (0.5*strip_ampl[i][j+1]) - strip_ampl[i][j] + (0.5*strip_ampl[i][j-1]);
					//double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					double b = (0.5*(strip_ampl[i][j+1] - strip_ampl[i][j-1])) - (2*a*j);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
				// --
			}
			if(strip_ampl[i][j]>(sigma*(detector->get_RMS(i)))){
				current_strip.TOT++;
				current_strip.signal_sample[j] = true;
			}
			else current_strip.signal_sample[j] = false;
		}
		// time calculation with rising edge

		int k=current_strip.MaxSample;
		double mean_xx = 0;
		double mean_xy = 0;
		double mean_x = 0;
		double mean_y = 0;
		int point_n = 0;
		/*
		while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k--;
		}
		*/
		/*
		while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
			k--;
		}
		k++;
		while(point_n<4 && strip_ampl[i][k]<current_strip.MaxAmpl){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k++;
			point_n++;
		}
		*/
		while(k>=SampleMin && strip_ampl[i][k]>0.9*current_strip.MaxAmpl){
			k--;
		}
		while(k>=SampleMin && strip_ampl[i][k]>0.1*current_strip.MaxAmpl){
			mean_xx += k*k;
			mean_x += k;
			mean_xy += k*strip_ampl[i][k];
			mean_y += strip_ampl[i][k];
			k--;
			point_n++;
		}
		if(point_n>0){
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = k;
		
		// --
		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	vector<int> global_used_channel;
	unsigned int max_hole_size = detector->get_clustering_holes();
	while(global_used_channel.size()<channelOverThreshold.size()){
		pair<int,int> biggest_current_cluster(0,-1);
		vector<int> current_used_channel;
		for(int i=0;i<n;i++){
			if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[i])>0 && find(global_used_channel.begin(),global_used_channel.end(),MGv2_Detector::StripToChannel_a[i]) == global_used_channel.end()){
				pair<int,int> current_cluster(i,i);
				vector<int> used_channel;
				used_channel.push_back(MGv2_Detector::StripToChannel_a[i]);
				map<int,int> hole_channel;
				unsigned int noisy_channel_n = 0;
				for(int j=i+1;j<n;j++){
					if(/*find(used_channel.begin(),used_channel.end(),MGv2_Detector::StripToChannel_a[j]) == used_channel.end() &&*/ find(global_used_channel.begin(),global_used_channel.end(),MGv2_Detector::StripToChannel_a[j]) == global_used_channel.end() /*&& !(hole_channel.count(MGv2_Detector::StripToChannel_a[j])>0)*/){
						if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[j])>0){
							current_cluster.second = j;
							used_channel.push_back(MGv2_Detector::StripToChannel_a[j]);
						}
						else if(noisyChannels.count(MGv2_Detector::StripToChannel_a[j])>0){
							noisy_channel_n++;
							hole_channel.insert(pair<int,int>(MGv2_Detector::StripToChannel_a[j],j));
						}
						else if(hole_channel.size()<(max_hole_size + noisy_channel_n)){
							hole_channel.insert(pair<int,int>(MGv2_Detector::StripToChannel_a[j],j));
						}
						else break;
					}		
					else break;
				}
				//if(noisy_channel_n>0) cout << "make cluster with " << noisy_channel_n << " noisy channels" << endl;
				/*
				map<int,int>::iterator hole_it = hole_channel.begin();
				while(hole_it!=hole_channel.end()){
					if(hole_it->second > current_cluster.second){
						hole_channel.erase(hole_it);
						hole_it = hole_channel.begin();
					}
					else{
						++hole_it;
					}
				}
				for(hole_it = hole_channel.begin();hole_it!=hole_channel.end();++hole_it){
					used_channel.push_back(hole_it->first);
				}
				*/
				//--
				for(map<int,int>::iterator hole_it = hole_channel.begin();hole_it!=hole_channel.end();++hole_it){
					if(hole_it->second < current_cluster.second) used_channel.push_back(hole_it->first);
				}
				//--
				if((current_cluster.second - current_cluster.first)>(biggest_current_cluster.second - biggest_current_cluster.first)){
					biggest_current_cluster = current_cluster;
					current_used_channel = used_channel;
				}
			}
		}
		if(biggest_current_cluster.second != -1){
			cluster_list.push_back(biggest_current_cluster);
			global_used_channel.insert(global_used_channel.end(),current_used_channel.begin(),current_used_channel.end());
		}
	}
	//third loop : store the clusters and their caracteristics
	//NClus = cluster_list.size();
	for(unsigned int i=0;i<cluster_list.size();i++){
		/*
		TF1 * SRFfit;
		TGraphErrors * SRFgraph;
		int graph_point_n = 0;
		if(use_srf){
			SRFfit = new TF1("SRFfit",dynamic_cast<MGv2_Detector*>(detector),&MGv2_Detector::SRF_fit,0,1037,2,"MGv2_Detector","SRF_fit");
			SRFgraph = new TGraphErrors();
		}
		*/
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		//Micro TPC
		/*
		vector<double> pos_TPC(Tomography::get_instance()->get_Nsample(),0);
		vector<double> ampl_TPC(Tomography::get_instance()->get_Nsample(),0);
		vector<bool> used_sample_TPC(Tomography::get_instance()->get_Nsample(),false);
		double first_time = numeric_limits<double>::max();
		*/
		//--
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MGv2_Detector::StripToChannel_a[j]];
			double effective_ampl = current_strip.MaxAmpl/count(global_used_channel.begin(),global_used_channel.end(),MGv2_Detector::StripToChannel_a[j]);
			ClusPos = (ClusPos*ClusAmpl + j*effective_ampl)/(ClusAmpl + effective_ampl);
			ClusAmpl += effective_ampl;

			//Micro TPC
			/*
			if(current_strip.Time < first_time && current_strip.TOT>0) first_time = current_strip.Time;

			for(int k=SampleMin;k<SampleMax;k++){
				if(!(current_strip.signal_sample[k])) continue;
				double current_ampl = strip_ampl[MGv2_Detector::StripToChannel_a[j]][k];
				pos_TPC[k] = (pos_TPC[k]*ampl_TPC[k] + j*current_ampl)/(ampl_TPC[k]+current_ampl);
				ampl_TPC[k] += current_ampl;
				used_sample_TPC[k] = true;
			}
			*/
			// --

			if(effective_ampl>ClusMaxStripAmpl){
				ClusMaxStripAmpl = effective_ampl;
				ClusMaxSample = current_strip.MaxSample;
				ClusT = current_strip.Time;
				ClusMaxStrip = j;
				//ClusTOT[i] = current_strip.TOT;
			}
			if(current_strip.TOT>ClusTOT) ClusTOT = current_strip.TOT;
			/*
			if(use_srf){
				SRFgraph->SetPoint(graph_point_n,j,effective_ampl);
				SRFgraph->SetPointError(graph_point_n,0.5*MGv2_Detector::StripPitch,detector->get_RMS(MGv2_Detector::StripToChannel_a[j]));
				graph_point_n++;
			}
			*/
		}

		//Micro TPC
		/*
		double mean_xx = 0;
		double mean_xy = 0;
		double mean_x = 0;
		double mean_y = 0;
		unsigned short used_sample = 0;
		for(int k=SampleMin;k<SampleMax;k++){
			if(used_sample_TPC[k]){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*pos_TPC[k];
				mean_y += pos_TPC[k];
				used_sample++;
			}
		}
		if(used_sample>0){
			mean_xx /= used_sample;
			mean_xy /= used_sample;
			mean_y /= used_sample;
			mean_x /= used_sample;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			ClusPos = mean_y + (first_time - mean_x)*(mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			ClusT = first_time;
		}
		*/
		//--

		/*
		if(graph_point_n>2 && use_srf){
			for(int iPoint=0;iPoint<graph_point_n;iPoint++){
				double x_current = -1;
				double y_current = -1;
				SRFgraph->GetPoint(iPoint,x_current,y_current);
				y_current /= ClusMaxStripAmpl;
				//y_current = strip_ampl[MGv2_Detector::StripToChannel_a[static_cast<int>(x_current)]][static_cast<int>(ClusMaxSample)]/ClusMaxStripAmpl;
				SRFgraph->SetPoint(iPoint,x_current,y_current);
				SRFgraph->SetPointError(iPoint,MGv2_Detector::StripPitch,(SRFgraph->GetEY())[iPoint]/ClusMaxStripAmpl);
			}
			SRFfit->SetParameter(0,ClusPos);
			SRFfit->SetParLimits(0,ClusPos-0.5*ClusSize,ClusPos+0.5*ClusSize);
			SRFfit->SetParameter(1,0.25*ClusSize);
			SRFfit->SetParLimits(1,0,ClusSize);
			SRFgraph->Fit(SRFfit,"QN");
			
			if(detector.get_n_in_tree() == 0 && i == 0){
				TCanvas * cSRF = new TCanvas();
				TLine * posLine = new TLine(ClusPos[i],0,ClusPos[i],1);
				posLine->SetLineStyle(2);
				posLine->SetLineColor(4);
				cSRF->cd();
				SRFgraph->Draw("AP");
				SRFfit->Draw("same");
				posLine->Draw();
				cSRF->Modified();
				cSRF->Update();
				cout << graph_point_n << " " << SRFfit->GetChisquare() << endl;
				gROOT->GetApplication()->Run(true);
				delete posLine;
				delete cSRF;
			}
			
			ClusPos = SRFfit->GetParameter(0);
			//ClusSize = SRFfit->GetParameter(1);
		}
		if(use_srf){
			delete SRFfit; delete SRFgraph;
		}
		*/
		clusters.push_back(new MGv2_Cluster(detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
		//clusters.push_back(MGv2_Cluster(&detector,i,pos_TPC,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MGv2_Event::ConvCluster(){
	const double sigma_gaus = 2;
	TH1D * convHist = new TH1D("convHist","convHist",MGv2_Detector::Nstrip,0,MGv2_Detector::Nstrip);
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	map<int,StripInfo> allChannels;
	for(int i=0;i<MGv2_Detector::Nchannel;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				/*
				if(j>0 && j<31){
					double a = 0.5*strip_ampl[i][j+1] - 2*strip_ampl[i][j] + strip_ampl[i][j-1];
					double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
			}
			if(strip_ampl[i][j]>sigma*(detector->get_RMS(i))) current_strip.TOT++;
		}

		if(current_strip.TOT>2){
			int k=current_strip.MaxSample;
			double mean_xx = 0;
			double mean_xy = 0;
			double mean_x = 0;
			double mean_y = 0;
			while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*strip_ampl[i][k];
				mean_y += strip_ampl[i][k];
				k--;
			}
			int point_n = current_strip.MaxSample - k;
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = 0;

		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}


	for(int i=0;i<MGv2_Detector::Nstrip;i++){
		for(int j=0;j<MGv2_Detector::Nstrip;j++){
			convHist->Fill(i,(allChannels[MGv2_Detector::StripToChannel_a[j]].MaxAmpl)*TMath::Gaus(j,i,sigma_gaus));
		}
	}
	double current_derivative = 0;
	double old_derivative = 0;
	int clus_n = 0;
	for(int i=1;i<MGv2_Detector::Nstrip;i++){
		current_derivative = convHist->GetBinContent(i+1) - convHist->GetBinContent(i);
		if(old_derivative>0 && current_derivative<0){
			clusters.push_back(new MGv2_Cluster(detector,clus_n,i-1,sigma_gaus,convHist->GetBinContent(i),allChannels[MGv2_Detector::StripToChannel_a[i-1]].MaxSample,allChannels[MGv2_Detector::StripToChannel_a[i-1]].MaxAmpl,allChannels[MGv2_Detector::StripToChannel_a[i-1]].TOT,allChannels[MGv2_Detector::StripToChannel_a[i-1]].Time,i-1));
			clus_n++;
		}
		old_derivative = current_derivative;
	}
	delete convHist;

}
void MGv2_Event::HoughCluster(int hole_nb){
	for(vector<Cluster*>::iterator clus_it = clusters.begin();clus_it != clusters.end();++clus_it){
		delete *clus_it;
	}
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int TOTCut = Tomography::get_instance()->get_TOTCut();
	int p = 61;
	int n = 1037;
	map<int,bool> channelOverThreshold;
	map<int,StripInfo> allChannels;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				/*
				if(j>0 && j<31){
					double a = 0.5*strip_ampl[i][j+1] - 2*strip_ampl[i][j] + strip_ampl[i][j-1];
					double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				*/
			}
			if(strip_ampl[i][j]>sigma*(detector->get_RMS(i))) current_strip.TOT++;
		}

		if(current_strip.TOT>2){
			int k=current_strip.MaxSample;
			double mean_xx = 0;
			double mean_xy = 0;
			double mean_x = 0;
			double mean_y = 0;
			while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector->get_RMS(i)))){
				mean_xx += k*k;
				mean_x += k;
				mean_xy += k*strip_ampl[i][k];
				mean_y += strip_ampl[i][k];
				k--;
			}
			int point_n = current_strip.MaxSample - k;
			mean_xx /= point_n;
			mean_xy /= point_n;
			mean_y /= point_n;
			mean_x /= point_n;
			//double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
			//double intercept = mean_y - slope*mean_x;
			current_strip.Time = mean_x - mean_y*(mean_xx - mean_x*mean_x)/(mean_xy - mean_x*mean_y);
		}
		else current_strip.Time = 0;

		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	int k = 0;
	while(k<n){
		if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[k])>0){
			pair<int,int> current_cluster(k,k);
			unsigned int current_hole_nb = 0;
			for(int j=k+1;j<n;j++){
				if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[j])>0){
					current_cluster.second = j;
				}
				else if(current_hole_nb<hole_nb){
					current_hole_nb++;
				}
				else break;
			}
			if((1 + current_cluster.second - current_cluster.first) > (dynamic_cast<MGv2_Detector*>(detector)->get_ClusSizeCut_Min())) cluster_list.push_back(current_cluster);
			k = current_cluster.second+1;
		}
		else k++;
	}
	k=0;
	if(cluster_list.empty()){
		while(k<n){
			if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[k])>0){
				pair<int,int> current_cluster(k,k);
				unsigned int current_hole_nb = 0;
				for(int j=k+1;j<n;j++){
					if(channelOverThreshold.count(MGv2_Detector::StripToChannel_a[j])>0){
						current_cluster.second = j;
					}
					else if(current_hole_nb<hole_nb){
						current_hole_nb++;
					}
					else break;
				}
				cluster_list.push_back(current_cluster);
				k = current_cluster.second+1;
			}
			else k++;
		}
	}
	//third loop : store the clusters and their caracteristics
	//NClus = cluster_list.size();
	for(unsigned int i=0;i<cluster_list.size();i++){
		/*
		TF1 * SRFfit = new TF1("SRFfit",dynamic_cast<MGv2_Detector*>(detector),&MGv2_Detector::SRF_fit,0,1037,2,"MGv2_Detector","SRF_fit");
		TGraphErrors * SRFgraph = new TGraphErrors();
		int graph_point_n = 0;
		*/
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MGv2_Detector::StripToChannel_a[j]];
			double effective_ampl = current_strip.MaxAmpl;
			ClusPos = (ClusPos*ClusAmpl + j*effective_ampl)/(ClusAmpl + effective_ampl);
			ClusAmpl += effective_ampl;
			if(effective_ampl>ClusMaxStripAmpl){
				ClusMaxStripAmpl = effective_ampl;
				ClusMaxSample = current_strip.MaxSample;
				ClusT = current_strip.Time;
				ClusMaxStrip = j;
				//ClusTOT[i] = current_strip.TOT;
			}
			if(current_strip.TOT>ClusTOT) ClusTOT = current_strip.TOT;
			/*
			SRFgraph->SetPoint(graph_point_n,j,effective_ampl);
			SRFgraph->SetPointError(graph_point_n,0.5*MGv2_Detector::StripPitch,detector->get_RMS(MGv2_Detector::StripToChannel_a[j]));
			graph_point_n++;
			*/
		}
		/*
		if(graph_point_n>2 && use_srf){
			for(int iPoint=0;iPoint<graph_point_n;iPoint++){
				double x_current = -1;
				double y_current = -1;
				SRFgraph->GetPoint(iPoint,x_current,y_current);
				y_current /= ClusMaxStripAmpl;
				//y_current = strip_ampl[MGv2_Detector::StripToChannel_a[static_cast<int>(x_current)]][static_cast<int>(ClusMaxSample)]/ClusMaxStripAmpl;
				SRFgraph->SetPoint(iPoint,x_current,y_current);
				SRFgraph->SetPointError(iPoint,MGv2_Detector::StripPitch,(SRFgraph->GetEY())[iPoint]/ClusMaxStripAmpl);
			}
			SRFfit->SetParameter(0,ClusPos);
			SRFfit->SetParLimits(0,ClusPos-0.5*ClusSize,ClusPos+0.5*ClusSize);
			SRFfit->SetParameter(1,0.25*ClusSize);
			SRFfit->SetParLimits(1,0,ClusSize);
			SRFgraph->Fit(SRFfit,"QN");
			
			if(detector.get_n_in_tree() == 0 && i == 0){
				TCanvas * cSRF = new TCanvas();
				TLine * posLine = new TLine(ClusPos[i],0,ClusPos[i],1);
				posLine->SetLineStyle(2);
				posLine->SetLineColor(4);
				cSRF->cd();
				SRFgraph->Draw("AP");
				SRFfit->Draw("same");
				posLine->Draw();
				cSRF->Modified();
				cSRF->Update();
				cout << graph_point_n << " " << SRFfit->GetChisquare() << endl;
				gROOT->GetApplication()->Run(true);
				delete posLine;
				delete cSRF;
			}
			
			ClusPos = SRFfit->GetParameter(0);
			ClusSize = SRFfit->GetParameter(1);
		}
		delete SRFfit; delete SRFgraph;
		*/
		clusters.push_back(new MGv2_Cluster(detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MGv2_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=MGv2_Detector::Nchannel){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
TH1D * MGv2_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_MGv2_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	TH1D * histo = new TH1D(name.str().c_str(),name.str().c_str(),1037,detector->get_offset() - MGv2_Detector::size/2.,detector->get_offset() + MGv2_Detector::size/2.);
	vector<pair<int,int> > cluster_edges;
	for(vector<Cluster*>::const_iterator it=clusters.begin();it!=clusters.end();++it){
		cluster_edges.push_back(pair<int,int>(FloorNint((*it)->get_pos()-(*it)->get_size()),CeilNint((*it)->get_pos()+(*it)->get_size())));
	}
	vector<bool> is_used(1037,false);
	for(vector<pair<int,int> >::iterator it = cluster_edges.begin();it!=cluster_edges.end();++it){
		if(it->first<0) it->first = 0;
		if(it->second<0) it->second = 0;
		if(it->first>1036) it->first = 1036;
		if(it->second>1036) it->second = 1036;
		for(int strip=it->first;strip<=it->second;strip++){
			if(is_used[strip]) continue;
			int channel = MGv2_Detector::StripToChannel_a[strip];
			histo->Fill(strip*MGv2_Detector::StripPitch + detector->get_offset() - MGv2_Detector::size/2.,*max_element(strip_ampl[channel].begin(),strip_ampl[channel].end()));
			is_used[strip] = true;
		}
	}
	return histo;
}
TH1D * MGv2_Event::get_TOT_hist() const{
	ostringstream name;
	name << "TOT_MGv2_det_" << detector->get_n_in_tree() << "_evn_" << evn;
	TH1D * histo = new TH1D(name.str().c_str(),name.str().c_str(),1037,0,1037);
	vector<pair<int,int> > cluster_edges;
	
	double sigma = Tomography::get_instance()->get_sigma();
	int SampleMin = Tomography::get_instance()->get_SampleMin();
	int SampleMax = Tomography::get_instance()->get_SampleMax();
	int p = 61;
	int n = 1037;
	map<int,int> Channels_TOT;
	for(int i=0;i<p;i++){
		StripInfo current_strip;
		current_strip.MaxAmpl = 0;
		current_strip.MaxSample = 0;
		current_strip.TOT = 0;
		current_strip.Time = 0;
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>(sigma*(detector->get_RMS(i)))){
				current_strip.TOT++;
			}
		}
		Channels_TOT.insert(pair<int,int>(i,current_strip.TOT));
	}
	for(int i=0;i<n;i++){
		histo->Fill(i,Channels_TOT[MG_Detector::StripToChannel_a[i]]);
	}
	return histo;
}
Event * MGv2_Event::Clone() const{
	return new MGv2_Event(*this);
}
MGv2_Event::~MGv2_Event(){

}

CosmicBenchEvent::CosmicBenchEvent(){
	rayPairs.clear();
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	detectors = NULL;
	events.clear();
	evn = -1;
	evttime = 0;
}
CosmicBenchEvent::CosmicBenchEvent(const CosmicBenchEvent& other){
	evn = other.evn;
	evttime = other.evttime;
	detectors = other.detectors;
	rayPairs.clear();
	rayPairs.assign(other.rayPairs.begin(),other.rayPairs.end());
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
	for(vector<Event*>::const_iterator it = other.events.begin();it!= other.events.end();++it){
		events.push_back((*it)->Clone());
	}
}
CosmicBenchEvent& CosmicBenchEvent::operator=(const CosmicBenchEvent& other){
	evn = other.evn;
	evttime = other.evttime;
	detectors = other.detectors;
	rayPairs.clear();
	rayPairs.assign(other.rayPairs.begin(),other.rayPairs.end());
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
	for(vector<Event*>::const_iterator it = other.events.begin();it!= other.events.end();++it){
		events.push_back((*it)->Clone());
	}
	return *this;
}
CosmicBenchEvent::CosmicBenchEvent(const CosmicBench * const detectors_, Tanalyse_R * treeObject, long entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	detectors = detectors_;
	evn = treeObject->evn;
	evttime = treeObject->evttime;
	rayPairs.clear();
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
	int det_N = detectors->get_det_N_tot();
	for(int i=0;i<det_N;i++){
		events.push_back(detectors->get_detector(i)->build_event(treeObject));
	}
}
CosmicBenchEvent::CosmicBenchEvent(const CosmicBench * const detectors_, const Tanalyse_R * const treeObject){
	evn = treeObject->evn;
	evttime = treeObject->evttime;
	rayPairs.clear();
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
	int det_N = detectors->get_det_N_tot();
	for(int i=0;i<det_N;i++){
		events.push_back(detectors->get_detector(i)->build_event(treeObject));
	}
	detectors = detectors_;
}
CosmicBenchEvent::CosmicBenchEvent(const CosmicBench * const detectors_,const vector<Event*> events_){
	rayPairs.clear();
	events.clear();
	unsigned int det_N = detectors->get_det_N_tot();
	if(events_.size()!=det_N){
		cout << "problem in event size" << endl;
		return;
	}
	evn = (events_.front())->get_evn();
	evttime = (events_.front())->get_evttime();
	set<pair<Tomography::det_type,int> > det_is_used;
	for(vector<Event*>::const_iterator it = events_.begin();it!=events_.end();++it){
		if((*it)->get_evn() != evn){
			cout << "attempt to merge event with different number in cosmicbench event" << endl;
			return;
		}
		if(det_is_used.count(pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree()))){
			cout << "problem in events" << endl;
			return;
		}
		det_is_used.insert(pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree()));
		events.push_back((*it)->Clone());
	}
	detectors = detectors_;
}
CosmicBenchEvent::~CosmicBenchEvent(){
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
	rayPairs.clear();
}
void CosmicBenchEvent::Demux_CM(){
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		if((*it)->get_type() == Tomography::CM){
			CM_Demux_Event * currentDemux = new CM_Demux_Event(*dynamic_cast<CM_Event*>(*it));
			delete *it;
			*it = new CM_Demux_Event(*currentDemux);
			delete currentDemux;
		}
	}
}
void CosmicBenchEvent::createPairs(){
	this->Demux_CM();
	double chiSquare_threshold = numeric_limits<double>::max();
	map<bool, map<bool, map<int,vector<Cluster*> > > > currentClusters;
	map<bool, map<bool, map<int,int> > > sizes;
	double max_z = numeric_limits<double>::max();
	double min_z = numeric_limits<double>::min();
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		double z;
		bool is_up;
		bool is_X;
		int layer;
		if((*it)->get_is_ref()){
			vector<Cluster*> temp_clusters((*it)->get_clusters());
			for(vector<Cluster*>::iterator jt=temp_clusters.begin();jt!=temp_clusters.end();++jt){
				z = (*jt)->get_z();
				is_X = (*jt)->get_is_X();
				layer = (*jt)->get_layer();
				is_up = Tomography::get_instance()->get_is_up(layer);
				if(is_up && z<max_z) max_z = z;
				if((!is_up) && z>min_z) min_z = z;
				currentClusters[is_up][is_X][layer].push_back((*jt)->Clone());
				sizes[is_up][is_X][layer]++;
				delete *jt;
			}
		}
	}

	if(currentClusters[true][true].size() == 2 && currentClusters[true][false].size() == 2 && currentClusters[false][true].size() == 2 && currentClusters[false][false].size() == 2){
		vector<RayPair> suitableRays;
		map<bool, map<bool, vector<map<int,int> > > > comb;
		comb[true][true] = combinaisons(sizes[true][true]);
		comb[true][false] = combinaisons(sizes[true][false]);
		comb[false][true] = combinaisons(sizes[false][true]);
		comb[false][false] = combinaisons(sizes[false][false]);
		map<bool,map<bool,map<int,Ray_2D> > > possibleRay_2D;
		map<bool,map<bool,multimap<pair<int,int>,int> > > clus_to_ray;
		map<bool,map<bool,multimap<int,pair<int,int> > > > ray_to_clus;
		for(map<bool, map<bool, vector<map<int,int> > > >::iterator it = comb.begin();it!=comb.end();++it){
			for(map<bool, vector<map<int,int> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
				for(unsigned int i=0;i<(jt->second).size();i++){
					possibleRay_2D[it->first][jt->first][i] = Ray_2D(Ray_2D((jt->first) ? 'X' : 'Y'));
					for(map<int,vector<Cluster*> >::iterator clus_it = currentClusters[it->first][jt->first].begin();clus_it!= currentClusters[it->first][jt->first].end();++clus_it){
						possibleRay_2D[it->first][jt->first][i].add_cluster(clus_it->second[(jt->second)[i][clus_it->first]]);
					}
					possibleRay_2D[it->first][jt->first][i].process();
					if(possibleRay_2D[it->first][jt->first][i].get_chiSquare()<0 || possibleRay_2D[it->first][jt->first][i].get_chiSquare()>chiSquare_threshold){
						possibleRay_2D[it->first][jt->first].erase(i);
						continue;
					}
					for(map<int,int>::iterator id_it = (jt->second)[i].begin();id_it!=(jt->second)[i].end();++id_it){
						clus_to_ray[it->first][jt->first].insert(pair<pair<double,int>,int>(*id_it,i));
						ray_to_clus[it->first][jt->first].insert(pair<int, pair<double,int> >(i,*id_it));
					}
				}
			}
		}

		bool b=true;
		while(b){
			b = false;
			double bestDoca = 1000000;
			double current_chiSquare = chiSquare_threshold;
			map<bool, map<bool, int > > best_comb;
			RayPair bestRay = RayPair();

			
			for(map<int,Ray_2D>::iterator HX_it = possibleRay_2D[true][true].begin();HX_it!=possibleRay_2D[true][true].end();++HX_it){
				for(map<int,Ray_2D>::iterator HY_it = possibleRay_2D[true][false].begin();HY_it!=possibleRay_2D[true][false].end();++HY_it){
					for(map<int,Ray_2D>::iterator BX_it = possibleRay_2D[false][true].begin();BX_it!=possibleRay_2D[false][true].end();++BX_it){
						for(map<int,Ray_2D>::iterator BY_it = possibleRay_2D[false][false].begin();BY_it!=possibleRay_2D[false][false].end();++BY_it){
							RayPair currentRayPair = RayPair(Ray(HX_it->second,HY_it->second),Ray(BX_it->second,BY_it->second));
							currentRayPair.process();
							double currentDoca = currentRayPair.get_doca();
							Point currentPoCA = currentRayPair.get_PoCA();
							double distance_z = 0;
							double distance_x = 0;
							double distance_y = 0;
							//if(currentPoCA.get_Z()>(max_z+10) || currentPoCA.get_Z()<(min_z-10)) continue;
							if(currentPoCA.get_Z()>max_z) distance_z = currentPoCA.get_Z() - max_z;
							if(currentPoCA.get_Z()<min_z) distance_z = -(currentPoCA.get_Z()) + min_z;
							else distance_z= 0;
							//if(currentPoCA.get_X()>10.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_X()<-10.*Tomography::get_instance()->get_XY_size()/10.) continue;
							if(currentPoCA.get_X()>(Tomography::get_instance()->get_XY_size()/2.)) distance_x = currentPoCA.get_X() - (Tomography::get_instance()->get_XY_size()/2.);
							else if(currentPoCA.get_X()<(-Tomography::get_instance()->get_XY_size()/2.)) distance_x = -(currentPoCA.get_X()) - (Tomography::get_instance()->get_XY_size()/2.);
							else distance_x = 0;
							//if(currentPoCA.get_Y()>10.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_Y()<-10.*Tomography::get_instance()->get_XY_size()/10.) continue;
							if(currentPoCA.get_Y()>(Tomography::get_instance()->get_XY_size()/2.)) distance_y = currentPoCA.get_Y() - (Tomography::get_instance()->get_XY_size()/2.);
							else if(currentPoCA.get_Y()<(-Tomography::get_instance()->get_XY_size()/2.)) distance_y = -(currentPoCA.get_Y()) - (Tomography::get_instance()->get_XY_size()/2.);
							else distance_y = 0;
							currentDoca += Sqrt((distance_x*distance_x)+(distance_y*distance_y)+(distance_z*distance_z));
							if(currentDoca<bestDoca){
								bestDoca = currentDoca;
								bestRay = currentRayPair;
								best_comb[true][true] = HX_it->first;
								best_comb[true][false] = HY_it->first;
								best_comb[false][true] = BX_it->first;
								best_comb[false][false] = BY_it->first;
								b = true;
							}
						}
					}
				}
			}

			if(b){
				suitableRays.push_back(bestRay);
				for(map<bool, map<bool, int> >::iterator it = best_comb.begin();it!=best_comb.end();++it){
					for(map<bool, int>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
						pair<multimap<int,pair<int,int> >::iterator,multimap<int,pair<int,int> >::iterator > clus_of_best_ray = ray_to_clus[it->first][jt->first].equal_range(best_comb[it->first][jt->first]);
						for(multimap<int,pair<int,int> >::iterator clus_it=clus_of_best_ray.first;clus_it!=clus_of_best_ray.second;++clus_it){
							pair<multimap<pair<int,int>,int >::iterator,multimap<pair<int,int>,int >::iterator > ray_using_best_clus = clus_to_ray[it->first][jt->first].equal_range(clus_it->second);
							for(multimap<pair<int,int>,int>::iterator ray_it=ray_using_best_clus.first;ray_it!=ray_using_best_clus.second;++ray_it){
								if(possibleRay_2D[it->first][jt->first].find(ray_it->second)!=possibleRay_2D[it->first][jt->first].end()) possibleRay_2D[it->first][jt->first].erase(ray_it->second);
							}
							clus_to_ray[it->first][jt->first].erase(ray_using_best_clus.first,ray_using_best_clus.second);
						}
						ray_to_clus[it->first][jt->first].erase(clus_of_best_ray.first,clus_of_best_ray.second);
					}
				}
			}
			if(possibleRay_2D[true][true].size() < 1) b = false;
			if(possibleRay_2D[true][false].size() < 1) b = false;
			if(possibleRay_2D[false][true].size() < 1) b = false;
			if(possibleRay_2D[false][false].size() < 1) b = false;
		}
		rayPairs = suitableRays;
	}
	//---------------------------------------------------------------------------------------------------------------------
	else{
		map<bool, map<bool, vector<Ray_2D> > > suitableRays;
		//compute for both and down sensitive areas
		for(map<bool, map<bool, map<int,vector<Cluster*> > > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
			//compute for both coordinates
			for(map<bool, map<int,vector<Cluster*> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
				bool b = true;
				
				//get size
				//for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				//	if((kt->second).size()==0) b = false;
				//	sizes[it->first][jt->first][kt->first] = (kt->second).size();
				//}
				
				//find the biggest number of good clusters combinaisons
				while(b && (jt->second).size()>1){
					b = false;
					//find best combinaison of clusters
					vector<map<int,int> > comb = combinaisons(sizes[it->first][jt->first]);
					double current_chiSquare = chiSquare_threshold;
					map<int,int> best_comb;
					char coord = (jt->first) ? 'X' : 'Y';
					Ray_2D bestRay = Ray_2D(coord);
					for(vector<map<int,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
						//try a comb
						Ray_2D currentRay = Ray_2D(coord);
						for(map<int,vector<Cluster*> >::iterator nt = (jt->second).begin();nt!= (jt->second).end();++nt){
							currentRay.add_cluster(nt->second[(*kt)[nt->first]]);
						}
						currentRay.process();
						if(currentRay.get_chiSquare()<current_chiSquare && currentRay.get_chiSquare()>-1){
							bestRay = Ray_2D(currentRay);
							current_chiSquare = currentRay.get_chiSquare();
							best_comb = *kt;
							b = true;
						}
					}
					if(b){
						suitableRays[it->first][jt->first].push_back(Ray_2D(bestRay));
						for(map<int,int>::iterator kt = best_comb.begin();kt!=best_comb.end();++kt){
							delete (jt->second)[kt->first][kt->second];
							(jt->second)[kt->first].erase((jt->second)[kt->first].begin()+kt->second);
							sizes[it->first][jt->first][kt->first]--;
							if(sizes[it->first][jt->first][kt->first]<1){
								//b = false;
								(jt->second).erase((jt->second).find(kt->first));
								sizes[it->first][jt->first].erase(sizes[it->first][jt->first].find(kt->first));
							}
						}
					}
				}
			}
		}

		//cout << Max(Max(suitableRays[true][true].size(),suitableRays[true][false].size()),Max(suitableRays[false][true].size(),suitableRays[false][false].size())) << endl;
		int min_size = Min(Min(suitableRays[true][true].size(),suitableRays[true][false].size()),Min(suitableRays[false][true].size(),suitableRays[false][false].size()));
		
		while(min_size>0){
			double bestDoca = 50;//numeric_limits<double>::max();
			vector<Ray_2D>::iterator best_it = suitableRays[true][true].end();
			vector<Ray_2D>::iterator best_jt = suitableRays[true][false].end();
			vector<Ray_2D>::iterator best_kt = suitableRays[false][true].end();
			vector<Ray_2D>::iterator best_lt = suitableRays[false][false].end();
			for(vector<Ray_2D>::iterator it = suitableRays[true][true].begin(); it!= suitableRays[true][true].end();++it){
				for(vector<Ray_2D>::iterator jt = suitableRays[true][false].begin(); jt!= suitableRays[true][false].end();++jt){
					for(vector<Ray_2D>::iterator kt = suitableRays[false][true].begin(); kt!= suitableRays[false][true].end();++kt){
						for(vector<Ray_2D>::iterator lt = suitableRays[false][false].begin(); lt!= suitableRays[false][false].end();++lt){
							RayPair currentRayPair(Ray(*it,*jt),Ray(*kt,*lt));
							double currentDoca = currentRayPair.get_doca();
							Point currentPoCA = currentRayPair.get_PoCA();
							if(currentPoCA.get_Z()>max_z || currentPoCA.get_Z()<min_z) continue;
							if(currentPoCA.get_X()>6.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_X()<-6.*Tomography::get_instance()->get_XY_size()/10.) continue;
							if(currentPoCA.get_Y()>6.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_Y()<-6.*Tomography::get_instance()->get_XY_size()/10.) continue;
							
							if(currentDoca<bestDoca){
								bestDoca = currentDoca;
								best_it = vector<Ray_2D>::iterator(it);
								best_jt = vector<Ray_2D>::iterator(jt);
								best_kt = vector<Ray_2D>::iterator(kt);
								best_lt = vector<Ray_2D>::iterator(lt);
							}
						}
					}
				}
			}
			if(best_it == suitableRays[true][true].end() || best_jt == suitableRays[true][false].end() || best_kt == suitableRays[false][true].end() || best_lt == suitableRays[false][false].end()) break;
			rayPairs.push_back(RayPair(Ray(*best_it,*best_jt),Ray(*best_kt,*best_lt)));
			rayPairs.back().upRay.angle_correction();
			rayPairs.back().downRay.angle_correction();
			suitableRays[true][true].erase(best_it);
			suitableRays[true][false].erase(best_jt);
			suitableRays[false][true].erase(best_kt);
			suitableRays[false][false].erase(best_lt);
			min_size--;
		}
	}

	for(map<bool, map<bool, map<int,vector<Cluster*> > > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		for(map<bool, map<int,vector<Cluster*> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(map<int,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				for(unsigned int i=0;i<(kt->second).size();i++){
					delete kt->second[i];
				}
			}
		}
	}
}
/*
void CosmicBenchEvent::createPairs(){
	this->Demux_CM();
	double chiSquare_threshold = numeric_limits<double>::max();
	map<bool, map<bool, map<double,vector<Cluster*> > > > currentClusters;
	map<bool, map<bool, map<double,int> > > sizes;
	double max_z = numeric_limits<double>::min();
	double min_z = numeric_limits<double>::max();
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		double z;
		bool is_up;
		bool is_X;
		if((*it)->get_is_ref()){
			vector<Cluster*> temp_clusters((*it)->get_clusters());
			for(vector<Cluster*>::iterator jt=temp_clusters.begin();jt!=temp_clusters.end();++jt){
				z = (*jt)->get_z();
				is_up = (*jt)->get_is_up();
				is_X = (*jt)->get_is_X();
				if(z>max_z) max_z = z;
				if(z<min_z) min_z = z;
				currentClusters[is_up][is_X][z].push_back((*jt)->Clone());
				sizes[is_up][is_X][z]++;
				delete *jt;
			}
		}
	}
	map<bool, map<bool, vector<Ray_2D> > > suitableRays;
	//compute for both and down sensitive areas
	for(map<bool, map<bool, map<double,vector<Cluster*> > > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		//compute for both coordinates
		for(map<bool, map<double,vector<Cluster*> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			bool b = true;
			
			//get size
			//for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			//	if((kt->second).size()==0) b = false;
			//	sizes[it->first][jt->first][kt->first] = (kt->second).size();
			//}
			
			//find the biggest number of good clusters combinaisons
			while(b && (jt->second).size()>1){
				b = false;
				//find best combinaison of clusters
				vector<map<double,int> > comb = combinaisons(sizes[it->first][jt->first]);
				double current_chiSquare = chiSquare_threshold;
				map<double,int> best_comb;
				char coord = (jt->first) ? 'X' : 'Y';
				Ray_2D bestRay = Ray_2D(coord);
				for(vector<map<double,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
					//try a comb
					Ray_2D currentRay = Ray_2D(coord);
					for(map<double,vector<Cluster*> >::iterator nt = (jt->second).begin();nt!= (jt->second).end();++nt){
						currentRay.add_cluster(nt->second[(*kt)[nt->first]]);
					}
					currentRay.process();
					if(currentRay.get_chiSquare()<current_chiSquare && currentRay.get_chiSquare()>-1){
						bestRay = Ray_2D(currentRay);
						current_chiSquare = currentRay.get_chiSquare();
						best_comb = *kt;
						b = true;
					}
				}
				if(b){
					suitableRays[it->first][jt->first].push_back(Ray_2D(bestRay));
					for(map<double,int>::iterator kt = best_comb.begin();kt!=best_comb.end();++kt){
						delete (jt->second)[kt->first][kt->second];
						(jt->second)[kt->first].erase((jt->second)[kt->first].begin()+kt->second);
						sizes[it->first][jt->first][kt->first]--;
						if(sizes[it->first][jt->first][kt->first]<1){
							//b = false;
							(jt->second).erase((jt->second).find(kt->first));
							sizes[it->first][jt->first].erase(sizes[it->first][jt->first].find(kt->first));
						}
					}
				}
			}
		}
	}

	//cout << Max(Max(suitableRays[true][true].size(),suitableRays[true][false].size()),Max(suitableRays[false][true].size(),suitableRays[false][false].size())) << endl;
	int min_size = Min(Min(suitableRays[true][true].size(),suitableRays[true][false].size()),Min(suitableRays[false][true].size(),suitableRays[false][false].size()));
	
	while(min_size>0){
		double bestDoca = 50;//numeric_limits<double>::max();
		vector<Ray_2D>::iterator best_it = suitableRays[true][true].end();
		vector<Ray_2D>::iterator best_jt = suitableRays[true][false].end();
		vector<Ray_2D>::iterator best_kt = suitableRays[false][true].end();
		vector<Ray_2D>::iterator best_lt = suitableRays[false][false].end();
		for(vector<Ray_2D>::iterator it = suitableRays[true][true].begin(); it!= suitableRays[true][true].end();++it){
			for(vector<Ray_2D>::iterator jt = suitableRays[true][false].begin(); jt!= suitableRays[true][false].end();++jt){
				for(vector<Ray_2D>::iterator kt = suitableRays[false][true].begin(); kt!= suitableRays[false][true].end();++kt){
					for(vector<Ray_2D>::iterator lt = suitableRays[false][false].begin(); lt!= suitableRays[false][false].end();++lt){
						RayPair currentRayPair(Ray(*it,*jt),Ray(*kt,*lt));
						double currentDoca = currentRayPair.get_doca();
						Point currentPoCA = currentRayPair.get_PoCA();
						if(currentPoCA.get_Z()>max_z || currentPoCA.get_Z()<min_z) continue;
						if(currentPoCA.get_X()>6.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_X()<-6.*Tomography::get_instance()->get_XY_size()/10.) continue;
						if(currentPoCA.get_Y()>6.*Tomography::get_instance()->get_XY_size()/10. || currentPoCA.get_Y()<-6.*Tomography::get_instance()->get_XY_size()/10.) continue;
						
						if(currentDoca<bestDoca){
							bestDoca = currentDoca;
							best_it = vector<Ray_2D>::iterator(it);
							best_jt = vector<Ray_2D>::iterator(jt);
							best_kt = vector<Ray_2D>::iterator(kt);
							best_lt = vector<Ray_2D>::iterator(lt);
						}
					}
				}
			}
		}
		if(best_it == suitableRays[true][true].end() || best_jt == suitableRays[true][false].end() || best_kt == suitableRays[false][true].end() || best_lt == suitableRays[false][false].end()) break;
		rayPairs.push_back(RayPair(Ray(*best_it,*best_jt),Ray(*best_kt,*best_lt)));
		rayPairs.back().upRay.angle_correction();
		rayPairs.back().downRay.angle_correction();
		suitableRays[true][true].erase(best_it);
		suitableRays[true][false].erase(best_jt);
		suitableRays[false][true].erase(best_kt);
		suitableRays[false][false].erase(best_lt);
		min_size--;
	}

	for(map<bool, map<bool, map<double,vector<Cluster*> > > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		for(map<bool, map<double,vector<Cluster*> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				for(unsigned int i=0;i<(kt->second).size();i++){
					delete kt->second[i];
				}
			}
		}
	}
}
*/

vector<RayPair> CosmicBenchEvent::get_rayPairs() const{
	return rayPairs;
}

RayPair CosmicBenchEvent::get_rayPair(unsigned int i) const{
	return rayPairs[i];
}
unsigned int CosmicBenchEvent::get_rayPairs_N() const{
	return rayPairs.size();
}
unsigned int CosmicBenchEvent::get_event_N() const{
	return events.size();
}
unsigned int CosmicBenchEvent::get_clus_N() const{
	unsigned int clus_N = 0;
	for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
		clus_N += (*it)->get_NClus();
	}
	return clus_N;
}
unsigned int CosmicBenchEvent::get_clus_N_by_det(const Detector * const det) const{
	for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
		Detector * currentDet = (*it)->get_det();
		if((*currentDet) == (*det)){
			delete currentDet;
			return (*it)->get_NClus();
		}
		delete currentDet;
	}
	return -1;
}
vector<Ray> CosmicBenchEvent::get_absorption_rays(double chiSquare_threshold){
	this->Demux_CM();
	if(chiSquare_threshold<0) chiSquare_threshold = Tomography::get_instance()->get_chisquare_threshold();
	map<bool, map<int,vector<Cluster*> > > currentClusters;
	map<bool, map<int,int> > sizes;
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		if((*it)->get_is_ref()){
			vector<Cluster*> tempClusters = (*it)->get_clusters();
			for(vector<Cluster*>::iterator jt=tempClusters.begin();jt!=tempClusters.end();++jt){
				currentClusters[(*jt)->get_is_X()][(*jt)->get_layer()].push_back((*jt)->Clone());
				sizes[(*jt)->get_is_X()][(*jt)->get_layer()]++;
				delete *jt;
			}
		}
	}
	map<bool, vector<Ray_2D> > suitableRays;
	vector<Ray> returnRays;
	
	if(sizes[true].size()<2 || sizes[false].size()<2){
		for(map<bool, map<int,vector<Cluster*> > >::iterator jt = currentClusters.begin();jt!=currentClusters.end();++jt){
			for(map<int,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				for(unsigned int i=0;i<(kt->second).size();i++){
					delete ((kt->second)[i]);
				}
			}
		}
		return returnRays;
	}
	
	//compute for both acoordinates
	for(map<bool, map<int,vector<Cluster*> > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		if((it->second).size()<3){
			Ray_2D unique_ray = Ray_2D((it->first) ? 'X' : 'Y');
			for(map<int,vector<Cluster*> >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
				vector<Cluster*>::iterator best_clust = (jt->second).end();
				double max_ampl = 0;
				for(vector<Cluster*>::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
					if((*kt)->get_ampl()>max_ampl){
						best_clust = kt;
						max_ampl = (*kt)->get_ampl();
					}
				}
				unique_ray.add_cluster(*best_clust);
			}
			unique_ray.process();
			suitableRays[it->first].push_back(unique_ray);
		}
		else{
			bool b = true;
			/*
			//get size
			for(map<double,vector<Cluster*> >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
				if((jt->second).size()==0) b = false;
				sizes[it->first][jt->first] = (jt->second).size();
			}
			*/
			//find the biggest number of good clusters combinaisons
			//if no ray found, allow to drop a detector
			for(int drop = 0;drop<2;drop++){
				if(suitableRays[it->first].size() > 0) break;
				// you can adjust the size to require more or less hit
				while(b && (it->second).size()>2){
					b = false;
					//find best combinaison of clusters
					vector<map<int,int> > comb = combinaisons(sizes[it->first], (drop > 0));
					double current_chiSquare = numeric_limits<double>::max();
					map<int,int> best_comb;
					char coord = (it->first) ? 'X' : 'Y';
					Ray_2D bestRay = Ray_2D(coord);
					for(vector<map<int,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
						//try a comb
						Ray_2D currentRay = Ray_2D(coord);
						/*
						bool has_up = false;
						bool has_down = false;
						*/
						for(map<int,int>::iterator nt = kt->begin();nt!= kt->end();++nt){
							currentRay.add_cluster(it->second[nt->first][nt->second]);
							/*
							if(nt->second[(*kt)[nt->first]]->get_is_up()) has_up = true;
							else has_down = true;
							*/
						}
						//if(!(has_up && has_down)) continue;
						currentRay.process();
						/*double sigma = currentRay.get_t_sigma();*/
						if(currentRay.get_chiSquare()<current_chiSquare && currentRay.get_chiSquare()>-1 && (currentRay.get_chiSquare()/currentRay.get_clus_n())<(2.*chiSquare_threshold)){
							bestRay = currentRay;
							current_chiSquare = currentRay.get_chiSquare();//sigma;
							best_comb = *kt;
							b = true;
						}
					}
					if(b){
						suitableRays[it->first].push_back(Ray_2D(bestRay));
						for(map<int,int>::iterator kt = best_comb.begin();kt!=best_comb.end();++kt){
							delete ((it->second)[kt->first][kt->second]);
							(it->second)[kt->first].erase((it->second)[kt->first].begin()+kt->second);
							sizes[it->first][kt->first]--;
							if(sizes[it->first][kt->first]<1){
								//b = false;
								(it->second).erase((it->second).find(kt->first));
								sizes[it->first].erase(sizes[it->first].find(kt->first));
							}
						}
					}
				}
			}
		}
	}
	/*
	unsigned int min_size = (suitableRays.size()>1) ? numeric_limits<unsigned int>::max() : 0;
	for(map<bool, vector<Ray_2D> >::iterator it=suitableRays.begin();it!=suitableRays.end();++it){
		if((it->second).size()<min_size) min_size = (it->second).size();
	}
	for(unsigned int i = 0;i<min_size;i++){
		returnRays.push_back(Ray(suitableRays[true][i],suitableRays[false][i]));
		returnRays.back().angle_correction();
	}
	*/
	for(vector<Ray_2D>::iterator jt = suitableRays[true].begin(); (jt!=suitableRays[true].end()) && (!suitableRays[false].empty());++jt){
		vector<Cluster*> x_clus = jt->get_clus();
		map<int,pair<Tomography::det_type,int> > x_dets;
		for(vector<Cluster*>::iterator kt = x_clus.begin();kt!=x_clus.end();++kt){
			 x_dets.insert(pair<int,pair<Tomography::det_type,int> >((*kt)->get_layer(),pair<Tomography::det_type,int>((*kt)->get_type(),(*kt)->get_n_in_tree())));
			 delete *kt;
		}
		for(vector<Ray_2D>::iterator kt = suitableRays[false].begin(); kt!=suitableRays[false].end();++kt){
			vector<Cluster*> y_clus = kt->get_clus();
			bool mismatch = false;
			for(vector<Cluster*>::iterator nt = y_clus.begin();nt!=y_clus.end();++nt){
				Detector * y_det = detectors->get_detector(detectors->find_det(*nt));
				map<int,pair<Tomography::det_type,int> >::iterator in_layer = x_dets.find((*nt)->get_layer());
				if(in_layer != x_dets.end())
				{
					if(((in_layer->second).first != y_det->get_perp_type()) ||  ((in_layer->second).second != y_det->get_perp_n())) mismatch = true;
				}
				delete *nt;
			}
			if(!mismatch){
				returnRays.push_back(Ray(*jt,*kt));
				returnRays.back().angle_correction();
				suitableRays[false].erase(kt);
				break;
			}
		}
	}
	for(map<bool, map<int,vector<Cluster*> > >::iterator jt = currentClusters.begin();jt!=currentClusters.end();++jt){
		for(map<int,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			for(unsigned int i=0;i<(kt->second).size();i++){
				delete ((kt->second)[i]);
			}
		}
	}
	return returnRays;
}

vector<Ray> CosmicBenchEvent::get_hough_rays(double chiSquare_threshold){
	if(chiSquare_threshold<0) chiSquare_threshold = Tomography::get_instance()->get_chisquare_threshold();
	vector<Ray> returnRays;
}

void CosmicBenchEvent::EventDisplay(TCanvas * c1){
	gStyle->SetOptStat(false);
	double chisquare_threshold = 100;
	int detector_color = 1;
	int ray_color = 2;
	int pos_color = 3;
	vector<Ray> eventRays = this->get_absorption_rays(chisquare_threshold);
	vector<Ray>::iterator rays_it = eventRays.begin();
	while(rays_it!=eventRays.end()){
		if(((rays_it->get_chiSquare_X()+rays_it->get_chiSquare_Y())/rays_it->get_clus_n())>chisquare_threshold){
			eventRays.erase(rays_it);
			rays_it = eventRays.begin();
		}
		else{
			++rays_it;
		}
	}
	bool is_null = (c1==0);
	if(is_null) c1 = new TCanvas();
	TCanvas * cDisplay = c1;
	ostringstream title;
	title << "Event " << evn;
	cDisplay->SetTitle(title.str().c_str());
	cDisplay->Clear();
	cDisplay->Divide(2);
	map<double,vector<TH1D*> > ampl_hists_X;
	map<double,vector<TH1D*> > ampl_hists_Y;
	map<double,vector<double> > clus_pos_X;
	map<double,vector<double> > clus_pos_Y;
	//map<double,TCanvas*> cHist;
	double min_dist = 10000;
	vector<TLine*> det_X;
	vector<TLine*> det_Y;
	for(vector<Event*>::iterator event_it = events.begin();event_it!=events.end();++event_it){
		if((*event_it)->get_is_X()){
			for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
				if(Abs(map_it->first - (*event_it)->get_z())<min_dist) min_dist = Abs(map_it->first - (*event_it)->get_z());
			}
			ampl_hists_X[(*event_it)->get_z()].push_back((*event_it)->get_ampl_hist());
			/*
			if(!ampl_hists_X.insert(pair<double,TH1D*>((*event_it)->get_z(),(*event_it)->get_ampl_hist())).second){
				cout << "problem in events" << endl;
				return;
			}
			*/
			vector<Cluster*> current_clusters = (*event_it)->get_clusters();
			for(vector<Cluster*>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				//clus_pos_X[(*event_it)->get_z()].push_back((*clus_it)->get_pos()*(*event_it)->get_StripPitch() - Tomography::get_instance()->get_XY_size()/2.);
				clus_pos_X[(*event_it)->get_z()].push_back((*clus_it)->get_pos_mm());
				delete *clus_it;
			}
			det_X.push_back((*event_it)->get_det()->get_line_display());
			(det_X.back())->SetLineColor(detector_color);
		}
		else{
			for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
				if(Abs(map_it->first - (*event_it)->get_z())<min_dist) min_dist = Abs(map_it->first - (*event_it)->get_z());
			}
			ampl_hists_Y[(*event_it)->get_z()].push_back((*event_it)->get_ampl_hist());
			/*
			if(!ampl_hists_Y.insert(pair<double,TH1D*>((*event_it)->get_z(),(*event_it)->get_ampl_hist())).second){
				cout << "problem in events" << endl;
				return;
			}
			*/
			vector<Cluster*> current_clusters = (*event_it)->get_clusters();
			for(vector<Cluster*>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				//clus_pos_Y[(*event_it)->get_z()].push_back((*clus_it)->get_pos()*(*event_it)->get_StripPitch() - Tomography::get_instance()->get_XY_size()/2.);
				clus_pos_Y[(*event_it)->get_z()].push_back((*clus_it)->get_pos_mm());
				delete *clus_it;
			}
			det_Y.push_back((*event_it)->get_det()->get_line_display());
			(det_Y.back())->SetLineColor(detector_color);
		}
	}
	double scale_factor = 1;
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			//cout << map_it->first << " : " << (*hist_it)->GetMaximum() << endl;
			if((*hist_it)->GetMaximum()>scale_factor) scale_factor = (*hist_it)->GetMaximum();
		}
	}
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			//cout << map_it->first << " : " << (*hist_it)->GetMaximum() << endl;
			if((*hist_it)->GetMaximum()>scale_factor) scale_factor = (*hist_it)->GetMaximum();
		}
	}
	double scale = 1.1;
	scale_factor *= scale;
	scale_factor = min_dist/scale_factor;
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			//(*hist_it)->Scale(scale_factor);
			if((*hist_it)->GetMaximum() > 0) (*hist_it)->Scale(min_dist/(scale*(*hist_it)->GetMaximum()));
			TF1 * offset = new TF1("offset","[0]",-Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.);
			offset->SetParameter(0,map_it->first);
			(*hist_it)->Add(offset);
			delete offset;
			//cHist[map_it->first]->cd(1);
			//(*hist_it)->Draw();
		}
	}
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			//(*hist_it)->Scale(scale_factor);
			if((*hist_it)->GetMaximum() > 0) (*hist_it)->Scale(min_dist/(scale*(*hist_it)->GetMaximum()));
			TF1 * offset = new TF1("offset","[0]",-Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.);
			offset->SetParameter(0,map_it->first);
			(*hist_it)->Add(offset);
			delete offset;
			//cHist[map_it->first]->cd(1);
			//(*hist_it)->Draw();
		}
	}
	vector<TLine*> clus_X;
	vector<TLine*> clus_Y;
	double min_z = Min(ampl_hists_X.begin()->first,ampl_hists_Y.begin()->first);
	double max_z = Max((--(ampl_hists_X.end()))->first,(--(ampl_hists_Y.end()))->first);
	double diff_z = max_z-min_z;
	min_z -= 0.1*diff_z;
	max_z += 0.1*diff_z;
	max_z += min_dist;
	for(map<double,vector<double> >::iterator it = clus_pos_X.begin();it!=clus_pos_X.end();++it){
		for(vector<double>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			clus_X.push_back(new TLine(*jt,it->first,*jt,Min(it->first + (min_dist/scale),max_z)));
			(clus_X.back())->SetLineColor(pos_color);
			(clus_X.back())->SetLineStyle(2);
		}
	}
	for(map<double,vector<double> >::iterator it = clus_pos_Y.begin();it!=clus_pos_Y.end();++it){
		for(vector<double>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			clus_Y.push_back(new TLine(*jt,it->first,*jt,Min(it->first + (min_dist/scale),max_z)));
			(clus_Y.back())->SetLineColor(pos_color);
			(clus_Y.back())->SetLineStyle(2);
		}
	}
	vector<TLine*> rays_X;
	vector<TLine*> rays_Y;
	for(rays_it = eventRays.begin();rays_it!=eventRays.end();++rays_it){
		rays_X.push_back(new TLine(rays_it->eval_X(min_z),min_z,rays_it->eval_X(max_z),max_z));
		rays_Y.push_back(new TLine(rays_it->eval_Y(min_z),min_z,rays_it->eval_Y(max_z),max_z));
		(rays_X.back())->SetLineColor(ray_color);
		(rays_Y.back())->SetLineColor(ray_color);
	}
	/*
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		det_X.push_back(new TLine(-Tomography::get_instance()->get_XY_size()/2.,map_it->first,Tomography::get_instance()->get_XY_size()/2.,map_it->first));
		det_Y.push_back(new TLine(-Tomography::get_instance()->get_XY_size()/2.,map_it->first,Tomography::get_instance()->get_XY_size()/2.,map_it->first));
		(det_X.back())->SetLineColor(detector_color);
		(det_Y.back())->SetLineColor(detector_color);

	}
	*/
	TH1D * bg_X = new TH1D("XZ plane","XZ plane",2,-6*Tomography::get_instance()->get_XY_size()/10.,6*Tomography::get_instance()->get_XY_size()/10.);
	TH1D * bg_Y = new TH1D("YZ plane","YZ plane",2,-6*Tomography::get_instance()->get_XY_size()/10.,6*Tomography::get_instance()->get_XY_size()/10.);
	//bg_X->Fill(-Tomography::get_instance()->get_XY_size()/2.,min_z);
	//bg_X->Fill(Tomography::get_instance()->get_XY_size()/2.,max_z);
	bg_X->SetAxisRange(min_z,max_z,"Y");
	bg_X->GetXaxis()->SetTitle("X [mm]");
	bg_X->GetYaxis()->SetTitle("Z [mm]");
	bg_X->SetDirectory(0);
	//bg_Y->Fill(-Tomography::get_instance()->get_XY_size()/2.,min_z);
	//bg_Y->Fill(Tomography::get_instance()->get_XY_size()/2.,max_z);
	bg_Y->SetAxisRange(min_z,max_z,"Y");
	bg_Y->GetXaxis()->SetTitle("Y [mm]");
	bg_Y->GetYaxis()->SetTitle("Z [mm]");
	bg_Y->SetDirectory(0);
	cDisplay->cd(1);
	bg_X->Draw("AXIS");
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			(*hist_it)->SetDirectory(0);
			(*hist_it)->Draw("SAME ][");
		}
	}
	for(vector<TLine*>::iterator line_it = clus_X.begin();line_it!=clus_X.end();++line_it){
		(*line_it)->Draw();
	}
	for(vector<TLine*>::iterator line_it = rays_X.begin();line_it!=rays_X.end();++line_it){
		(*line_it)->Draw();
	}
	for(vector<TLine*>::iterator line_it = det_X.begin();line_it!=det_X.end();++line_it){
		(*line_it)->Draw();
	}
	cDisplay->cd(2);
	bg_Y->Draw("AXIS");
	for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
			(*hist_it)->SetDirectory(0);
			(*hist_it)->Draw("SAME ][");
		}
	}
	for(vector<TLine*>::iterator line_it = clus_Y.begin();line_it!=clus_Y.end();++line_it){
		(*line_it)->Draw();
	}
	for(vector<TLine*>::iterator line_it = rays_Y.begin();line_it!=rays_Y.end();++line_it){
		(*line_it)->Draw();
	}
	for(vector<TLine*>::iterator line_it = det_Y.begin();line_it!=det_Y.end();++line_it){
		(*line_it)->Draw();
	}
	TText * canvas_title = new TText(0.5,0.965,title.str().c_str());
	canvas_title->SetTextAlign(22);
	cDisplay->cd();
	canvas_title->Draw();
	TText * left_pad_title = new TText(0.5,0.96,bg_X->GetTitle());
	left_pad_title->SetNDC(true);
	left_pad_title->SetTextAlign(22);
	cDisplay->cd(1);
	left_pad_title->Draw();
	TText * right_pad_title = new TText(0.5,0.96,bg_Y->GetTitle());
	right_pad_title->SetNDC(true);
	right_pad_title->SetTextAlign(22);
	cDisplay->cd(2);
	right_pad_title->Draw();
	cDisplay->GetPad(1)->SetFillStyle(0);
	cDisplay->GetPad(2)->SetFillStyle(0);
	cDisplay->Modified();
	cDisplay->Update();
	if(!is_null){
		for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
			for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
				delete (*hist_it);
			}
		}
		for(map<double,vector<TH1D*> >::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
			for(vector<TH1D*>::iterator hist_it = (map_it->second).begin();hist_it != (map_it->second).end();++hist_it){
				delete (*hist_it);
			}
		}
		delete bg_X; delete bg_Y; delete canvas_title; delete left_pad_title; delete right_pad_title;
		for(vector<TLine*>::iterator line_it = clus_X.begin();line_it!=clus_X.end();++line_it){
			delete (*line_it);
		}
		for(vector<TLine*>::iterator line_it = rays_X.begin();line_it!=rays_X.end();++line_it){
			delete (*line_it);
		}
		for(vector<TLine*>::iterator line_it = det_X.begin();line_it!=det_X.end();++line_it){
			delete (*line_it);
		}
		for(vector<TLine*>::iterator line_it = clus_Y.begin();line_it!=clus_Y.end();++line_it){
			delete (*line_it);
		}
		for(vector<TLine*>::iterator line_it = rays_Y.begin();line_it!=rays_Y.end();++line_it){
			delete (*line_it);
		}
		for(vector<TLine*>::iterator line_it = det_Y.begin();line_it!=det_Y.end();++line_it){
			delete (*line_it);
		}
	}
}

void CosmicBenchEvent::do_cuts(){
	Demux_CM();
	for(vector<Event*>::iterator it = events.begin(); it!=events.end(); ++it){
		(*it)->do_cuts();
	}
}
int CosmicBenchEvent::get_evn() const{
	return (events.front())->get_evn();
}
double CosmicBenchEvent::get_evttime() const{
	return (events.front())->get_evttime();
}
