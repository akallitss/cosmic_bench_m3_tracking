#define event_cpp
#include "event.h"

#include "T.h"
#include "detector.h"
#include "cluster.h"
#include "ray.h"

#include <string>
#include <map>
#include <vector>
#include <limits>
#include <iostream>
#include <utility>
#include <sstream>
#include <algorithm>

#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>

using std::pair;
using std::string;
using std::map;
using std::vector;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::flush;
using std::ostringstream;
using std::max_element;

using TMath::Min;
using TMath::Max;
using TMath::Abs;
using TMath::FloorNint;
using TMath::CeilNint;
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

Event::Event(){
	n_in_tree = -1;
	has_spark = true;
	is_ref = false;
	z = -1;
	NClus = -1;
	use_srf = false;
	is_X = false;
}
Event::Event(const Event& other){
	evn = other.evn;
	type = other.type;
	n_in_tree = other.n_in_tree;
	has_spark = other.has_spark;
	is_ref = other.is_ref;
	z = other.z;
	NClus = other.NClus;
	use_srf = other.use_srf;
	strip_ampl = other.strip_ampl;
	is_X = other.is_X;
}
Event& Event::operator=(const Event& other){
	evn = other.evn;
	type = other.type;
	n_in_tree = other.n_in_tree;
	has_spark = other.has_spark;
	is_ref = other.is_ref;
	z = other.z;
	NClus = other.NClus;
	use_srf = other.use_srf;
	strip_ampl = other.strip_ampl;
	is_X = other.is_X;
	return *this;
}
Event::Event(T * treeObject, bool use_srf_,long entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	has_spark = true;
	is_ref = false;
	z = -1;
	NClus = -1;
	use_srf = use_srf_;
	is_X = false;
}
Event::~Event(){

}
int Event::get_evn() const{
	return evn;
}
Tomography::det_type Event::get_type() const{
	return type;
}
int Event::get_n_in_tree() const{
	return n_in_tree;
}
bool Event::get_is_ref() const{
	return is_ref;
}
bool Event::get_is_X() const{
	return is_X;
}
double Event::get_z() const{
	return z;
}
int Event::get_NClus() const{
	return NClus;
}

CM_Event::CM_Event(): Event(){
	type = Tomography::CM;
	detector = CM_Detector();
}
CM_Event::CM_Event(const CM_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::CM;
	detector = other.detector;
}
CM_Event& CM_Event::operator=(const CM_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::CM;
	detector = other.detector;
	return *this;
}
CM_Event::CM_Event(T * treeObject,CM_Detector * det, bool use_srf_,long entry): Event(treeObject,use_srf_,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	clusters.clear();
	for(int i=0;i<treeObject->CM_NClus[det->get_cm_n_in_tree()];i++){
		if(CM_Cluster::is_suitable(treeObject,i,det,-1)){
			clusters.push_back(CM_Cluster(treeObject,i,det,-1));
		}
	}
	n_in_tree = det->get_cm_n_in_tree();
	use_thin_strip = det->get_use_thin_strip();
	z = det->get_z();
	is_X = det->get_is_X();
	has_spark = (treeObject->CM_Spark[n_in_tree]==1) ? true : false;
	is_ref = det->get_is_ref();
	detector = *det;
	type = Tomography::CM;
}
CM_Event::CM_Event(CM_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_, int evn_): Event(){
	if(strip_ampl_.size()!=64){
		cout << "problem in size" << endl;
		detector = CM_Detector();
		return;
	}
	n_in_tree = detector_.get_cm_n_in_tree();
	is_ref = detector_.get_is_ref();
	z = detector_.get_z();
	is_X = detector_.get_is_X();
	strip_ampl = strip_ampl_;
	detector = detector_;
	use_srf = use_srf_;
	evn = evn_;
	type = Tomography::CM;
}
void CM_Event::MultiCluster(){
	// TODO : implement multicluster for CM
}
void CM_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=64){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
vector<CM_Cluster> CM_Event::get_clusters() const{
	return clusters;
}
void CM_Event::do_cuts(){
	vector<CM_Cluster>::iterator clus_it = clusters.begin();
	while(clus_it!= clusters.end()){
		if(clus_it->is_suitable(&detector)) ++clus_it;
		else clus_it = clusters.erase(clus_it);
	}
}
CM_Event::~CM_Event(){
	clusters.clear();
}

CM_Demux_Event::CM_Demux_Event(): Event(){
	type = Tomography::CM_Demux;
}
CM_Demux_Event::CM_Demux_Event(const CM_Demux_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::CM_Demux;
}
CM_Demux_Event& CM_Demux_Event::operator=(const CM_Demux_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::CM_Demux;
	return *this;
}
CM_Demux_Event::CM_Demux_Event(const CM_Event& rawEvent){
	clusters.clear();
	if(rawEvent.use_thin_strip){
		for(vector<CM_Cluster>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(it->get_strip_type() == Tomography::Wide){
				for(vector<CM_Cluster>::const_iterator jt = rawEvent.clusters.begin();jt!=rawEvent.clusters.end();++jt){
					if(jt->get_strip_type() == Tomography::Thin){
						clusters.push_back(CM_Demux_Cluster(*jt,*it));
					}
				}
			}
		}
	}
	else{
		for(vector<CM_Cluster>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(it->get_strip_type() == Tomography::Wide){
				clusters.push_back(CM_Demux_Cluster(*it));
			}
		}
	}
	n_in_tree = rawEvent.n_in_tree;
	has_spark = rawEvent.has_spark;
	is_ref = rawEvent.is_ref;
	z = rawEvent.z;
	is_X = rawEvent.is_X;
	strip_ampl = rawEvent.strip_ampl;
	type = Tomography::CM_Demux;
}
vector<CM_Demux_Cluster> CM_Demux_Event::get_clusters() const{
	return clusters;
}
void CM_Demux_Event::MultiCluster(){
	// TODO : implement multicluster for CM
}
void CM_Demux_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=64){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
void CM_Demux_Event::do_cuts(){
	cout << "you shouldn't call this function" << endl;
}
CM_Demux_Event::~CM_Demux_Event(){
	clusters.clear();
}

MG_Event::MG_Event(): Event(){
	type = Tomography::MG;
	detector = MG_Detector();
}
MG_Event::MG_Event(const MG_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::MG;
	detector = other.detector;
}
MG_Event& MG_Event::operator=(const MG_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = Tomography::MG;
	detector = other.detector;
	return *this;
}
MG_Event::MG_Event(T * treeObject,MG_Detector * det, bool use_srf_,long entry): Event(treeObject,use_srf_,entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	clusters.clear();
	for(int i=0;i<treeObject->MG_NClus[det->get_mg_n_in_tree()];i++){
		if(MG_Cluster::is_suitable(treeObject,i,det,-1)) clusters.push_back(MG_Cluster(treeObject,i,det,-1));
	}
	n_in_tree = det->get_mg_n_in_tree();
	has_spark = (treeObject->MG_Spark[n_in_tree]==1) ? true : false;
	is_ref = det->get_is_ref();
	z = det->get_z();
	is_X = det->get_is_X();
	detector = *det;
	type = Tomography::MG;
}
MG_Event::MG_Event(MG_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_, int evn_): Event(){
	if(strip_ampl_.size()!=61){
		cout << "problem in size" << endl;
		detector = MG_Detector();
		return;
	}
	n_in_tree = detector_.get_mg_n_in_tree();
	is_ref = detector_.get_is_ref();
	z = detector_.get_z();
	is_X = detector_.get_is_X();
	strip_ampl = strip_ampl_;
	detector = detector_;
	use_srf = use_srf_;
	evn = evn_;
	type = Tomography::MG;
}
vector<MG_Cluster> MG_Event::get_clusters() const{
	return clusters;
}

void MG_Event::MultiCluster(){
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::sigma;
	int SampleMin = Tomography::SampleMin;
	int SampleMax = Tomography::SampleMax;
	int TOTCut = Tomography::TOTCut;
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
		for(int j=SampleMax;j<Tomography::Nsample;j++){
			current_strip.signal_sample[j] = false;
		}
		for(int j=SampleMin;j<SampleMax;j++){
			if(strip_ampl[i][j]>current_strip.MaxAmpl){
				current_strip.MaxAmpl = strip_ampl[i][j];
				current_strip.MaxSample = j;
				// time calculation with maximum
				
				if(j>0 && j<31){
					double a = (0.5*strip_ampl[i][j+1]) - strip_ampl[i][j] + (0.5*strip_ampl[i][j-1]);
					//double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					double b = (0.5*(strip_ampl[i][j+1] - strip_ampl[i][j-1])) - (2*a*j);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
				
				// --
			}
			if(strip_ampl[i][j]>(sigma*(detector.get_RMS(i)))){
				current_strip.TOT++;
				current_strip.signal_sample[j] = true;
			}
			else current_strip.signal_sample[j] = false;
		}
		// time calculation with rising edge
		/*
		if(current_strip.TOT>2){
			int k=current_strip.MaxSample;
			TGraph * rising_edge = new TGraph();
			while(k>=SampleMin && strip_ampl[i][k]>(sigma*(detector.get_RMS(i)))){
				rising_edge->SetPoint(current_strip.MaxSample - k,k,strip_ampl[i][k]);
				k--;
			}
			TF1 * rising_fit = new TF1("rising_fit","pol1(0)",k-1,current_strip.MaxSample +1);
			rising_fit->SetParameters(k,Tomography::ADC_max/(2.*Tomography::Nsample));
			rising_fit->SetParLimits(0,0,Tomography::Nsample);
			rising_fit->SetParLimits(1,0,Tomography::ADC_max);
			rising_edge->Fit(rising_fit,"QN");
			current_strip.Time = rising_fit->GetParameter(0);
			delete rising_edge; delete rising_fit;
		}
		else current_strip.Time = 0;
		*/
		// --
		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	vector<int> global_used_channel;
	unsigned int max_hole_size = detector.get_clustering_holes();
	while(global_used_channel.size()<channelOverThreshold.size()){
		pair<int,int> biggest_current_cluster(0,-1);
		vector<int> current_used_channel;
		for(int i=0;i<n;i++){
			if(channelOverThreshold.count(MG_Detector::StripToChannel[i])>0 && find(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel[i]) == global_used_channel.end()){
				pair<int,int> current_cluster(i,i);
				vector<int> used_channel;
				used_channel.push_back(MG_Detector::StripToChannel[i]);
				map<int,int> hole_channel;
				for(int j=i+1;j<n;j++){
					if(/*find(used_channel.begin(),used_channel.end(),MG_Detector::StripToChannel[j]) == used_channel.end() &&*/ find(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel[j]) == global_used_channel.end() /*&& !(hole_channel.count(MG_Detector::StripToChannel[j])>0)*/){
						if(channelOverThreshold.count(MG_Detector::StripToChannel[j])>0){
							current_cluster.second = j;
							used_channel.push_back(MG_Detector::StripToChannel[j]);
						}
						else if(hole_channel.size()<max_hole_size){
							hole_channel.insert(pair<int,int>(MG_Detector::StripToChannel[j],j));
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
	NClus = cluster_list.size();
	MG_Detector * detector_p = &detector;
	for(unsigned int i=0;i<cluster_list.size();i++){
		TF1 * SRFfit = new TF1("SRFfit",detector_p,&MG_Detector::SRF_fit,0,1024,2,"MG_Detector","SRF_fit");
		TGraphErrors * SRFgraph = new TGraphErrors();
		int graph_point_n = 0;
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		double tot_ampl = 0;
		double pos_TPC = 0;
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MG_Detector::StripToChannel[j]];
			double effective_ampl = current_strip.MaxAmpl/count(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel[j]);
			ClusPos = (ClusPos*ClusAmpl + j*effective_ampl)/(ClusAmpl + effective_ampl);
			ClusAmpl += effective_ampl;

			//Micro TPC
			/*
			for(int k=0;k<Tomography::Nsample;k++){
				if(!(current_strip.signal_sample[k])) continue;
				double current_tot_ampl = strip_ampl[MG_Detector::StripToChannel[j]][k]/count(global_used_channel.begin(),global_used_channel.end(),MG_Detector::StripToChannel[j]);
				pos_TPC = (pos_TPC*tot_ampl + j*current_tot_ampl)/(tot_ampl + current_tot_ampl);
				tot_ampl += current_tot_ampl;
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
			SRFgraph->SetPoint(graph_point_n,j,effective_ampl);
			SRFgraph->SetPointError(graph_point_n,0.5*MG_Detector::StripPitch,detector.get_RMS(MG_Detector::StripToChannel[j]));
			graph_point_n++;
		}

		if(graph_point_n>2 && use_srf){
			for(int iPoint=0;iPoint<graph_point_n;iPoint++){
				double x_current = -1;
				double y_current = -1;
				SRFgraph->GetPoint(iPoint,x_current,y_current);
				y_current /= ClusMaxStripAmpl;
				//y_current = strip_ampl[MG_Detector::StripToChannel[static_cast<int>(x_current)]][static_cast<int>(ClusMaxSample)]/ClusMaxStripAmpl;
				SRFgraph->SetPoint(iPoint,x_current,y_current);
				SRFgraph->SetPointError(iPoint,MG_Detector::StripPitch,(SRFgraph->GetEY())[iPoint]/ClusMaxStripAmpl);
			}
			SRFfit->SetParameter(0,ClusPos);
			SRFfit->SetParLimits(0,ClusPos-0.5*ClusSize,ClusPos+0.5*ClusSize);
			SRFfit->SetParameter(1,0.25*ClusSize);
			SRFfit->SetParLimits(1,0,ClusSize);
			SRFgraph->Fit(SRFfit,"QN");
			/*
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
			*/
			ClusPos = SRFfit->GetParameter(0);
			//ClusSize = SRFfit->GetParameter(1);
		}
		delete SRFfit; delete SRFgraph;
		clusters.push_back(MG_Cluster(&detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
		//clusters.push_back(MG_Cluster(&detector,i,pos_TPC,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MG_Event::HoughCluster(){
	clusters.clear();
	//first loop : find channels with signal and store them with their caracteristics
	double sigma = Tomography::sigma;
	int SampleMin = Tomography::SampleMin;
	int SampleMax = Tomography::SampleMax;
	int TOTCut = Tomography::TOTCut;
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
				if(j>0 && j<31){
					double a = 0.5*strip_ampl[i][j+1] - 2*strip_ampl[i][j] + strip_ampl[i][j-1];
					double b = strip_ampl[i][j] - strip_ampl[i][j-1] - a*((2*j)-1);
					current_strip.Time = -0.5*b/a;
				}
				else current_strip.Time = 0;
			}
			if(strip_ampl[i][j]>sigma*(detector.get_RMS(i))) current_strip.TOT++;
		}
		if(current_strip.TOT>TOTCut) channelOverThreshold.insert(pair<int,bool>(i,true));
		allChannels.insert(pair<int,StripInfo>(i,current_strip));
	}
	//second loop : group the channels in clusters
	vector<pair<int,int> > cluster_list;
	int k = 0;
	while(k<n){
		if(channelOverThreshold.count(MG_Detector::StripToChannel[k])>0){
			pair<int,int> current_cluster(k,k);
			for(int j=k+1;j<n;j++){
				if(channelOverThreshold.count(MG_Detector::StripToChannel[j])>0){
					current_cluster.second = j;
				}
				else break;
			}
			cluster_list.push_back(current_cluster);
			k = current_cluster.second+1;
		}
		else k++;
	}
	//third loop : store the clusters and their caracteristics
	NClus = cluster_list.size();
	MG_Detector * detector_p = &detector;
	for(unsigned int i=0;i<cluster_list.size();i++){
		TF1 * SRFfit = new TF1("SRFfit",detector_p,&MG_Detector::SRF_fit,0,1024,2,"MG_Detector","SRF_fit");
		TGraphErrors * SRFgraph = new TGraphErrors();
		int graph_point_n = 0;
		double ClusSize = 1 + cluster_list[i].second - cluster_list[i].first;
		double ClusPos = 0;
		double ClusAmpl = 0;
		double ClusMaxStripAmpl = 0;
		double ClusMaxSample = 0;
		int ClusMaxStrip = -1;
		double ClusT = 0;
		double ClusTOT = 0;
		for(int j = cluster_list[i].first;j<((cluster_list[i].second)+1);j++){
			StripInfo current_strip = allChannels[MG_Detector::StripToChannel[j]];
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
			SRFgraph->SetPoint(graph_point_n,j,effective_ampl);
			SRFgraph->SetPointError(graph_point_n,0.5*MG_Detector::StripPitch,detector.get_RMS(MG_Detector::StripToChannel[j]));
			graph_point_n++;
		}

		if(graph_point_n>2 && use_srf){
			for(int iPoint=0;iPoint<graph_point_n;iPoint++){
				double x_current = -1;
				double y_current = -1;
				SRFgraph->GetPoint(iPoint,x_current,y_current);
				y_current /= ClusMaxStripAmpl;
				//y_current = strip_ampl[MG_Detector::StripToChannel[static_cast<int>(x_current)]][static_cast<int>(ClusMaxSample)]/ClusMaxStripAmpl;
				SRFgraph->SetPoint(iPoint,x_current,y_current);
				SRFgraph->SetPointError(iPoint,MG_Detector::StripPitch,(SRFgraph->GetEY())[iPoint]/ClusMaxStripAmpl);
			}
			SRFfit->SetParameter(0,ClusPos);
			SRFfit->SetParLimits(0,ClusPos-0.5*ClusSize,ClusPos+0.5*ClusSize);
			SRFfit->SetParameter(1,0.25*ClusSize);
			SRFfit->SetParLimits(1,0,ClusSize);
			SRFgraph->Fit(SRFfit,"QN");
			/*
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
			*/
			ClusPos = SRFfit->GetParameter(0);
			ClusSize = SRFfit->GetParameter(1);
		}
		delete SRFfit; delete SRFgraph;
		clusters.push_back(MG_Cluster(&detector,i,ClusPos,ClusSize,ClusAmpl,ClusMaxSample,ClusMaxStripAmpl,ClusTOT,ClusT,ClusMaxStrip));
	}
}
void MG_Event::set_strip_ampl(vector<vector<double> > strip_ampl_){
	if(strip_ampl_.size()!=61){
		cout << "problem in size" << endl;
		return;
	}
	strip_ampl = strip_ampl_;
}
TH1D * MG_Event::get_ampl_hist() const{
	ostringstream name;
	name << "ampl_MG_det_" << n_in_tree << "_evn_" << evn;
	TH1D * histo = new TH1D(name.str().c_str(),name.str().c_str(),1024,-MG_Detector::size/2.,MG_Detector::size/2.);
	vector<pair<int,int> > cluster_edges;
	for(vector<MG_Cluster>::const_iterator it=clusters.begin();it!=clusters.end();++it){
		cluster_edges.push_back(pair<int,int>(FloorNint(it->get_pos()-it->get_size()),CeilNint(it->get_pos()+it->get_size())));
	}
	vector<bool> is_used(1024,false);
	for(vector<pair<int,int> >::iterator it = cluster_edges.begin();it!=cluster_edges.end();++it){
		if(it->first<0) it->first = 0;
		if(it->second<0) it->second = 0;
		if(it->first>1023) it->first = 1023;
		if(it->second>1023) it->second = 1023;
		for(int strip=it->first;strip<=it->second;strip++){
			if(is_used[strip]) continue;
			int channel = MG_Detector::StripToChannel[strip];
			histo->Fill(strip*MG_Detector::StripPitch - MG_Detector::size/2.,*max_element(strip_ampl[channel].begin(),strip_ampl[channel].end()));
			is_used[strip] = true;
		}
	}
	return histo;
}
void MG_Event::do_cuts(){
	vector<MG_Cluster>::iterator clus_it = clusters.begin();
	while(clus_it!= clusters.end()){
		if(clus_it->is_suitable(&detector)) ++clus_it;
		else clus_it = clusters.erase(clus_it);
	}
}
MG_Event::~MG_Event(){
	clusters.clear();
}

CosmicBenchEvent::CosmicBenchEvent(){
	rayPairs.clear();
	events.clear();
}
CosmicBenchEvent::CosmicBenchEvent(const CosmicBenchEvent& other){
	evn = other.evn;
	rayPairs.clear();
	rayPairs.assign(other.rayPairs.begin(),other.rayPairs.end());
	events.clear();
	for(vector<Event*>::const_iterator it = other.events.begin();it!= other.events.end();++it){
		if((*it)->get_type() == Tomography::CM) events.push_back(new CM_Event(*(dynamic_cast<CM_Event*>(*it))));
		else if((*it)->get_type() == Tomography::CM_Demux) events.push_back(new CM_Demux_Event(*(dynamic_cast<CM_Demux_Event*>(*it))));
		else if((*it)->get_type() == Tomography::MG) events.push_back(new MG_Event(*(dynamic_cast<MG_Event*>(*it))));
	}
}
CosmicBenchEvent& CosmicBenchEvent::operator=(const CosmicBenchEvent& other){
	evn = other.evn;
	rayPairs.clear();
	rayPairs.assign(other.rayPairs.begin(),other.rayPairs.end());
	events.clear();
	for(vector<Event*>::const_iterator it = other.events.begin();it!= other.events.end();++it){
		if((*it)->get_type() == Tomography::CM) events.push_back(new CM_Event(*(dynamic_cast<CM_Event*>(*it))));
		else if((*it)->get_type() == Tomography::CM_Demux) events.push_back(new CM_Demux_Event(*(dynamic_cast<CM_Demux_Event*>(*it))));
		else if((*it)->get_type() == Tomography::MG) events.push_back(new MG_Event(*(dynamic_cast<MG_Event*>(*it))));
	}
	return *this;
}
CosmicBenchEvent::CosmicBenchEvent(CosmicBench * detectors, T * treeObject, bool use_srf_, long entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	rayPairs.clear();
	events.clear();
	int det_N = detectors->get_CM_N() + detectors->get_MG_N();
	for(int i=0;i<det_N;i++){
		if((detectors->get_detector(i))->get_type() == Tomography::CM){
			events.push_back(new CM_Event(treeObject,(dynamic_cast<CM_Detector*>(detectors->get_detector(i))),use_srf_,-1));
		}
		else if((detectors->get_detector(i))->get_type() == Tomography::MG){
			events.push_back(new MG_Event(treeObject,(dynamic_cast<MG_Detector*>(detectors->get_detector(i))),use_srf_,-1));
		}
	}
}
CosmicBenchEvent::CosmicBenchEvent(CosmicBench * detectors, vector<Event*> events_){
	rayPairs.clear();
	events.clear();
	unsigned int det_N = detectors->get_CM_N() + detectors->get_MG_N();
	if(events_.size()!=det_N){
		cout << "problem in event number" << endl;
		return;
	}
	map<int,bool> CM_is_used;
	map<int,bool> MG_is_used;
	for(unsigned int i=0;i<det_N;i++){
		Detector * current_det = detectors->get_detector(i);
		if(current_det->get_type() == Tomography::MG){
			MG_is_used[dynamic_cast<MG_Detector*>(current_det)->get_mg_n_in_tree()] = false;
		}
		else if(current_det->get_type() == Tomography::CM){
			CM_is_used[dynamic_cast<CM_Detector*>(current_det)->get_cm_n_in_tree()] = false;
		}
	}
	evn = (events_.front())->get_evn();
	for(vector<Event*>::iterator it = events_.begin();it!=events_.end();++it){
		if((*it)->get_evn() != evn){
			cout << "attempt to merge event with different number in cosmicbench event" << endl;
			return;
		}
		if((*it)->get_type() == Tomography::MG){
			if(MG_is_used[(*it)->get_n_in_tree()]){
				cout << "problem in events" << endl;
				return;
			}
			else{
				MG_is_used[(*it)->get_n_in_tree()] = true;
				events.push_back(new MG_Event(*dynamic_cast<MG_Event*>(*it)));
			}
		}
		else if((*it)->get_type() == Tomography::CM){
			if(CM_is_used[(*it)->get_n_in_tree()]){
				cout << "problem in events" << endl;
				return;
			}
			else{
				CM_is_used[(*it)->get_n_in_tree()] = true;
				events.push_back(new CM_Event(*dynamic_cast<CM_Event*>(*it)));
			}
		}
		else if((*it)->get_type() == Tomography::CM_Demux){
			if(CM_is_used[(*it)->get_n_in_tree()]){
				cout << "problem in events" << endl;
				return;
			}
			else{
				CM_is_used[(*it)->get_n_in_tree()] = true;
				events.push_back(new CM_Demux_Event(*dynamic_cast<CM_Demux_Event*>(*it)));
			}
		}
	}
}
CosmicBenchEvent::~CosmicBenchEvent(){
	for(unsigned int i=0;i<events.size();i++){
		delete events[i];
	}
	events.clear();
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
	map<bool, map<bool, map<double,vector<Cluster*> > > > currentClusters;
	map<bool, map<bool, map<double,int> > > sizes;
	double max_z = numeric_limits<double>::min();
	double min_z = numeric_limits<double>::max();
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		double z;
		bool is_up;
		bool is_X;
		if((*it)->get_is_ref()){
			Tomography::det_type det_type = (*it)->get_type();
			if(det_type == Tomography::CM_Demux){
				//CM_Demux_Event * currentEvent = dynamic_cast<CM_Demux_Event*>(*it);
				vector<CM_Demux_Cluster> currentCluster = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
				for(vector<CM_Demux_Cluster>::iterator jt=currentCluster.begin();jt!=currentCluster.end();++jt){
					z = jt->get_z();
					is_up = jt->get_is_up();
					is_X = jt->get_is_X();
					if(z>max_z) max_z = z;
					if(z<min_z) min_z = z;
					currentClusters[is_up][is_X][z].push_back(new CM_Demux_Cluster(*jt));
					sizes[is_up][is_X][z]++;
				}
			}
			else if(det_type == Tomography::MG){
				//MG_Event * currentEvent = dynamic_cast<MG_Event*>(*it);
				vector<MG_Cluster> currentCluster = (dynamic_cast<MG_Event*>(*it))->get_clusters();
				for(vector<MG_Cluster>::iterator jt=currentCluster.begin();jt!=currentCluster.end();++jt){
					z = jt->get_z();
					is_up = jt->get_is_up();
					is_X = jt->get_is_X();
					if(z>max_z) max_z = z;
					if(z<min_z) min_z = z;
					currentClusters[is_up][is_X][z].push_back(new MG_Cluster(*jt));
					sizes[is_up][is_X][z]++;
				}
			}
			else{
				continue;
			}
		}
	}
	map<bool, map<bool, vector<Ray_2D> > > suitableRays;
	//compute for both and down sensitive areas
	for(map<bool, map<bool, map<double,vector<Cluster*> > > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		//compute for both coordinates
		for(map<bool, map<double,vector<Cluster*> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			bool b = true;
			/*
			//get size
			for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				if((kt->second).size()==0) b = false;
				sizes[it->first][jt->first][kt->first] = (kt->second).size();
			}
			*/
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
						if(currentPoCA.get_X()>600 || currentPoCA.get_X()<-100) continue;
						if(currentPoCA.get_Y()>600 || currentPoCA.get_Y()<-100) continue;
						
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
		if((*it)->get_type() == Tomography::CM){
			clus_N += ((dynamic_cast<CM_Event*>(*it))->get_clusters()).size();
		}
		else if((*it)->get_type() == Tomography::CM_Demux){
			clus_N += ((dynamic_cast<CM_Demux_Event*>(*it))->get_clusters()).size();
		}
		else if((*it)->get_type() == Tomography::MG){
			clus_N += ((dynamic_cast<MG_Event*>(*it))->get_clusters()).size();
		}
	}
	return clus_N;
}
unsigned int CosmicBenchEvent::get_clus_N_by_det(Detector * det) const{
	unsigned int clus_N = 0;
	if(det->get_type() == Tomography::CM){
		for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
			if((*it)->get_type() == Tomography::CM){
				if((*it)->get_n_in_tree() == (dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree()){
					clus_N += ((dynamic_cast<CM_Event*>(*it))->get_clusters()).size();
				}
			}
			else if((*it)->get_type() == Tomography::CM_Demux){
				if((*it)->get_n_in_tree() == (dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree()){
					clus_N += ((dynamic_cast<CM_Demux_Event*>(*it))->get_clusters()).size();
				}
			}
		}
	}
	else if(det->get_type() == Tomography::MG){
		for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
			if((*it)->get_type() == Tomography::MG){
				if((*it)->get_n_in_tree() == (dynamic_cast<MG_Detector*>(det))->get_mg_n_in_tree()){
					clus_N += ((dynamic_cast<MG_Event*>(*it))->get_clusters()).size();
				}
			}
		}
	}
	return clus_N;
}
vector<Ray> CosmicBenchEvent::get_absorption_rays(double chiSquare_threshold){
	this->Demux_CM();
	map<bool, map<double,vector<Cluster*> > > currentClusters;
	map<bool, map<double,int> > sizes;
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		if((*it)->get_is_ref()){
			Tomography::det_type det_type = (*it)->get_type();
			if(det_type == Tomography::CM_Demux){
				//CM_Demux_Event * currentEvent = dynamic_cast<CM_Demux_Event*>(*it);
				vector<CM_Demux_Cluster> currentCluster = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
				for(vector<CM_Demux_Cluster>::iterator jt=currentCluster.begin();jt!=currentCluster.end();++jt){
					currentClusters[jt->get_is_X()][jt->get_z()].push_back(new CM_Demux_Cluster(*jt));
					sizes[jt->get_is_X()][jt->get_z()]++;
				}
			}
			else if(det_type == Tomography::MG){
				//MG_Event * currentEvent = dynamic_cast<MG_Event*>(*it);
				vector<MG_Cluster> currentCluster = (dynamic_cast<MG_Event*>(*it))->get_clusters();
				for(vector<MG_Cluster>::iterator jt=currentCluster.begin();jt!=currentCluster.end();++jt){
					currentClusters[jt->get_is_X()][jt->get_z()].push_back(new MG_Cluster(*jt));
					sizes[jt->get_is_X()][jt->get_z()]++;
				}
			}
			else continue;
		}
	}
	map<bool, vector<Ray_2D> > suitableRays;
	vector<Ray> returnRays;
	if(sizes[true].size()<3 || sizes[false].size()<3){
		for(map<bool, map<double,vector<Cluster*> > >::iterator jt = currentClusters.begin();jt!=currentClusters.end();++jt){
			for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				for(unsigned int i=0;i<(kt->second).size();i++){
					delete ((kt->second)[i]);
				}
			}
		}
		return returnRays;
	}
	//compute for both acoordinates
	for(map<bool, map<double,vector<Cluster*> > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
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
			if(suitableRays[it->first].size() > 0) continue;
			// you can adjust the size to require more or less hit
			while(b && (it->second).size()>2){
				b = false;
				//find best combinaison of clusters
				vector<map<double,int> > comb = combinaisons(sizes[it->first], (drop>0));
				double current_chiSquare = numeric_limits<double>::max();
				map<double,int> best_comb;
				char coord = (it->first) ? 'X' : 'Y';
				Ray_2D bestRay = Ray_2D(coord);
				for(vector<map<double,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
					//try a comb
					Ray_2D currentRay = Ray_2D(coord);
					/*
					bool has_up = false;
					bool has_down = false;
					*/
					for(map<double,int>::iterator nt = kt->begin();nt!= kt->end();++nt){
						currentRay.add_cluster(it->second[nt->first][nt->second]);
						/*
						if(nt->second[(*kt)[nt->first]]->get_is_up()) has_up = true;
						else has_down = true;
						*/
					}
					//if(!(has_up && has_down)) continue;
					currentRay.process();
					/*double sigma = currentRay.get_t_sigma();*/
					if(currentRay.get_chiSquare()<current_chiSquare && currentRay.get_chiSquare()>-1 && (currentRay.get_chiSquare()/currentRay.get_clus_n())<(10.*chiSquare_threshold)){
						bestRay = currentRay;
						current_chiSquare = currentRay.get_chiSquare();//sigma;
						best_comb = *kt;
						b = true;
					}
				}
				if(b){
					suitableRays[it->first].push_back(Ray_2D(bestRay));
					for(map<double,int>::iterator kt = best_comb.begin();kt!=best_comb.end();++kt){
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
	unsigned int min_size = (suitableRays.size()>1) ? numeric_limits<unsigned int>::max() : 0;
	for(map<bool, vector<Ray_2D> >::iterator it=suitableRays.begin();it!=suitableRays.end();++it){
		if((it->second).size()<min_size) min_size = (it->second).size();
	}
	for(unsigned int i = 0;i<min_size;i++){
		returnRays.push_back(Ray(suitableRays[true][i],suitableRays[false][i]));
		returnRays.back().angle_correction();
	}
	for(map<bool, map<double,vector<Cluster*> > >::iterator jt = currentClusters.begin();jt!=currentClusters.end();++jt){
		for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			for(unsigned int i=0;i<(kt->second).size();i++){
				delete ((kt->second)[i]);
			}
		}
	}
	return returnRays;
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
	cDisplay->Clear();
	cDisplay->Divide(2);
	map<double,TH1D*> ampl_hists_X;
	map<double,TH1D*> ampl_hists_Y;
	map<double,vector<double> > clus_pos_X;
	map<double,vector<double> > clus_pos_Y;
	//map<double,TCanvas*> cHist;
	double min_dist = 10000;
	for(vector<Event*>::iterator event_it = events.begin();event_it!=events.end();++event_it){
		if((*event_it)->get_type() != Tomography::MG){
			cout << "no display for CM, sorry :(" << endl;
				return;
		}
		if((*event_it)->get_is_X()){
			for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
				if(Abs(map_it->first - (*event_it)->get_z())<min_dist) min_dist = Abs(map_it->first - (*event_it)->get_z());
			}
			if(!ampl_hists_X.insert(pair<double,TH1D*>((*event_it)->get_z(),dynamic_cast<MG_Event*>(*event_it)->get_ampl_hist())).second){
				cout << "problem in events" << endl;
				return;
			}
			vector<MG_Cluster> current_clusters = dynamic_cast<MG_Event*>(*event_it)->get_clusters();
			for(vector<MG_Cluster>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				clus_pos_X[(*event_it)->get_z()].push_back(clus_it->get_pos()*MG_Detector::StripPitch - Tomography::XY_size/2.);
			}
		}
		else{
			for(map<double,TH1D*>::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
				if(Abs(map_it->first - (*event_it)->get_z())<min_dist) min_dist = Abs(map_it->first - (*event_it)->get_z());
			}
			if(!ampl_hists_Y.insert(pair<double,TH1D*>((*event_it)->get_z(),dynamic_cast<MG_Event*>(*event_it)->get_ampl_hist())).second){
				cout << "problem in events" << endl;
				return;
			}
			vector<MG_Cluster> current_clusters = dynamic_cast<MG_Event*>(*event_it)->get_clusters();
			for(vector<MG_Cluster>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				clus_pos_Y[(*event_it)->get_z()].push_back(clus_it->get_pos()*MG_Detector::StripPitch - Tomography::XY_size/2.);
			}
		}
	}
	double scale_factor = 1;
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		//cout << map_it->first << " : " << map_it->second->GetMaximum() << endl;
		if(map_it->second->GetMaximum()>scale_factor) scale_factor = map_it->second->GetMaximum();
		//cHist[map_it->first] = new TCanvas();
		//cHist[map_it->first]->Divide(2);
	}
	for(map<double,TH1D*>::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		//cout << map_it->first << " : " << map_it->second->GetMaximum() << endl;
		if(map_it->second->GetMaximum()>scale_factor) scale_factor = map_it->second->GetMaximum();
	}
	double scale = 1.1;
	scale_factor *= scale;
	scale_factor = min_dist/scale_factor;
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		//map_it->second->Scale(scale_factor);
		if(map_it->second->GetMaximum() > 0) map_it->second->Scale(min_dist/(scale*map_it->second->GetMaximum()));
		TF1 * offset = new TF1("offset","[0]",-Tomography::XY_size/2.,Tomography::XY_size/2.);
		offset->SetParameter(0,map_it->first);
		map_it->second->Add(offset);
		delete offset;
		//cHist[map_it->first]->cd(1);
		//map_it->second->Draw();
	}
	for(map<double,TH1D*>::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		//map_it->second->Scale(scale_factor);
		if(map_it->second->GetMaximum() > 0) map_it->second->Scale(min_dist/(scale*map_it->second->GetMaximum()));
		TF1 * offset = new TF1("offset","[0]",-Tomography::XY_size/2.,Tomography::XY_size/2.);
		offset->SetParameter(0,map_it->first);
		map_it->second->Add(offset);
		delete offset;
		//cHist[map_it->first]->cd(2);
		//map_it->second->Draw();
		//cHist[map_it->first]->Modified();
		//cHist[map_it->first]->Update();
	}
	vector<TLine*> clus_X;
	vector<TLine*> clus_Y;
	for(map<double,vector<double> >::iterator it = clus_pos_X.begin();it!=clus_pos_X.end();++it){
		for(vector<double>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			clus_X.push_back(new TLine(*jt,it->first,*jt,it->first + (min_dist/scale)));
			(clus_X.back())->SetLineColor(pos_color);
			(clus_X.back())->SetLineStyle(2);
		}
	}
	for(map<double,vector<double> >::iterator it = clus_pos_Y.begin();it!=clus_pos_Y.end();++it){
		for(vector<double>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			clus_Y.push_back(new TLine(*jt,it->first,*jt,it->first + (min_dist/scale)));
			(clus_Y.back())->SetLineColor(pos_color);
			(clus_Y.back())->SetLineStyle(2);
		}
	}
	double min_z = Min(ampl_hists_X.begin()->first,ampl_hists_Y.begin()->first);
	double max_z = Max((--(ampl_hists_X.end()))->first,(--(ampl_hists_Y.end()))->first);
	double temp = min_z;
	min_z -= 0.1*(max_z-min_z);
	max_z += 0.1*(max_z-temp);
	vector<TLine*> rays_X;
	vector<TLine*> rays_Y;
	for(rays_it = eventRays.begin();rays_it!=eventRays.end();++rays_it){
		rays_X.push_back(new TLine(rays_it->eval_X(min_z),min_z,rays_it->eval_X(max_z),max_z));
		rays_Y.push_back(new TLine(rays_it->eval_Y(min_z),min_z,rays_it->eval_Y(max_z),max_z));
		(rays_X.back())->SetLineColor(ray_color);
		(rays_Y.back())->SetLineColor(ray_color);
	}
	vector<TLine*> det_X;
	vector<TLine*> det_Y;
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		det_X.push_back(new TLine(-Tomography::XY_size/2.,map_it->first,Tomography::XY_size/2.,map_it->first));
		det_Y.push_back(new TLine(-Tomography::XY_size/2.,map_it->first,Tomography::XY_size/2.,map_it->first));
		(det_X.back())->SetLineColor(detector_color);
		(det_Y.back())->SetLineColor(detector_color);

	}
	TH1D * bg_X = new TH1D("XZ plane","XZ plane",2,-6*Tomography::XY_size/10.,6*Tomography::XY_size/10.);
	TH1D * bg_Y = new TH1D("YZ plane","YZ plane",2,-6*Tomography::XY_size/10.,6*Tomography::XY_size/10.);
	bg_X->Fill(-Tomography::XY_size/2.,min_z);
	bg_X->Fill(Tomography::XY_size/2.,max_z);
	bg_X->SetAxisRange(min_z,max_z,"Y");
	bg_Y->Fill(-Tomography::XY_size/2.,min_z);
	bg_Y->Fill(Tomography::XY_size/2.,max_z);
	bg_Y->SetAxisRange(min_z,max_z,"Y");
	cDisplay->cd(1);
	bg_X->Draw("AXIS");
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		map_it->second->DrawCopy("SAME ][");
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
	for(map<double,TH1D*>::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		map_it->second->DrawCopy("SAME ][");
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
	cDisplay->Modified();
	cDisplay->Update();
	for(map<double,TH1D*>::iterator map_it = ampl_hists_X.begin();map_it!=ampl_hists_X.end();++map_it){
		delete map_it->second;
	}
	for(map<double,TH1D*>::iterator map_it = ampl_hists_Y.begin();map_it!=ampl_hists_Y.end();++map_it){
		delete map_it->second;
	}
	if(!is_null){
		delete bg_X; delete bg_Y;
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