#define event_cpp
#include "event.h"
#include <string>
#include <map>
#include <vector>
#include <limits>
#include "T.h"
#include "detector.h"
#include "cluster.h"
#include "ray.h"
#include <iostream>
#include <TMath.h>

using std::string;
using std::map;
using std::vector;
using std::numeric_limits;
using std::cout;
using std::endl;

using TMath::Min;
using TMath::Max;

vector<map<double,int> > CosmicBenchEvent::combinaisons(map<double,int> sizes){
	map<double,int> partial_product;
	int current_product = 1;
	for(map<double,int>::iterator it = sizes.begin();it!=sizes.end();++it){
		partial_product[it->first] = current_product;
		current_product*=it->second;
	}
	vector<map<double,int> > result(current_product,map<double,int>());
	for(int i=0;i<current_product;i++){
		for(map<double,int>::iterator it = sizes.begin();it!=sizes.end();++it){
			result[i][it->first] = (i/partial_product[it->first]) % it->second;
		}
	}
	return result;
}


Event::Event(){
	type ="";
	n_in_tree = -1;
	has_spark = true;
	is_ref = false;
	z = -1;
}
Event::Event(const Event& other){
	evn = other.evn;
	type = other.type;
	n_in_tree = other.n_in_tree;
	has_spark = other.has_spark;
	is_ref = other.is_ref;
	z = other.z;
}
Event& Event::operator=(const Event& other){
	evn = other.evn;
	type = other.type;
	n_in_tree = other.n_in_tree;
	has_spark = other.has_spark;
	is_ref = other.is_ref;
	z = other.z;
	return *this;
}
Event::Event(T * treeObject,int entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	has_spark = true;
	is_ref = false;
	z = -1;
}
Event::~Event(){

}
int Event::get_evn() const{
	return evn;
}
string Event::get_type() const{
	return type;
}
int Event::get_n_in_tree() const{
	return n_in_tree;
}
bool Event::get_is_ref() const{
	return is_ref;
}
double Event::get_z() const{
	return z;
}

CM_Event::CM_Event(): Event(){
	type = "CM";
}
CM_Event::CM_Event(const CM_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "CM";
}
CM_Event& CM_Event::operator=(const CM_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "CM";
	return *this;
}
CM_Event::CM_Event(T * treeObject,CM_Detector * det,int entry): Event(treeObject,entry){
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
	has_spark = (treeObject->CM_Spark[n_in_tree]==1) ? true : false;
	is_ref = det->get_is_ref();
	type = "CM";
}
CM_Event::~CM_Event(){
	clusters.clear();
}

CM_Demux_Event::CM_Demux_Event(): Event(){
	type = "CM_Demux";
}
CM_Demux_Event::CM_Demux_Event(const CM_Demux_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "CM_Demux";
}
CM_Demux_Event& CM_Demux_Event::operator=(const CM_Demux_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "CM_Demux";
	return *this;
}
CM_Demux_Event::CM_Demux_Event(const CM_Event& rawEvent){
	clusters.clear();
	if(rawEvent.use_thin_strip){
		for(vector<CM_Cluster>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(it->get_strip_type() == "Wide"){
				for(vector<CM_Cluster>::const_iterator jt = rawEvent.clusters.begin();jt!=rawEvent.clusters.end();++jt){
					if(jt->get_strip_type() == "Thin"){
						clusters.push_back(CM_Demux_Cluster(*jt,*it));
					}
				}
			}
		}
	}
	else{
		for(vector<CM_Cluster>::const_iterator it = rawEvent.clusters.begin();it!=rawEvent.clusters.end();++it){
			if(it->get_strip_type() == "Wide"){
				clusters.push_back(CM_Demux_Cluster(*it));
			}
		}
	}
	n_in_tree = rawEvent.n_in_tree;
	has_spark = rawEvent.has_spark;
	is_ref = rawEvent.is_ref;
	z = rawEvent.z;
	type = "CM_Demux";
}
vector<CM_Demux_Cluster> CM_Demux_Event::get_clusters() const{
	return clusters;
}
CM_Demux_Event::~CM_Demux_Event(){
	clusters.clear();
}

MG_Event::MG_Event(): Event(){
	type = "MG";
}
MG_Event::MG_Event(const MG_Event& other): Event(other){
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "MG";
}
MG_Event& MG_Event::operator=(const MG_Event& other){
	Event::operator=(other);
	clusters.clear();
	clusters.assign(other.clusters.begin(),other.clusters.end());
	type = "MG";
	return *this;
}
MG_Event::MG_Event(T * treeObject,MG_Detector * det, int entry): Event(treeObject,entry){
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
	type = "MG";
}
vector<MG_Cluster> MG_Event::get_clusters() const{
	return clusters;
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
		if((*it)->type == "CM") events.push_back(new CM_Event(*(dynamic_cast<CM_Event*>(*it))));
		else if((*it)->type == "CM_Demux") events.push_back(new CM_Demux_Event(*(dynamic_cast<CM_Demux_Event*>(*it))));
		else if((*it)->type == "MG") events.push_back(new MG_Event(*(dynamic_cast<MG_Event*>(*it))));
	}
}
CosmicBenchEvent& CosmicBenchEvent::operator=(const CosmicBenchEvent& other){
	evn = other.evn;
	rayPairs.clear();
	rayPairs.assign(other.rayPairs.begin(),other.rayPairs.end());
	events.clear();
	for(vector<Event*>::const_iterator it = other.events.begin();it!= other.events.end();++it){
		if((*it)->type == "CM") events.push_back(new CM_Event(*(dynamic_cast<CM_Event*>(*it))));
		else if((*it)->type == "CM_Demux") events.push_back(new CM_Demux_Event(*(dynamic_cast<CM_Demux_Event*>(*it))));
		else if((*it)->type == "MG") events.push_back(new MG_Event(*(dynamic_cast<MG_Event*>(*it))));
	}
	return *this;
}
CosmicBenchEvent::CosmicBenchEvent(CosmicBench * detectors, T * treeObject, int entry){
	if(entry>-1){
		treeObject->LoadTree(entry);
		treeObject->GetEntry(entry);
	}
	evn = treeObject->evn;
	rayPairs.clear();
	events.clear();
	int det_N = detectors->get_CM_N() + detectors->get_MG_N();
	for(int i=0;i<det_N;i++){
		if((detectors->get_detector(i))->get_type() == "CM"){
			events.push_back(new CM_Event(treeObject,(dynamic_cast<CM_Detector*>(detectors->get_detector(i))),-1));
		}
		else if((detectors->get_detector(i))->get_type() == "MG"){
			events.push_back(new MG_Event(treeObject,(dynamic_cast<MG_Detector*>(detectors->get_detector(i))),-1));
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
		if((*it)->get_type() == "CM"){
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
			string det_type = (*it)->get_type();
			if(det_type == "CM_Demux"){
				CM_Demux_Event * currentEvent = dynamic_cast<CM_Demux_Event*>(*it);
				for(vector<CM_Demux_Cluster>::iterator jt=(currentEvent->clusters).begin();jt!=(currentEvent->clusters).end();++jt){
					z = jt->z;
					is_up = jt->is_up;
					is_X = jt->is_X;
					if(z>max_z) max_z = z;
					if(z<min_z) min_z = z;
					currentClusters[is_up][is_X][z].push_back(new CM_Demux_Cluster(*jt));
				}
			}
			else if(det_type == "MG"){
				MG_Event * currentEvent = dynamic_cast<MG_Event*>(*it);
				for(vector<MG_Cluster>::iterator jt=(currentEvent->clusters).begin();jt!=(currentEvent->clusters).end();++jt){
					z = jt->z;
					is_up = jt->is_up;
					is_X = jt->is_X;
					if(z>max_z) max_z = z;
					if(z<min_z) min_z = z;
					currentClusters[is_up][is_X][z].push_back(new MG_Cluster(*jt));
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
			//get size
			for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				if((kt->second).size()==0) b = false;
				sizes[it->first][jt->first][kt->first] = (kt->second).size();
			}
			//find the biggest number of good clusters combinaisons
			while(b){
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
						if(sizes[it->first][jt->first][kt->first]<1) b = false;
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
		if((*it)->type == "CM"){
			clus_N += ((dynamic_cast<CM_Event*>(*it))->clusters).size();
		}
		else if((*it)->type == "CM_Demux"){
			clus_N += ((dynamic_cast<CM_Demux_Event*>(*it))->clusters).size();
		}
		else if((*it)->type == "MG"){
			clus_N += ((dynamic_cast<MG_Event*>(*it))->clusters).size();
		}
	}
	return clus_N;
}
unsigned int CosmicBenchEvent::get_clus_N_by_det(Detector * det) const{
	unsigned int clus_N = 0;
	if(det->get_type() == "CM"){
		for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
			if((*it)->type == "CM"){
				if((*it)->n_in_tree == (dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree()){
					clus_N += ((dynamic_cast<CM_Event*>(*it))->clusters).size();
				}
			}
			else if((*it)->type == "CM_Demux"){
				if((*it)->n_in_tree == (dynamic_cast<CM_Detector*>(det))->get_cm_n_in_tree()){
					clus_N += ((dynamic_cast<CM_Demux_Event*>(*it))->clusters).size();
				}
			}
		}
	}
	else if(det->get_type() == "MG"){
		for(vector<Event*>::const_iterator it=events.begin();it!=events.end();++it){
			if((*it)->type == "MG"){
				if((*it)->n_in_tree == (dynamic_cast<MG_Detector*>(det))->get_mg_n_in_tree()){
					clus_N += ((dynamic_cast<MG_Event*>(*it))->clusters).size();
				}
			}
		}
	}
	return clus_N;
}
vector<Ray> CosmicBenchEvent::get_absorption_rays(){
	this->Demux_CM();
	double chiSquare_threshold = numeric_limits<double>::max();
	map<bool, map<double,vector<Cluster*> > > currentClusters;
	map<bool, map<double,int> > sizes;
	for(vector<Event*>::iterator it=events.begin();it!=events.end();++it){
		if((*it)->get_is_ref()){
			string det_type = (*it)->get_type();
			if(det_type == "CM_Demux"){
				CM_Demux_Event * currentEvent = dynamic_cast<CM_Demux_Event*>(*it);
				for(vector<CM_Demux_Cluster>::iterator jt=(currentEvent->clusters).begin();jt!=(currentEvent->clusters).end();++jt){
					currentClusters[jt->is_X][jt->z].push_back(new CM_Demux_Cluster(*jt));
				}
			}
			else if(det_type == "MG"){
				MG_Event * currentEvent = dynamic_cast<MG_Event*>(*it);
				for(vector<MG_Cluster>::iterator jt=(currentEvent->clusters).begin();jt!=(currentEvent->clusters).end();++jt){
					currentClusters[jt->is_X][jt->z].push_back(new MG_Cluster(*jt));
				}
			}
			else continue;
		}
	}
	map<bool, vector<Ray_2D> > suitableRays;
	//compute for both acoordinates
	for(map<bool, map<double,vector<Cluster*> > >::iterator it = currentClusters.begin();it!=currentClusters.end();++it){
		bool b = true;
		//get size
		for(map<double,vector<Cluster*> >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			if((jt->second).size()==0) b = false;
			sizes[it->first][jt->first] = (jt->second).size();
		}
		//find the biggest number of good clusters combinaisons
		while(b && (it->second).size()>1){
			b = false;
			//find best combinaison of clusters
			vector<map<double,int> > comb = combinaisons(sizes[it->first]);
			double current_chiSquare = chiSquare_threshold;
			map<double,int> best_comb;
			char coord = (it->first) ? 'X' : 'Y';
			Ray_2D bestRay = Ray_2D(coord);
			for(vector<map<double,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
				//try a comb
				Ray_2D currentRay = Ray_2D(coord);
				for(map<double,vector<Cluster*> >::iterator nt = (it->second).begin();nt!= (it->second).end();++nt){
					currentRay.add_cluster(nt->second[(*kt)[nt->first]]);
				}
				currentRay.process();
				/*double sigma = currentRay.get_t_sigma();*/
				if(currentRay.get_chiSquare()<current_chiSquare/*sigma<current_chiSquare*/ && currentRay.get_chiSquare()>-1){
					bestRay = currentRay;
					current_chiSquare = currentRay.get_chiSquare();//sigma;
					best_comb = *kt;
					b = true;
				}
			}
			if(b){
				suitableRays[it->first].push_back(Ray_2D(bestRay));
				for(map<double,int>::iterator kt = best_comb.begin();kt!=best_comb.end();++kt){
					delete (it->second)[kt->first][kt->second];
					(it->second)[kt->first].erase((it->second)[kt->first].begin()+kt->second);
					sizes[it->first][kt->first]--;
					if(sizes[it->first][kt->first]<1) b = false;
				}
			}
		}
	}
	unsigned int min_size = (suitableRays.size()>1) ? numeric_limits<unsigned int>::max() : 0;
	for(map<bool, vector<Ray_2D> >::iterator it=suitableRays.begin();it!=suitableRays.end();++it){
		if((it->second).size()<min_size) min_size = (it->second).size();
	}
	vector<Ray> returnRays;
	for(unsigned int i = 0;i<min_size;i++){
		returnRays.push_back(Ray(suitableRays[true][i],suitableRays[false][i]));
	}
	for(map<bool, map<double,vector<Cluster*> > >::iterator jt = currentClusters.begin();jt!=currentClusters.end();++jt){
		for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			for(unsigned int i=0;i<(kt->second).size();i++){
				delete kt->second[i];
			}
		}
	}
	return returnRays;
}