#include "tomography.h"
#include "datareader.h"
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <TROOT.h>
#include <TCanvas.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

using std::map;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::string;

using boost::property_tree::ptree;

bool Tomography::can_continue = true;
bool Tomography::is_batch = gROOT->IsBatch();
bool Tomography::live_graphic_display = !Tomography::is_batch;

ostream& Tomography::operator<<(ostream& os, const det_type& det){
	switch(det){
		case CM : os << "CM"; break;
		case CM_Demux : os << "CM_Demux"; break;
		case MG : os << "MG"; break;
		default : os << "unknown det";
	}
	return os;
}
ostream& Tomography::operator<<(ostream& os, const strip_type& strip){
	switch(strip){
		case Wide : os << "Wide"; break;
		case Thin : os << "Thin"; break;
		case Demux : os << "Demux"; break;
		default : os << "unknown strip type";
	}
	return os;
}
ostream& Tomography::operator<<(ostream& os, const elec_type& elec){
	switch(elec){
		case Dream : os << "Dream"; break;
		case Feminos : os << "Feminos"; break;
		default : os << "unknown electronic type";
	}
	return os;
}

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_){
	if(map_.size()<1){
		os << "[]";
		return os;
	}
	os << "[ ";
	typename map<T,R>::const_iterator it=map_.begin();
	typename map<T,R>::const_iterator jt = map_.end();
	jt--;
	while(it!=jt){
		os << it->first << " -> " << it->second << " ; ";
		++it;
	}
	os << it->first << " -> " << it->second << " ]";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_){
	if(vec_.empty()){
		os << "{}";
		return os;
	}
	os << "{ ";
	typename vector<T>::const_iterator it = vec_.begin();
	typename vector<T>::const_iterator jt = vec_.end();
	jt--;
	while(it!=jt){
		os << *it << " ; ";
		++it;
	}
	os << *it << " }";
	return os;
}

//template ostream& operator<<(ostream& os, const map<int,int>& map_);
template ostream& operator<<(ostream& os, const map<int,double>& map_);
template ostream& operator<<(ostream& os, const map<double,int>& map_);
template ostream& operator<<(ostream& os, const map<int,int>& map_);
template ostream& operator<<(ostream& os, const map<bool,map<int,int> >& map_);

void Tomography::signal_handler(int s){
	cout << "\nCaught signal " << s << endl;
	cout << endl;
	can_continue = false;
}

Tomography::elec_type Tomography::str_to_elec(string str){
	elec_type return_value = unknown_elec;
	if(str == "dream") return_value = Dream;
	else if(str == "feminos") return_value = Feminos;
	return return_value;
}
static map<Tomography::det_type,int> CMN_div_build(){
	map<Tomography::det_type,int> return_map;
	return_map[Tomography::CM] = 2;
	return_map[Tomography::MG] = 2;
	return return_map;
}
const map<Tomography::det_type,int> Tomography::CMN_div = CMN_div_build();
/*
void Tomography::process_elec_files(ptree config_tree){
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	int CM_N = 0;
	int MG_N = 0;
	map<int,Tomography::det_type> det_type_by_asic;
	map<int,int> det_n_by_asic;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::CM;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("cm_n");
		CM_N++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::MG;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("mg_n");
		MG_N++;
	}
	if((total_CM_N!=CM_N) || (total_MG_N!=MG_N)){
		cout << "problem in detectors number" << endl;
		return;
	}
	Tomography::elec_type electronic_type = Tomography::str_to_elec(config_tree.get<string>("electronic_type"));
	string data_file_basename = config_tree.get<string>("data_file_basename");
	string signalName = config_tree.get<string>("signal_file");
	string PedName = config_tree.get<string>("Ped");
	string RMSName = config_tree.get<string>("RMSPed");
	int max_event = config_tree.get<int>("max_event");
	ifstream file;
	file.open(signalName.c_str(),ifstream::in);
	bool exists = file.good();
	file.close();
	file.open(PedName.c_str(),ifstream::in);
	bool ped_done = file.good();
	file.close();
	file.open(RMSName.c_str(),ifstream::in);
	bool compute_rms = !(file.good());
	file.close();
	DataReader * blah = NULL;
	if(electronic_type == Tomography::Dream){
		blah = new DreamDataReader(signalName,PedName,RMSName,det_type_by_asic,det_n_by_asic,exists,exists,exists,max_event);
	}
	else if(electronic_type == Tomography::Feminos){
		blah = new FeminosDataReader(signalName,PedName,RMSName,det_type_by_asic,det_n_by_asic,exists,exists,exists,max_event);
	}
	else{
		cout << "electronic type : " << electronic_type << " unknown !" << endl;
		return;
	}
	int first_data_file = config_tree.get<int>("data_file_first");
	int last_data_file = config_tree.get<int>("data_file_last") + 1;
	for(int i=first_data_file;i<last_data_file;i++){
		ostringstream dataFileName;
		dataFileName << data_file_basename << setw(3) << setfill('0') << i;
		if(electronic_type == Tomography::Dream) dataFileName << "_01.fdf";
		else if(electronic_type == Tomography::Feminos) dataFileName << ".aqs";
		else dataFileName << ".txt";
		blah->add_file_to_process(dataFileName.str());
	}
	blah->process();
	if(!ped_done) blah->compute_ped();
	blah->do_ped_sub();
	blah->do_common_noise_sub();
	if(compute_rms) blah->compute_RMSPed();
}
*/
void Tomography::save_canvases(){
	time_t current_time = time(NULL);
	char buffer[100];
	strftime(buffer,100,"%y%m%d_%HH%M",localtime(&current_time));
	string base_name = "canvas_";
	base_name += buffer;
	base_name += "_";
	for(int i=0;i<gROOT->GetListOfCanvases()->GetSize();i++){
		TCanvas * current_canvas = dynamic_cast<TCanvas*>(gROOT->GetListOfCanvases()->At(i));
		current_canvas->SaveAs((base_name + current_canvas->GetName() + ".png").c_str());
		current_canvas->SaveAs((base_name + current_canvas->GetName() + ".C").c_str());
	}

}