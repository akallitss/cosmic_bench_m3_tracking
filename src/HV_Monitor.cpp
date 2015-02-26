#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <sys/timex.h>
#include <sstream>
#include <fstream>
#include <cstdio>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "CAEN_comm.h"

#include <TFile.h>
#include <TTree.h>

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::map;
using std::vector;
using std::setw;
using std::setprecision;
using std::left;
using std::showpos;
using std::noshowpos;
using std::find;
using std::ostringstream;
using std::ofstream;

using boost::property_tree::ptree;

int main(int argc, char ** argv){

	if(argc<2){
		cout << "Indicate a config file name" << endl;
		return 1;
	}
	ptree config_tree;
	read_json(argv[1], config_tree);
	vector<string> channel_names;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("Channels")){
		channel_names.push_back(child.second.data());
	}
	CAEN_Comm * blah = new CAEN_Comm(config_tree.get<string>("IP").c_str(),config_tree.get<string>("username").c_str(),config_tree.get<string>("password").c_str());
	map<int,map<int,string> > chombier = blah->get_Ch_name();
	map<int,vector<int> > used_channel;
	int total_channel = 0;
	for(map<int,map<int,string> >::iterator it = chombier.begin();it!=chombier.end();++it){
		for(map<int,string>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			if(find(channel_names.begin(),channel_names.end(),jt->second) != channel_names.end()){
				used_channel[it->first].push_back(jt->first);
				total_channel++;
			}
		}
	}
	char ** chan_names = new char*[total_channel];
	int current_channel = 0;
	for(map<int,vector<int> >::iterator map_it=used_channel.begin();map_it!=used_channel.end();++map_it){
		for(vector<int>::iterator vec_it=(map_it->second).begin();vec_it!=(map_it->second).end();++vec_it){
			chan_names[current_channel] = const_cast<char*>(chombier[map_it->first][*vec_it].c_str());
			current_channel++;
		}
	}

	map<CAEN_Ch::Param,map<int,map<int,CAEN_Ch::param_value> > > channel_values;
	struct timespec wait_time;
	struct timespec remaining_time;
	wait_time.tv_sec = (config_tree.get<int>("tick_length")>0 ) ? (config_tree.get<int>("tick_length") - 1) : 0;
	struct ntptimeval current_time;
	ntp_gettime(&current_time);
	cout << "\rt = " << current_time.time.tv_sec << " | " << setw(4) << 0 << "%" << flush;

	TFile * fOut = new TFile(config_tree.get<string>("filename").c_str(),"RECREATE");
	TTree * outTree = new TTree("T","HV");
	outTree->Branch("time",&current_time.time.tv_sec,"time/L");
	ostringstream branch_name;
	branch_name << "name[" << total_channel << "]/C";
	outTree->Branch("name",chan_names,branch_name.str().c_str());

	map<CAEN_Ch::Param,CAEN_Ch::param_value*> channel_params;
	bool has_IMon = false;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("Params")){
		CAEN_Ch::Param current_param = CAEN_Ch::get_Param(child.second.data());
		if(channel_params.count(current_param)>0) continue;
		channel_params[current_param] = new CAEN_Ch::param_value[total_channel];
		ostringstream current_branch_name;
		current_branch_name << "HV_" << CAEN_Ch::get_Param_name(current_param);
		ostringstream current_branch_name_root;
		current_branch_name_root << current_branch_name.str() << "[" << total_channel << "]/";
		if(CAEN_Ch::Param_is_float(current_param)){
			current_branch_name_root << "F";
		}
		else{
			current_branch_name_root << "I";
		}
		outTree->Branch(current_branch_name.str().c_str(),channel_params[current_param],current_branch_name_root.str().c_str());
		if(current_param == CAEN_Ch::IMon) has_IMon = true;
	}

	ofstream IMon_csv;
	if(has_IMon){
		remove("tmp_IMon.csv");
		IMon_csv.open("tmp_IMon.csv");
		IMon_csv << "time;";
		for(int i=0;i<total_channel;i++){
			IMon_csv << chan_names[i] << ";";
		}
		IMon_csv << endl;
	}

	unsigned int duration = config_tree.get<int>("duration");
	for(unsigned int i=0;i<duration;i++){
		ntp_gettime(&current_time);
		wait_time.tv_nsec = 1000*(1000000-current_time.time.tv_usec);
		nanosleep(&wait_time,&remaining_time);
		for(map<CAEN_Ch::Param,CAEN_Ch::param_value*>::iterator params_it = channel_params.begin();params_it!=channel_params.end();++params_it){
			channel_values[params_it->first] = blah->get_Ch_param(used_channel,params_it->first);

		}
		if(has_IMon) IMon_csv << current_time.time.tv_sec << ";";
		for(map<CAEN_Ch::Param,map<int,map<int,CAEN_Ch::param_value> > >::iterator values_it = channel_values.begin();values_it!=channel_values.end();++values_it){
			current_channel = 0;
			for(map<int,vector<int> >::iterator map_it = used_channel.begin();map_it!=used_channel.end();++map_it){
				for(vector<int>::iterator vec_it = (map_it->second).begin();vec_it!=(map_it->second).end();++vec_it){
					CAEN_Ch::param_value current_value = (values_it->second)[map_it->first][*vec_it];
					channel_params[values_it->first][current_channel] = current_value;
					current_channel++;
					if(values_it->first == CAEN_Ch::IMon) IMon_csv << current_value.real << ";";
				}
			}
		}
		if(has_IMon) IMon_csv << endl;
		outTree->Fill();
		cout << "\rt = " << current_time.time.tv_sec << " | " << setw(6) << setprecision(4) << 100.*(i+1)/duration << "%" << flush;
	}
	cout << endl;
	outTree->Write();
	fOut->Close();
	delete blah;
	if(has_IMon){
		IMon_csv.close();
		remove("tmp_IMon.csv");
	}
	return 0;
}