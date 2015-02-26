#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <sys/timex.h>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "CAEN_comm.h"

#include <TFile.h>
#include <TTree.h>

using std::cout;
using std::endl;
using std::string;
using std::map;
using std::vector;
using std::setw;
using std::left;
using std::showpos;
using std::noshowpos;
using std::find;
using std::ostringstream;

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
	CAEN_Comm * blah = new CAEN_Comm("132.166.15.231","admin","mcube");
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
	wait_time.tv_sec = 0;
	struct ntptimeval current_time;

	TFile * fOut = new TFile(config_tree.get<string>("filename").c_str(),"RECREATE");
	TTree * outTree = new TTree("T","HV");
	outTree->Branch("time",&current_time.time.tv_sec,"time/L");
	ostringstream branch_name;
	branch_name << "name[" << total_channel << "]/C";
	outTree->Branch("name",chan_names,branch_name.str().c_str());

	map<CAEN_Ch::Param,CAEN_Ch::param_value*> channel_params;
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
	}

	for(int i=0;i<config_tree.get<int>("duration");i++){
		ntp_gettime(&current_time);
		wait_time.tv_nsec = 1000*(1000000-current_time.time.tv_usec);
		nanosleep(&wait_time,&remaining_time);
		for(map<CAEN_Ch::Param,CAEN_Ch::param_value*>::iterator params_it = channel_params.begin();params_it!=channel_params.end();++params_it){
			channel_values[params_it->first] = blah->get_Ch_param(used_channel,params_it->first);

		}
		for(map<CAEN_Ch::Param,map<int,map<int,CAEN_Ch::param_value> > >::iterator values_it = channel_values.begin();values_it!=channel_values.end();++values_it){
			current_channel = 0;
			for(map<int,vector<int> >::iterator map_it = used_channel.begin();map_it!=used_channel.end();++map_it){
				for(vector<int>::iterator vec_it = (map_it->second).begin();vec_it!=(map_it->second).end();++vec_it){
					channel_params[values_it->first][current_channel] = (values_it->second)[map_it->first][*vec_it];
					current_channel++;
				}
			}
		}
		outTree->Fill();
	}
	outTree->Write();
	fOut->Close();
	delete blah;
	return 0;
}