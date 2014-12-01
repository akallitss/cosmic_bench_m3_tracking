#include "datareader.h"
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

using std::string;
using std::map;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::ifstream;
using std::cout;
using std::endl;

using boost::property_tree::ptree;

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the Run caracs" << endl;
		return 1;
	}

	string config_file = argv[1];
	ptree config_tree;
	read_json(config_file, config_tree);
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	int CM_N = 0;
	int MG_N = 0;
	map<int,string> det_type_by_asic;
	map<int,int> det_n_by_asic;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		det_type_by_asic[child.second.get<int>("asic_n")] = "CM";
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("cm_n");
		CM_N++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		det_type_by_asic[child.second.get<int>("asic_n")] = "MG";
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("mg_n");
		MG_N++;
	}
	if((total_CM_N!=CM_N) || (total_MG_N!=MG_N)){
		cout << "problem in detectors number" << endl;
		return 1;
	}
	string electronic_type = config_tree.get<string>("electronic_type");
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
	if(electronic_type == "dream"){
		blah = new DreamDataReader(signalName,PedName,RMSName,det_type_by_asic,det_n_by_asic,exists,exists,exists,max_event);
	}
	else if(electronic_type == "feminos"){
		blah = new FeminosDataReader(signalName,PedName,RMSName,det_type_by_asic,det_n_by_asic,exists,exists,exists,max_event);
	}
	else{
		cout << "electronic type : " << electronic_type << " unknown !" << endl;
		return 1;
	}
	int first_data_file = config_tree.get<int>("data_file_first");
	int last_data_file = config_tree.get<int>("data_file_last") + 1;
	for(int i=first_data_file;i<last_data_file;i++){
		ostringstream dataFileName;
		dataFileName << data_file_basename << setw(3) << setfill('0') << i;
		if(electronic_type == "dream") dataFileName << ".fdf";
		else if(electronic_type == "feminos") dataFileName << ".aqs";
		else dataFileName << ".txt";
		blah->add_file_to_process(dataFileName.str());
	}
	blah->process();
	if(!ped_done) blah->compute_ped();
	blah->do_ped_sub();
	blah->do_common_noise_sub();
	if(compute_rms) blah->compute_RMSPed();
	return 0;
}