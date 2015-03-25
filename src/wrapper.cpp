#include <string>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cstdio>

#include "tomography.h"
#include "signal.h"

using std::string;
using std::cout;
using std::endl;
using std::remove;

using boost::property_tree::ptree;

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the data run and ped run carac" << endl;
		return 1;
	}

	string config_file_wrapper = argv[1];
	ptree config_tree_wrapper;
	read_json(config_file_wrapper, config_tree_wrapper);
	ptree config_tree_bench;
	read_json("config_default.cfg", config_tree_bench);
	config_tree_bench.put<string>("Ped",config_tree_wrapper.get<string>("Ped"));
	config_tree_bench.put<string>("RMSPed",config_tree_wrapper.get<string>("RMSPed"));
	remove(config_tree_bench.get<string>("signal_file"));
	config_tree_bench.put<string>("data_file_basename",config_tree_wrapper.get<string>("pedrun_name"));
	config_tree_bench.put<int>("data_file_first",config_tree_wrapper.get<int>("pedrun_min"));
	config_tree_bench.put<int>("data_file_last",config_tree_wrapper.get<int>("pedrun_max"));
	Tomography::process_elec_files(config_tree_bench);
	remove(config_tree_bench.get<string>("signal_file"));
	config_tree_bench.put<string>("data_file_basename",config_tree_wrapper.get<string>("datarun_name"));
	config_tree_bench.put<int>("data_file_first",config_tree_wrapper.get<int>("datarun_min"));
	config_tree_bench.put<int>("data_file_last",config_tree_wrapper.get<int>("datarun_max"));
	Signal * blah = new Signal(config_tree_bench);
	blah->ElecToRays(config_tree_wrapper.get<string>("outTree"));
	delete blah;
	return 0;
}