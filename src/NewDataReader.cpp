#include <string>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "datareader.h"

using std::string;
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
	DataReader blah(config_tree,true);
	blah.process();
	blah.compute_ped();
	blah.read_ped();
	blah.do_ped_sub();
	blah.do_common_noise_sub();
	blah.compute_RMSPed();
	return 0;
}