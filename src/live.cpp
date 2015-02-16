#include "liveDisplay.h"

#include <TROOT.h>
#include <TRint.h>

#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <iostream>

using std::cout;
using std::endl;
using std::ostringstream;

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the Cosmic Bench caracs" << endl;
		return 1;
	}
	int n = 100;
	char path[n];
	ostringstream config_file;
	config_file << getcwd(path,n) << "/" << argv[1];
	int argcR = 1;
	char * argvR[1];
	argvR[0] = argv[0];
	TRint * theApp = new TRint("Rint",&argcR,argvR,0,0,true);
	liveDisplay * this_display = new liveDisplay(config_file.str());
	this_display->flux_map(1550);
	theApp->Run(true);
	return 0;

}
