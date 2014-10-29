#include "signal.h"
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include "TRint.h"

using std::ostringstream;
using std::cout;
using std::endl;


int main(int argc, char ** argv){
	if(argc<3){
		cout << "You must indicate a config file which contains the Cosmic Bench caracs and what you want to plot" << endl;
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
	Signal * blah = new Signal(config_file.str());
	string multicluster = "multicluster";
	string hough = "hough";
	if(argv[2] == multicluster){
		blah->MultiCluster();
	}
	else if(argv[2] == hough){
		if(argc<4){
			cout << "you must indicate an event number" << endl;
			return 1;
		}
		int event_nb = atoi(argv[3]);
		blah->HoughTracking(event_nb);
		theApp->Run(true);
	}
	//delete blah; delete theApp;
	return 0;
}