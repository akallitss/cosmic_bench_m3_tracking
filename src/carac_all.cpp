#include "carac.h"
//#include <TROOT.h>
#include <TRint.h>
#include <iostream>
#include "tomography.h"
#include <csignal>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using std::cout;
using std::endl;
using std::ostringstream;

int main(int argc, char ** argv){
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = Tomography::signal_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);
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
	TRint * theApp;
	if(!Tomography::is_batch) theApp = new TRint("Rint",&argcR,argvR,0,0,true);
	Carac * blah = new Carac(config_file.str());//"/home/irfulx176/mnt/sbouteil/Documents/deviation/config.cfg");
	string residus = "residus";
	if(argv[2] == residus){
		//blah->Residus();
		blah->Residus_ref();
	}
	else{
		cout << "function not found" << endl;
		delete blah; delete theApp;
		return 1;
	}
	if(Tomography::is_batch) Tomography::save_canvases();
	else {
		theApp->Run(true);
		delete theApp;
	}
	delete blah;
	return 0;
}
