#include "Signal.h"
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <csignal>
#include "TRint.h"

using std::ostringstream;
using std::cout;
using std::endl;


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
	sigIntHandler.sa_handler = Tomography::signal_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);
	Signal * blah = new Signal(config_file.str());
	string multicluster = "multicluster";
	string hough = "hough";
	string SoB = "SoB";
	//string EtA = "EtA";
	string display = "display";
	string dispersion = "dispersion";
	bool is_interactive = true;
	if(argv[2] == multicluster){
		blah->MultiCluster();
		is_interactive = false;
	}
	/*
	else if(argv[2] == EtA){
		blah->ElecToAnalyse();
	}
	*/
	else if(argv[2] == SoB){
		blah->SignalOverNoiseDisplay();
	}
	else if(argv[2] == hough){
		if(argc<4){
			cout << "you must indicate an event number" << endl;
			return 1;
		}
		int event_nb = atoi(argv[3]);
		blah->HoughTracking(event_nb);
	}
	else if(argv[2] == display){
		if(argc<4){
			cout << "you must indicate an event number" << endl;
			return 1;
		}
		int event_nb = atoi(argv[3]);
		blah->EventDisplay(0,event_nb);
	}
	else if(argv[2] == dispersion){
		blah->SignalDispersion();
	}
	if(Tomography::is_batch) Tomography::save_canvases();
	else {
		if(is_interactive) theApp->Run(true);
		delete theApp;
	}
	delete blah;
	return 0;
}