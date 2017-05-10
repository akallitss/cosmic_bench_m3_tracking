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
	if(argc<3){
		cout << "You must indicate a config file which contains the Cosmic Bench caracs and what you want to plot" << endl;
		return 1;
	}
	const int n = 100;
	char path[n];
	ostringstream config_file;
	config_file << getcwd(path,n) << "/" << argv[1];
	Tomography::Init(config_file.str());
	Signal * blah = new Signal(config_file.str());
	string multicluster = "multicluster";
	string multicluster_raw = "rawcluster";
	string hough = "hough";
	string SoB = "SoB";
	//string EtA = "EtA";
	string display = "display";
	string dispersion = "dispersion";
	string debug = "debug";
	string conv = "conv";
	string noiselevels = "noiselevels";
	bool is_interactive = true;
	if(argv[2] == multicluster){
		blah->MultiCluster();
		is_interactive = false;
	}
	else if(argv[2] == multicluster_raw){
		blah->MultiCluster_raw();
		is_interactive = false;
	}
	else if(argv[2] == noiselevels){
		blah->NoiseLevels();
	}
	else if(argv[2] == SoB){
		blah->SignalOverNoiseDisplay();
	}
	else if(argv[2] == conv){
		blah->ConvClusterTest();
	}
	else if(argv[2] == hough){
		if(argc<4){
			cout << "you must indicate an event number" << endl;
			return 1;
		}
		int event_nb = atoi(argv[3]);
		blah->HoughTracking(event_nb);
	}
	else if(argv[2] == debug){
		if(argc<4){
			cout << "you must indicate an event number" << endl;
			return 1;
		}
		int event_nb = atoi(argv[3]);
		blah->DebugHoles(event_nb);
	}
	else if(argv[2] == display){
		if(argc<6){
			cout << "you must indicate the signal you want to plot (raw,ped,corr) and a start and stop event number" << endl;
			return 1;
		}
		string signal_correction_str = argv[3];
		Tomography::signal_type signal_correction;
		if(signal_correction_str == "raw"){
			signal_correction = Tomography::raw;
		}
		else if(signal_correction_str == "ped"){
			signal_correction = Tomography::ped;
		}
		else signal_correction = Tomography::corr;
		int event_nb_start = atoi(argv[4]);
		int event_nb_stop = atoi(argv[5]);
		blah->EventDisplay(event_nb_start,event_nb_stop,signal_correction);
	}
	else if(argv[2] == dispersion){
		blah->SignalDispersion();
	}
	if(is_interactive) Tomography::get_instance()->Run();
	Tomography::Quit();
	delete blah;
	return 0;
}
