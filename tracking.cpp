#include "analyse.h"
#include <TROOT.h>
#include <TRint.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <TH2D.h>

using std::cout;
using std::endl;
using std::ostringstream;

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
	Analyse * blah = new Analyse(config_file.str());//"/home/irfulx176/mnt/sbouteil/Documents/deviation/config.cfg");
	string efficacity = "efficacity";
	string residus = "residus";
	string fluxMap = "fluxmap";
	string raypairs = "raypairs";
	string bugtest = "bugtest";
	if(argv[2] == efficacity){
		blah->Efficacity();
	}
	else if(argv[2] == residus){
		//blah->Residus();
		blah->Residus_ref();
	}
	else if(argv[2] == fluxMap){
		if(argc<4){
			cout << "you must indicate the altitude of the flux map" << endl;
			delete blah; delete theApp;
			return 1;
		}
		else{
			float z = atof(argv[3]);
			blah->AbsorptionFluxMap(z);
		}
	}
	else if(argv[2] == raypairs){
		if(argc<4){
			cout << "you must indicate the output file name" << endl;
			delete blah; delete theApp;
			return 1;
		}
		else{
			string outName = argv[3];
			blah->StoreRayPairs(outName);
		}
	}
	else if(argv[2] == bugtest){
		blah->bugtest();
	}
	else{
		cout << "function not found" << endl;
		cout << blah->get_z_Up() << endl;
		cout << blah ->get_z_Down() << endl;
		//blah->MultiGenDebug(0);
		delete blah; delete theApp;
		return 1;
	}
	theApp->Run(true);
	delete blah; delete theApp;
	return 0;
}
