#include "analyse.h"
#include <TROOT.h>
#include <TRint.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TCanvas.h>
#include "acceptanceFunction.h"

using std::cout;
using std::endl;
using std::flush;
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
	string eff2D = "eff2D";
	string residus = "residus";
	string fluxMap = "fluxmap";
	string tomoAbs = "tomoAbs";
	string raypairs = "raypairs";
	string bugtest = "bugtest";
	string srf = "srf";
	string correlation = "correlation";
	string eventdisplay = "eventdisplay";
	string eventdisplaymult = "eventdisplaymult";
	string SoN = "SoN";
	string acceptance = "acceptance";
	if(argv[2] == efficacity){
		blah->Efficacity();
	}
	else if(argv[2] == residus){
		//blah->Residus();
		blah->Residus_ref();
	}
	else if(argv[2] == eff2D){
		blah->Residus_ref_2D();
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
	else if(argv[2] == tomoAbs){
		if(argc<4){
			cout << "you must indicate the altitude of the flux map" << endl;
			delete blah; delete theApp;
			return 1;
		}
		else{
			float z = atof(argv[3]);
			blah->AbsorptionFluxMapNormTheo(z);
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
	else if(argv[2] == srf){
		if(argc<4){
			blah->CalcStripResponseFunction();
		}
		else{
			int i = atoi(argv[3]);
			blah->CalcStripResponseFunction(i);
		}
	}
	else if(argv[2] == eventdisplaymult){
		if(argc<4){
			cout << "you must indicate the number of the event you wanna see" << endl;
			delete blah; delete theApp;
			return 1;
		}
		else{
			int i = atoi(argv[3]);
			TCanvas * cDisplay = new TCanvas();
			for(int j=0;j<i;j++){
				cout << "\r" << "Event : " << j << flush;
				blah->EventDisplay(j, cDisplay);
				gSystem->Sleep(2000);
			}
			cout << endl;
		}
	}
	else if(argv[2] == eventdisplay){
		if(argc<4){
			cout << "you must indicate the index of the event you wanna see" << endl;
			delete blah; delete theApp;
			return 1;
		}
		else{
			int i = atoi(argv[3]);
			blah->EventDisplay(i);
		}
	}
	else if(argv[2] == correlation){
		blah->Correlation();
	}
	else if(argv[2] == SoN){
		blah->SignalOverNoise();
	}
	else if(argv[2] == acceptance){
		acceptanceFunction chombier(0,500,0,500,1270,0,0);
		chombier.plot_3D();
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
