#include "signal.h"
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include "TRint.h"
#include <TProfile.h>
#include <map>
#include <TCanvas.h>

using std::map;
using std::ostringstream;
using std::cout;
using std::endl;


int main(int /*argc*/, char ** argv){
	int argcR = 1;
	char * argvR[1];
	argvR[0] = argv[0];
	TRint * theApp = new TRint("Rint",&argcR,argvR,0,0,true);
	map<int,string> config_name;
	config_name[1] = "config/config_testCapa1.cfg";
	config_name[2] = "config/config_testCapa2.cfg";
	config_name[3] = "config/config_testCapa3.cfg";
	config_name[4] = "config/config_testCapa4.cfg";
	map<int,string> capacitance;
	capacitance[0] = "220pF";
	capacitance[1] = "470pF";
	capacitance[2] = "1nF";
	capacitance[3] = "2.2nF";
	map<int,TProfile*> SoN;
	Signal * blah = new Signal(config_name[1]);
	map<int,TProfile*> current_SoN = blah->SignalOverNoise();
	for(map<int,TProfile*>::iterator it = current_SoN.begin();it!=current_SoN.end();++it){
		SoN[it->first] = new TProfile(*(it->second));
		SoN[it->first]->SetDirectory(0);
	}
	delete blah;
	TCanvas * cDisplay = new TCanvas();
	cDisplay->Divide(4);
	for(int i=0;i<4;i++){
		cDisplay->cd(i+1);
		SoN[i]->SetTitle(capacitance[i].c_str());
		SoN[i]->Draw();
	}
	cDisplay->Modified();
	cDisplay->Update();
	for(int i=2;i<5;i++){
		blah = new Signal(config_name[i]);
		current_SoN = blah->SignalOverNoise();
		for(int j=0;j<4;j++){
			SoN[j]->Add(current_SoN[(9-i+j)%4]);
		}
		delete blah;
		for(int j=0;j<4;i++){
			cDisplay->cd(j+1);
			SoN[j]->Draw();
		}
		cDisplay->Modified();
		cDisplay->Update();
	}
	for(int i=0;i<4;i++){
		cDisplay->cd(i+1);
		SoN[i]->Draw();
	}
	cDisplay->Modified();
	cDisplay->Update();
	theApp->Run(true);
	return 0;
}