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
	if(argc<3){
		cout << "You must indicate a config file which contains the Cosmic Bench caracs and what you want to plot" << endl;
		return 1;
	}
	const int n = 100;
	char path[n];
	ostringstream config_file;
	config_file << getcwd(path,n) << "/" << argv[1];
	Tomography::Init(config_file.str());
	Carac * blah = new Carac(config_file.str());//"/home/irfulx176/mnt/sbouteil/Documents/deviation/config.cfg");
	string residus = "residus";
	if(argv[2] == residus){
		//blah->Residus();
		if(argc<3) blah->Residus_ref();
		else{
			double theta_max = atof(argv[3]);
			cout << "using track with tan(theta) < " << theta_max << endl;
			blah->Residus_ref(theta_max);
		}
	}
	else{
		cout << "function not found" << endl;
		Tomography::Quit();
		delete blah;
		return 1;
	}
	Tomography::get_instance()->Run();
	Tomography::Quit();
	delete blah;
	return 0;
}
