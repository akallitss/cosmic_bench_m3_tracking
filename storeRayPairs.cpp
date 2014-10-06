#include "analyse.h"
#include "TROOT.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <TH2D.h>
#include <TCanvas.h>

using std::cout;
using std::endl;
using std::ostringstream;

int main(int argc, char ** argv){
	Analyse * analysePb = new Analyse("/home/irfulx176/mnt/sbouteil/Documents/deviation/config_Pb.cfg");
	Analyse * analyseNoPb = new Analyse("/home/irfulx176/mnt/sbouteil/Documents/deviation/config_noPb.cfg");
	float z_Pb = 1553;
	analysePb->AbsorptionFluxMapNorm(z_Pb,analyseNoPb->AbsorptionFluxMap(z_Pb,60),60);
	theApp->Run(true);
	delete analysePb; delete analyseNoPb;
	return 0;
}