#include "analyse.h"
#include "tomography.h"
//#include <TROOT.h>
#include <TRint.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>

using std::cout;
using std::endl;
using std::ostringstream;
using TMath::Sqrt;

int main(int /*argc*/, char ** argv){
	int argcR = 1;
	char * argvR[1];
	argvR[0] = argv[0];
	TRint * theApp = new TRint("Rint",&argcR,argvR,0,0,true);
	Analyse * analysePb = new Analyse("/home/irfulx176/mnt/sbouteil/Documents/deviation/config_Pb.cfg");
	Analyse * analyseNoPb = new Analyse("/home/irfulx176/mnt/sbouteil/Documents/deviation/config_noPb.cfg");
	double z_Pb = 1553;
	int nbins = 50;
	TCanvas * cNoPb = new TCanvas("fluxMap_noPb","fluxMap_noPb");
	TH2D * fluxMapNoPb = analyseNoPb->AbsorptionFluxMap(z_Pb, cNoPb);
	fluxMapNoPb->SetTitle("fluxMapBg");
	fluxMapNoPb->SetName("fluxMapBg");
	cNoPb->cd();
	fluxMapNoPb->Draw("COLZ");
	cNoPb->Modified();
	cNoPb->Update();
	TCanvas * cPb = new TCanvas("fluxMap_Pb","fluxMap_Pb");
	TCanvas * cDiff = new TCanvas("fluxMap_Diff","fluxMap_Diff");
	TCanvas * cSigma = new TCanvas("fluxMap_Sigma","fluxMap_Sigma");
	analysePb->AbsorptionFluxMapNorm(z_Pb,fluxMapNoPb, nbins, cPb, cDiff, cSigma);
	
	theApp->Run(true);
	delete analysePb; delete analyseNoPb;
	delete fluxMapNoPb;
	delete cNoPb; delete cPb; delete cSigma; delete cDiff;
	delete theApp;
	return 0;
}
