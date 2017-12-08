#include "analyse.h"
//#include <TROOT.h>
#include <TRint.h>
#include <iostream>
#include <sstream>
#include <TH2D.h>
#include <TCanvas.h>
#include "acceptanceFunction.h"
#include "tomography.h"
#include <TMath.h>
#include <TEllipse.h>
#include <csignal>
#include <cstdlib>
//#include <stdio.h>
#include <unistd.h>

using std::cout;
using std::endl;
using std::flush;
using std::ostringstream;

using TMath::Pi;

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
	Analyse * blah = new Analyse(config_file.str());//"/home/irfulx176/mnt/sbouteil/Documents/deviation/config.cfg");
	string efficacity = "efficacity";
	string eff2D = "eff2D";
	string residus = "residus";
	string restime = "restime";
	string ampltime = "ampltime";
	string fluxMap = "fluxmap";
	string tomoAbs = "tomoAbs";
	string raypairs = "raypairs";
	string bugtest = "bugtest";
	string srf = "srf";
	string correlation = "correlation";
	string eventdisplay = "eventdisplay";
	string eventdisplaymult = "eventdisplaymult";
	string SoN = "SoN";
	//string acceptance = "acceptance";
	string watto = "watto";
	string scanpyramids = "scanpyramids";
	string rays = "rays";
	if(argv[2] == efficacity){
		blah->Efficacity();
	}
	else if(argv[2] == residus){
		blah->Residus_ref();
		//blah->Residus();
		//blah->Residus_ref_MT();
	}
	else if(argv[2] == restime){
		blah->Residus_time();
	}
	else if(argv[2] == ampltime){
		blah->Amplitude_time();
	}
	else if(argv[2] == eff2D){
		blah->Residus_ref_2D();
	}
	else if(argv[2] == rays){
		if(argc<4){
			cout << "you must indicate the ray file name" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		blah->ExportAbsorptionRays(argv[3]);
	}
	else if(argv[2] == fluxMap){
		if(argc<4){
			cout << "you must indicate the altitude of the flux map" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		else{
			float z = atof(argv[3]);
			float y_angle = 0;
			int nbins = -1;
			if(argc>4){
				y_angle = atof(argv[4]);
			}
			if(argc>5){
				nbins = atoi(argv[5]);
			}
			cout << "flux map param : " << endl;
			cout << "   z : " << z << "mm" << endl;
			cout << "   y_angle : " << y_angle << "°" << endl;
			cout << "   nbins : ";
			if(nbins>0) cout << nbins;
			else cout << "auto";
			cout << endl;
			blah->AbsorptionFluxMap(z,0,y_angle*Pi()/180.,nbins);
		}
	}
	else if(argv[2] == tomoAbs){
		if(argc<4){
			cout << "you must indicate the altitude of the flux map" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		else{
			float z = atof(argv[3]);
			double bench_angle = 0;
			if(argc>4){
				bench_angle = atof(argv[4]);
			}
			blah->AbsorptionFluxMapNormTheo(z,bench_angle*Pi()/180.);
		}
	}
	else if(argv[2] == scanpyramids){
		if(argc<5){
			cout << "you must indicate the angle of the telescope wrt to the ground and the multiplicity of the trigger" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		else{
			double bench_angle = atof(argv[3]);
			int mult = atoi(argv[4]);
			blah->AbsorptionFluxMapNormTheoAngle(bench_angle,mult);
		}
	}
	else if(argv[2] == raypairs){
		if(argc<4){
			cout << "you must indicate the output file name" << endl;
			Tomography::Quit();
			delete blah;
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
		if(argc<5){
			cout << "you must indicate the first and last event you wanna see" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		else{
			int i = atoi(argv[3]);
			int j = atoi(argv[4]);
			TCanvas * cDisplay = new TCanvas();
			for(int k=i;k<=j && Tomography::get_instance()->get_can_continue();k++){
				cout << "\r" << "Event : " << k << flush;
				blah->EventDisplay(k, cDisplay);
				sleep(2);
			}
			cout << endl;
		}
	}
	else if(argv[2] == eventdisplay){
		if(argc<4){
			cout << "you must indicate the index of the event you wanna see" << endl;
			Tomography::Quit();
			delete blah;
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
	/*
	else if(argv[2] == acceptance){
		acceptanceFunction chombier(0,500,0,500,1270,0,0);
		TH2D * plot = new TH2D(chombier.plot_XY(500,-200,700,500,-200,700,1600));
		plot->Draw("colz");
	}
	*/
	else if(argv[2] == watto){
		if(argc<9){
			cout << "you must indicate in order : flux map altitude, ellipse x center, ellipse y center, ellipse x axis length, ellipse y axis length, ellipse rotation angle" << endl;
			Tomography::Quit();
			delete blah;
			return 1;
		}
		else{
			float z = atof(argv[3]);
			float y_angle = 0;
			if(argc>9){
				y_angle = atof(argv[9]);
			}
			cout << "flux map param : " << endl;
			cout << "   z : " << z << "mm" << endl;
			cout << "   y_angle : " << y_angle << "°" << endl;
			blah->WatToFluxMap(z,TEllipse(atof(argv[4]),atof(argv[5]),atof(argv[6]),atof(argv[7]),0,360,atof(argv[8])),0,0,y_angle*Pi()/180.);
		}
	}
	else{
		cout << "function not found" << endl;
		cout << blah->get_z_Up() << endl;
		cout << blah ->get_z_Down() << endl;
		//blah->MultiGenDebug(0);
		Tomography::Quit();
		delete blah;
		return 1;
	}
	Tomography::get_instance()->Run();
	Tomography::Quit();
	delete blah;
	return 0;
}
