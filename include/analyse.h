#ifndef analyse_h
#define analyse_h

#include "Tanalyse_R.h"
#include "detector.h"

#include <TEllipse.h>

#include <string>
#include <boost/property_tree/ptree.hpp>

using std::string;
using boost::property_tree::ptree;

class TCanvas;
class TH2D;
class TFile;

//use analyse.root files to do analysis involving tracking
class Analyse: public Tanalyse_R, public CosmicBench{
	public:
		//build object given the file name containing the json config tree
		Analyse(string configFilePath);
		//build object given the config tree
		Analyse(ptree config_tree);
		~Analyse();
		//compute efficiency variations over time
		void Efficacity();
		//compute residual (x_det - x_track) for all detectors, probed detectors are included in the tracking
		void Residus();
		//compute alignment properties (residual in function of interesting variables) and efficiency map for non_ref detectors, probed detectors are excluded of the tracking
		void Residus_ref();
		//currently bugged do no use
		void Residus_ref_MT();
		//compute efficiency
		void Residus_time();
		//compute aplitude and trigger rate variations through time
		void Amplitude_time();
		//attempt to compute a missalignement representative variable
		double Residus_ref_cost();
		//compute 2D efficiency map for 2D detectors
		void Residus_ref_2D();
		//compute the muon flux map at the plane which cross the z axis (bench axis) at z with an angle around the y axis y_angle
		//the output canvas can be passed or automatically created
		//the resulting histogram bin number can be passed or computed using the number of tree entries
		TH2D * AbsorptionFluxMap(double z, TCanvas * c1 = 0, double y_angle = 0, int nbins = -1);
		//export straight muon tracks to the file outFileName, more details in Tray.h for the file format
		void ExportAbsorptionRays(string outFileName);
		//Similar to AbsorptionFluxMap but take a pre-computed background flux to normalize the computed flux
		void AbsorptionFluxMapNorm(double z,TH2D * background, int nbins = 100, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0);
		//Similar to AbsorptionFluxMap but also compute a theoretical background flux using acceptanceFunction (cf. acceptanceFunction.h)
		//bench_angle is the angle between the z axis and the zenith
		void AbsorptionFluxMapNormTheo(double z, double bench_angle = 0, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0, TCanvas * c4 = 0);
		//Similar to AbsorptionFluxMapNormTheo but output the resulting flux in angular coordinates
		void AbsorptionFluxMapNormTheoAngle(double bench_angle, int mult, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0, TCanvas * c4 = 0);
		//Used for WatTo data analysis, the ellipse correspond to the tank area. It shows the dynamic measurement of the tank absorption
		void WatToFluxMap(double z, TEllipse el, TCanvas * c1 = 0, TCanvas * c2 = 0, double y_angle = 0);
		//export muon tracks undergoing a scattering to the file outFileName by storing the information of both half-tracks
		void StoreRayPairs(string outFileName);
		//attempt to adapt the pad response function developped by ILC TPC
		void CalcStripResponseFunction(int bin_nb = 0);
		//compute correlations between the two coordinates of a 2D detector such as the MG2D
		void Correlation();
		//void MultiGenDebug(int i);
		void bugtest();
		//display bench with tracks and cluster positions
		void EventDisplay(long event_nb, TCanvas * c1 = 0);
		//compute signal and noise distribution to get the S/N ratio
		void SignalOverNoise();
	protected:
		string signal_file_name;
		long max_event;
		TFile * f;

};

#endif
