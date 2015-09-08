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

class Analyse: public Tanalyse_R, public CosmicBench{
	public:
		Analyse(string configFilePath);
		Analyse(ptree config_tree);
		~Analyse();
		void Efficacity();
		void Residus();
		void Residus_ref();
		double Residus_ref_cost();
		void Residus_ref_2D();
		TH2D * AbsorptionFluxMap(double z, TCanvas * c1 = 0, double y_angle = 0);
		void ExportAbsorptionRays(string outFileName);
		void AbsorptionFluxMapNorm(double z,TH2D * background, int nbins = 100, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0);
		void AbsorptionFluxMapNormTheo(double z, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0, TCanvas * c4 = 0);
		void WatToFluxMap(double z, TEllipse el, TCanvas * c1 = 0, TCanvas * c2 = 0, double y_angle = 0);
		void StoreRayPairs(string outFileName);
		void CalcStripResponseFunction(int bin_nb = 0);
		void Correlation();
		double get_z_Up() const;
		double get_z_Down() const;
		//void MultiGenDebug(int i);
		void bugtest();
		void EventDisplay(long event_nb, TCanvas * c1 = 0);
		void SignalOverNoise();
	protected:
		string signal_file_name;
		long max_event;
		TFile * f;

};

#endif