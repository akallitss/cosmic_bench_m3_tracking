#ifndef analyse_h
#define analyse_h
#include "T.h"
#include "detector.h"
#include <TTree.h>
#include <string>
#include <TH2D.h>
#include <TCanvas.h>

using std::string;

class Analyse: public T, public CosmicBench{
	public:
		Analyse(string configFilePath);
		~Analyse();
		void Efficacity();
		void Residus();
		void Residus_ref();
		void Residus_ref_2D();
		TH2D * AbsorptionFluxMap(double z, TCanvas * c1 = 0);
		void AbsorptionFluxMapNorm(double z,TH2D * background, int nbins = 100, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0);
		void AbsorptionFluxMapNormTheo(double z, TCanvas * c1 = 0, TCanvas * c2 = 0, TCanvas * c3 = 0, TCanvas * c4 = 0);
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

};

#endif