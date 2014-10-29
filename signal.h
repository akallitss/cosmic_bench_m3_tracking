#ifndef signal_h
#define signal_h
#include "Tsignal.h"
#include "detector.h"

#include <string>

using std::string;

class Signal: public Tsignal, public CosmicBench{
	public:
		Signal(string configFilePath);
		~Signal();
		void MultiCluster();
		void HoughTracking(int event_nb);
	protected:
		string analyseTree;
		bool use_srf;
};
#endif