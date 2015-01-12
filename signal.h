#ifndef signal_h
#define signal_h
#include "Tsignal.h"
#include "detector.h"

#include <TProfile.h>

#include <string>
#include <map>

using std::string;
using std::map;

class Signal: public Tsignal, public CosmicBench{
	public:
		Signal(string configFilePath);
		~Signal();
		void MultiCluster();
		void ElecToAnalyse();
		void HoughTracking(int event_nb);
		map<int,TProfile*> SignalOverNoise();
		void SignalOverNoiseDisplay();
	protected:
		string analyseTree;
		bool use_srf;
		string electronic_type;
		string data_file_basename;
		string signalName;
		string PedName;
		string RMSName;
		int max_event;
		int data_file_first;
		int data_file_last;
		map<int,string> det_type_by_asic;
		map<int,int> det_n_by_asic;
};
#endif