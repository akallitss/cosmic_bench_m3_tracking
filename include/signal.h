#ifndef signal_h
#define signal_h
#include "Tsignal_R.h"
#include "detector.h"
#include "tomography.h"

#include <TProfile.h>
#include <TFile.h>

#include <string>
#include <map>

#include <boost/property_tree/ptree.hpp>

using std::string;
using std::map;

using boost::property_tree::ptree;

class Signal: public Tsignal_R, public CosmicBench{
	public:
		Signal(string configFilePath);
		Signal(ptree config_tree);
		~Signal();
		void MultiCluster();
		//void ElecToAnalyse();
		//void ElecToRays(string outFileName);
		void EventDisplay(int evn_min = 0, int evn_max = 20);
		void HoughTracking(long event_nb);
		map<int,TProfile*> SignalOverNoise();
		void SignalOverNoiseDisplay();
		void SignalDispersion();
	protected:
		string analyseTree;
		bool use_srf;
		Tomography::elec_type electronic_type;
		string data_file_basename;
		string signalName;
		string PedName;
		string RMSName;
		long max_event;
		int data_file_first;
		int data_file_last;
		map<int,Tomography::det_type> det_type_by_asic;
		map<int,int> det_n_by_asic;
		TFile * fIn;
		bool exists;
};
#endif