#ifndef datareader_h
#define datareader_h
#include <string>
#include <vector>
#include <map>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <fstream>
#include "detector.h"
#include "tomography.h"

using std::string;
using std::vector;
using std::map;
using std::ifstream;

class DataReader{
	public:
		DataReader(string baseFileName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		DataReader(string signalName, string pedName, string RMSName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		virtual ~DataReader();
		void add_file_to_process(string inFileName);
		virtual void process() = 0;
		void compute_ped();
		void do_ped_sub(string ped_file = "");
		void do_common_noise_sub();
		void compute_RMSPed();
		static const int Nsample = Tomography::Nsample;
		static const int Nstrip_MG = MG_Detector::Nstrip;
		static const int Nstrip_CM = CM_Detector::Nstrip;
		virtual map<Tomography::det_type,vector<vector<vector<double> > > > read_event(ifstream * file,long& event_nb, bool fill_tree = true) = 0;
		map<Tomography::det_type,vector<vector<vector<double> > > > read_event(ifstream * file,int& event_nb, bool fill_tree = true);
	protected:
		void read_ped(string ped_file = "");
		void Fill();
		void Write();
		virtual void read_file(string file_name,long evn_offset) = 0;
		void reset_tree_leaf();
		vector<string> file_names;
		TFile * outFile;
		TTree * outTree;
		int Nevent;
		double evttime;
		int TsampleNum[Tomography::Nsample];
		float (*StripAmpl_MG)[61][Tomography::Nsample];
		float (*StripAmpl_CM)[64][Tomography::Nsample];
		float (*StripAmpl_MG_ped)[61][Tomography::Nsample];
		float (*StripAmpl_CM_ped)[64][Tomography::Nsample];
		float (*StripAmpl_MG_corr)[61][Tomography::Nsample];
		float (*StripAmpl_CM_corr)[64][Tomography::Nsample];
		float (*Pedestal_MG)[61];
		float (*Pedestal_CM)[64];
		unsigned int MG_N;
		unsigned int CM_N;
		map<int,Tomography::det_type> det_type_by_asic;
		map<int,int> det_n_by_asic;
		bool is_first;
		string outFileName;
		string PedFileName;
		string RMSPedFileName;
		bool exists;
		bool ped_done;
		bool cns_done;
		Tomography::elec_type DAQType;
		long max_event;
		long global_offset;
		TBranch * dumb_branch;
	private:
		void Init(string signalName, string pedName, string RMSName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);

};

class DreamDataReader: public DataReader{
	public:
		DreamDataReader(string baseFileName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		DreamDataReader(string signalName, string pedName, string RMSName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		~DreamDataReader();
		void process();
		map<Tomography::det_type,vector<vector<vector<double> > > > read_event(ifstream * file,long& event_nb, bool fill_tree = true);
		int get_first_event_nb(string file_name);
	protected:
		void read_file(string file_name,long evn_offset);
		int mapping(Tomography::det_type det_type, int channel);
};

class FeminosDataReader: public DataReader{
	public:
		FeminosDataReader(string baseFileName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		FeminosDataReader(string signalName, string pedName, string RMSName, map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, long max_event_ = -1);
		~FeminosDataReader();
		void process();
		map<Tomography::det_type,vector<vector<vector<double> > > > read_event(ifstream * file,long& event_nb, bool fill_tree = true);
		int get_first_event_nb(string file_name);
	protected:
		void read_file(string file_name,long evn_offset);
		int mapping(Tomography::det_type det_type, int channel);
};

#endif