#ifndef datareader_h
#define datareader_h
#include <string>
#include <vector>
#include <map>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>

using std::string;
using std::vector;
using std::map;

class DataReader{
	public:
		DataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, int max_event_ = -1);
		~DataReader();
		void add_file_to_process(string inFileName);
		virtual void process() = 0;
		void compute_ped();
		void do_ped_sub(string ped_file = "");
		void do_common_noise_sub();
		void compute_RMSPed();
		static const int Nsample;
		static const int Nstrip_MG;
		static const int Nstrip_CM;
	protected:
		void read_ped(string ped_file = "");
		void Fill();
		void Write();
		virtual void read_file(string file_name,int evn_offset) = 0;
		void reset_tree_leaf();
		vector<string> file_names;
		TFile * outFile;
		TTree * outTree;
		int Nevent;
		int TsampleNum[32];
		float (*StripAmpl_MG)[61][32];
		float (*StripAmpl_CM)[64][32];
		float (*StripAmpl_MG_ped)[61][32];
		float (*StripAmpl_CM_ped)[64][32];
		float (*StripAmpl_MG_corr)[61][32];
		float (*StripAmpl_CM_corr)[64][32];
		float (*Pedestal_MG)[61];
		float (*Pedestal_CM)[64];
		unsigned int MG_N;
		unsigned int CM_N;
		map<int,string> det_type_by_asic;
		map<int,int> det_n_by_asic;
		bool is_first;
		string outFileName;
		string PedFileName;
		string RMSPedFileName;
		bool exists;
		bool ped_done;
		bool cns_done;
		string DAQType;
		int max_event;
		int global_offset;
		TBranch * dumb_branch;

};

class DreamDataReader: public DataReader{
	public:
		DreamDataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, int max_event_ = -1);
		~DreamDataReader();
		void process();
	protected:
		void read_file(string file_name,int evn_offset);
		int mapping(string det_type, int channel);
};

class FeminosDataReader: public DataReader{
	public:
		FeminosDataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_=false,bool ped_done_=false,bool cns_done_=false, int max_event_ = -1);
		~FeminosDataReader();
		void process();
	protected:
		void read_file(string file_name,int evn_offset);
		int get_first_event_nb(string file_name);
		int mapping(string det_type, int channel);
};

#endif