#ifndef datareader_h
#define datareader_h
#include "tomography.h"
#include <vector>
#include <map>
#include <string>
#include <boost/property_tree/ptree.hpp>

using std::vector;
using std::map;
using std::string;

using boost::property_tree::ptree;

class Tsignal_W;
class ElecReader;

class DataReader{
	public:
		DataReader();
		DataReader(ptree config_tree, bool save_to_disk);
		//DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, map<int,int> feu_id_to_n_, int first_index_, int last_index_, string PedName_, string RMSName_);
		//DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, vector<int> fem_id,int first_index_, int last_index_, string PedName_, string RMSName_);
		~DataReader();
		void process();
		void process_event();
		void compute_ped();
		void read_ped();
		void compute_RMSPed();
		void do_ped_sub();
		void do_common_noise_sub();
		void do_ped_sub_event();
		void do_common_noise_sub_event();
		void do_ped_CMN_sub_event();
		long get_event_n() const;
		long get_max_event() const;
		double get_evttime() const;
		bool is_end() const;
		map<Tomography::det_type,vector<vector<vector<float> > > > get_data() const;
		map<Tomography::det_type,vector<vector<float> > > get_Ped() const;
		template<typename S, typename T> static map<Tomography::det_type,vector<vector<vector<S> > > > do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<T> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
	protected:
		//void Init(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string PedName_, string RMSName_, string outFileName = "", long max_event_ = -1);
		static int Dream_mapping(Tomography::det_type det,int channel);
		static int Feminos_mapping(Tomography::det_type det,int channel);
		int (*mapping)(Tomography::det_type,int);
		ElecReader * reader;
		Tsignal_W * outTree;
		Tomography::elec_type DAQtype;
		map<int,Tomography::det_type> det_type_by_asic;
		map<int,int> det_n_by_asic;
		map<Tomography::det_type,unsigned short> det_N;
		long max_event;
		long Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<float> > > > StripAmpl;
		map<Tomography::det_type,vector<vector<float> > > Ped;
		string PedName;
		string RMSName;

};

#endif