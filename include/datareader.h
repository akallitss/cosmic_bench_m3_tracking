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

//helper struct to store information about the detector connected to an ASIC and how it is connected
struct asic_carac{
	Tomography::det_type detector_type;
	int detector_n;
	int asic_n_in_det;
	bool connector_direction;
};

//class to read raw data comming from the electronics and build the signal.root files, the electronic specific part is abstracted in the ElecReader class (cf. ElecReader.h)
class DataReader{
	public:
		DataReader();
		//build object with config tree, is_live is to be used with the -q option of FeuUdpControl
		DataReader(ptree config_tree, bool save_to_disk, bool is_live = false);
		//build object to process an entire folder of fdf file, regrouping them and correcting the evttime with their timestamp
		DataReader(ptree config_tree, string multi_run_dir);
		//DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, map<int,int> feu_id_to_n_, int first_index_, int last_index_, string PedName_, string RMSName_);
		//DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, vector<int> fem_id,int first_index_, int last_index_, string PedName_, string RMSName_);
		~DataReader();
		//read all the fdf files to make a signal.root file
		void process();
		//process a single event
		void process_event();
		//compute the pedestal and output them in the Ped.dat file
		void compute_ped();
		//read the pedestal from the Ped.dat file in order to substract them
		void read_ped();
		//compute the pedestal RMS and output them in the RMSPed.dat file
		void compute_RMSPed();
		//using Ped.dat data, substract the pedestal from all the data in the signal.root file
		void do_ped_sub();
		//substract the common noise from all the pedestal subtracted data in the signal.root file
		void do_common_noise_sub();
		//using Ped.dat data, substract the pedestal from the loaded event data in the signal.root file
		void do_ped_sub_event();
		//substract the common noise from the loaded event pedestal subtracted data in the signal.root file
		void do_common_noise_sub_event();
		//simultaneously substract the pedestal and common noise from the loaded event data in the signal.root file and using the Ped.dat data
		void do_ped_CMN_sub_event();
		//get current event id
		long get_event_n() const;
		//get max event id to be processed
		long get_max_event() const;
		//get current event event time
		double get_evttime() const;
		//check if there is still data to be processed
		bool is_end() const;
		//retrieve current event data : amplitude by detector by channel by sample
		map<Tomography::det_type,vector<vector<vector<float> > > > get_data() const;
		//retrieve the pedestals by detector and by channel
		map<Tomography::det_type,vector<vector<float> > > get_Ped() const;
		//using the ped_in, substract the pedestal and common noise of data_in and output it
		template<typename S, typename T> static map<Tomography::det_type,vector<vector<vector<S> > > > do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<T> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
	protected:
		//void Init(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string PedName_, string RMSName_, string outFileName = "", long max_event_ = -1);
		//call the correct mapping fuction for dream electronics for the given detector
		static int Dream_mapping(Tomography::det_type det,int channel, bool direction);
		//call the correct mapping fuction for feminos electronics for the given detector
		static int Feminos_mapping(Tomography::det_type det,int channel, bool direction);
		//point to either Dream_mapping or Feminos_mapping according to the configuration tree
		int (*mapping)(Tomography::det_type,int,bool);
		ElecReader * reader;
		Tsignal_W * outTree;
		Tomography::elec_type DAQtype;
		//map<int,Tomography::det_type> det_type_by_asic;
		//map<int,int> det_n_by_asic;
		map<int,struct asic_carac> asic_list;
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
