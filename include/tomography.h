#ifndef tomography_h
#define tomography_h

#include <ostream>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <csignal>

#include <queue>
#include <pthread.h>

#include <TRint.h>

#include <boost/property_tree/ptree.hpp>

#include "ray.h"
//#include "MT_tomography.h"

using std::cout;
using std::endl;
using std::ostream;
using std::map;
using std::string;
using std::vector;

using std::queue;

using boost::property_tree::ptree;

class Detector;
class CosmicBenchEvent;
class Event;

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_);

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_);

class Tomography{
	public:

		enum det_type{
			unknown_det,
			CM,
			MG,
			MGv2
		};
		enum strip_type{
			unknown_strip = 10,
			Wide,
			Thin,
			Demux
		};
		enum elec_type{
			unknown_elec = 100,
			Dream,
			Feminos
		};
		enum signal_type{
			unknown_signal = 1000,
			raw,
			ped,
			corr

		};
		static elec_type str_to_elec(string str);

		static constexpr const double ADC_max = 4096;
		static constexpr const int Nchannel = 64;
		static constexpr const int Nasic_FEU = 8;
		static constexpr const int Nasic_Feminos = 4;
		static constexpr const int Max_Nsample = 512;
		static const string DreamExt;
		static const string FeminosExt;
		//TODO : make this map const
		static map<const Tomography::det_type,const Detector* const> Static_Detector;

		static Tomography* get_instance();
		static void Init(ptree config_tree_);
		static void Init(string config_tree_file);
		static void Quit();

		int get_Nsample() const;
		double get_XY_size() const;
		int get_SampleMin() const;
		int get_SampleMax() const;
		int get_TOTCut() const;
		bool get_is_up(int layer) const;
		double get_sigma() const;
		double get_noise_RMS() const;
		double get_chisquare_threshold() const;
		double get_clock_rate() const;
		bool get_live_graphic_display() const;
		bool get_is_batch() const;
		bool get_can_continue() const;

		unsigned short get_thread_number() const;

		void save_canvases();
		void Run();

	private:
		static Tomography* singleton_instance;
		Tomography();
		Tomography(ptree config_tree_);
		~Tomography();

		TRint * root_interpreter;

		ptree config_tree;
		int Nsample;
		double XY_size;
		int SampleMin;
		int SampleMax;
		double sigma;
		int TOTCut;
		int first_down_layer;
		double clock_rate;
		double noise_RMS;
		double chisquare_threshold;
		bool live_graphic_display; // toggle updating of canvas during calculation
		bool is_batch;
		unsigned short thread_number;
		unsigned long raw_data_treated;
		unsigned long ped_data_treated;
		unsigned long corr_data_treated;
		unsigned long event_treated;
		unsigned long ray_treated;
		unsigned long deviation_treated;

		struct sigaction sigIntHandler;
		static void signal_handler(int s);
		bool can_continue;


};
ostream& operator<<(ostream& os, const Tomography::det_type& det);
ostream& operator<<(ostream& os, const Tomography::strip_type& strip);
ostream& operator<<(ostream& os, const Tomography::elec_type& elec);

/*
namespace Tomography{

	enum det_type{
		unknown_det,
		CM,
		MG,
		MGv2
	};
	enum strip_type{
		unknown_strip = 10,
		Wide,
		Thin,
		Demux
	};
	enum elec_type{
		unknown_elec = 100,
		Dream,
		Feminos
	};

	//Utility function for these enums
	elec_type str_to_elec(string str);
	ostream& operator<<(ostream& os, const det_type& det);
	ostream& operator<<(ostream& os, const strip_type& strip);
	ostream& operator<<(ostream& os, const elec_type& elec);

	//Real const
	const double ADC_max = 4096;
	const int Nchannel = 64;
	const int Nasic_FEU = 8;
	const int Nasic_Feminos = 4;
	const string DreamExt = "fdf";
	const string FeminosExt = "aqs";

	//Run dependant
	const int Nsample = 32;
	const double XY_size = 500;
	const int SampleMin = 5;
	const int SampleMax = 21;
	const double sigma = 4.;
	const int TOTCut = 3;
	const double chisquare_threshold = 100;

	//List of dummy Detector to iterate over and/or get static members
	extern map<const Tomography::det_type,const Detector* const> Static_Detector;

	//Help to handle batch/graphic output
	extern bool live_graphic_display; // toggle updating of canvas during calculation
	extern bool is_batch;
	void save_canvases();

	//Handle Ctrl-C
	extern bool can_continue;// = true;
	inline void signal_handler(int s){
		cout << "\nCaught signal " << s << endl;
		cout << endl;
		can_continue = false;
	}

	//void process_elec_files(ptree config_tree);
}
*/
#endif
