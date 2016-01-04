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
#include "MT_tomography.h"

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
class Display_Thread;

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_);

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_);

struct raw_data;
struct ped_data;
struct event_data;
struct ray_data;
struct deviation_data;

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

		static const double ADC_max = 4096;
		static const int Nchannel = 64;
		static const int Nasic_FEU = 8;
		static const int Nasic_Feminos = 4;
		static const int Max_Nsample = 512;
		static const string DreamExt;
		static const string FeminosExt;
		//TODO : make this map const
		static map<const Tomography::det_type,const Detector* const> Static_Detector;
		static Display_Thread TomoCout;

		static Tomography* get_instance();
		static void Init(ptree config_tree_);
		static void Init(string config_tree_file);
		static void Quit();

		int get_Nsample() const;
		double get_XY_size() const;
		int get_SampleMin() const;
		int get_SampleMax() const;
		int get_TOTCut() const;
		double get_sigma() const;
		double get_chisquare_threshold() const;
		bool get_live_graphic_display() const;
		bool get_is_batch() const;
		bool get_can_continue() const;

		struct raw_data get_next_raw_data();
		void push_next_raw_data(struct raw_data new_data);
		bool is_raw_data_empty() const;
		struct ped_data get_next_ped_data();
		void push_next_ped_data(struct ped_data new_data);
		bool is_ped_data_empty() const;
		struct ped_data get_next_corr_data();
		void push_next_corr_data(struct ped_data new_data);
		bool is_corr_data_empty() const;
		struct event_data get_next_event_data();
		void push_next_event_data(struct event_data new_data);
		bool is_event_data_empty() const;
		struct ray_data get_next_ray_data();
		void push_next_ray_data(struct ray_data new_data);
		bool is_ray_data_empty() const;
		struct deviation_data get_next_deviation_data();
		void push_next_deviation_data(struct deviation_data new_data);
		bool is_deviation_data_empty() const;

		unsigned short get_thread_number() const;

		string init_count() const;
		string print_count() const;

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
		double chisquare_threshold;
		bool live_graphic_display; // toggle updating of canvas during calculation
		bool is_batch;
		unsigned short thread_number;
		queue<raw_data> raw_data_queue;
		pthread_mutex_t raw_data_mutex;
		unsigned long raw_data_treated;
		queue<ped_data> ped_data_queue;
		pthread_mutex_t ped_data_mutex;
		unsigned long ped_data_treated;
		queue<ped_data> corr_data_queue;
		pthread_mutex_t corr_data_mutex;
		unsigned long corr_data_treated;
		queue<event_data> event_queue;
		pthread_mutex_t event_mutex;
		unsigned long event_treated;
		queue<ray_data> ray_queue;
		pthread_mutex_t ray_mutex;
		unsigned long ray_treated;
		queue<deviation_data> deviation_queue;
		pthread_mutex_t deviation_mutex;
		unsigned long deviation_treated;

		struct sigaction sigIntHandler;
		static void signal_handler(int s);
		bool can_continue;


};
ostream& operator<<(ostream& os, const Tomography::det_type& det);
ostream& operator<<(ostream& os, const Tomography::strip_type& strip);
ostream& operator<<(ostream& os, const Tomography::elec_type& elec);

struct raw_data{
	int Nevent;
	double evttime;
	map<Tomography::det_type,vector<vector<vector<float> > > > strip_data;
	raw_data(): Nevent(-1), evttime(0), strip_data(){}
};
struct ped_data{
	int Nevent;
	double evttime;
	map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
	ped_data(): Nevent(-1), evttime(0), strip_data(){}
};
struct event_data{
	map<Tomography::det_type,vector<Event*> > det_data;
	int Nevent;
	double evttime;
	event_data(): det_data(), Nevent(-1), evttime(0){}
};
struct ray_data{
	CosmicBenchEvent * CBevent;
	vector<Ray> rays;
	ray_data(): CBevent(NULL), rays(){}
};
struct deviation_data{
	CosmicBenchEvent * CBevent;
	vector<RayPair> rays;
	deviation_data(): CBevent(NULL), rays(){}
};

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