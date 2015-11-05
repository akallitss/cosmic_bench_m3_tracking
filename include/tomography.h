#ifndef tomography_h
#define tomography_h

#include <ostream>
#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <boost/property_tree/ptree.hpp>

using std::cout;
using std::endl;
using std::ostream;
using std::map;
using std::string;
using std::vector;

using boost::property_tree::ptree;

class Detector;

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_);

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_);

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
	const int Nsample = 32;
	const int Nchannel = 64;
	const int Nasic_FEU = 8;
	const int Nasic_Feminos = 4;
	const string DreamExt = "fdf";
	const string FeminosExt = "aqs";

	//Run dependant
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

#endif