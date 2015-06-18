#ifndef tomography_h
#define tomography_h

#include <ostream>
#include <map>
#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>

using std::ostream;
using std::map;
using std::string;
using std::vector;

using boost::property_tree::ptree;

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_);

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_);

namespace Tomography{

	enum det_type{
		unknown_det,
		CM,
		CM_Demux,
		MG
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
	const double XY_size = 500;
	const double ADC_max = 4096;
	const int Nsample = 32;
	const int SampleMin = 5;
	const int SampleMax = 21;
	const double sigma = 4.;
	const int TOTCut = 3;
	const int Nchannel = 64;
	const int Nasic_FEU = 8;
	const int Nasic_Feminos = 4;
	extern const map<det_type,int> CMN_div;
	const string DreamExt = "fdf";
	const string FeminosExt = "aqs";
	const double chisquare_threshold = 100;
	extern bool live_graphic_display; // toggle updating of canvas during calculation
	extern bool is_batch;
	extern bool can_continue;
	void signal_handler(int s);
	elec_type str_to_elec(string str);
	ostream& operator<<(ostream& os, const det_type& det);
	ostream& operator<<(ostream& os, const strip_type& strip);
	ostream& operator<<(ostream& os, const elec_type& elec);
	//void process_elec_files(ptree config_tree);
	void save_canvases();
}

#endif