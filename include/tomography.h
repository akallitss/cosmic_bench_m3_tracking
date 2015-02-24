#ifndef tomography_h
#define tomography_h

#include <ostream>
#include <map>
#include <string>

using std::ostream;
using std::map;
using std::string;

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_);

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
	const double chisquare_threshold = 100;
	const bool live_graphic_display = true; // toggle updating of canvas during calculation
	elec_type str_to_elec(string str);
	ostream& operator<<(ostream& os, const det_type& det);
	ostream& operator<<(ostream& os, const strip_type& strip);
	ostream& operator<<(ostream& os, const elec_type& elec);
}

#endif