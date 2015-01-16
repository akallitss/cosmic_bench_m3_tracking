#ifndef tomography_h
#define tomography_h

#include <ostream>
#include <map>
#include <string>

using std::ostream;
using std::map;
using std::string;

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
	const int Nsample = 32;
	elec_type str_to_elec(string str);
	ostream& operator<<(ostream& os, const det_type& det);
	ostream& operator<<(ostream& os, const strip_type& strip);
	ostream& operator<<(ostream& os, const elec_type& elec);
}

#endif