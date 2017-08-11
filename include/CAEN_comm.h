#ifndef CAEN_comm_h
#define CAEN_comm_h
#include <string>
#include <map>
#include <vector>

using std::string;
using std::map;
using std::vector;

//wrapper around the CAEN library
namespace CAEN_Ch{

	//get/settable parameter for CAEN crates
	enum Param{
		V0Set = 0,
		I0Set,
		V1Set,
		I1Set,
		RUp,
		RDown,
		Trip,
		SVMax,
		VMon,
		IMon,
		Status = 100,
		Pw,
		POn,
		PDwn,
		TripInt,
		TripExt
	};

	//string to enum conversion
	string get_Param_name(Param param);
	Param get_Param(string name);
	//some parameters are linked to float data and others to integer data
	bool Param_is_float(Param param);
	//parameter representing a voltage can be negative
	bool Param_is_voltage(Param param);

	typedef union param_value param_value;
	union param_value{
		int discrete;
		float real;
	};

}

//store specific board (part of a crate) information
class Board_spec{
	public:
		Board_spec();
		//build object with the board index, the number of HV channel it carries and if these channels are positive ones
		Board_spec(int index_,int chan_nb_,bool is_pos_);
		~Board_spec();
		int get_index() const;
		int get_chan_nb() const;
		bool get_is_pos() const;
	protected:
		int index;
		int chan_nb;
		int is_pos;
};
//allow to communicate with the CAEN crate
class CAEN_Comm{
	public:
		CAEN_Comm();
		//build object by connecting to the crate at the given IP with the given credentials
		CAEN_Comm(string IP, string username, string passwd);
		~CAEN_Comm();
		//retrieve channel names from each channel of each board of the crate
		map<int,map<int,string> > get_Ch_name();
		//retrieve channel given parameter from each channel of each board of the crate
		map<int,map<int,CAEN_Ch::param_value> > get_Ch_param(map<int,vector<int> > Ch_list, CAEN_Ch::Param param);
	protected:
		int handle;
		map<int,Board_spec> boards;
};

#endif