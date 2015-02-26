#include <string>
#include <map>
#include <vector>

using std::string;
using std::map;
using std::vector;

namespace CAEN_Ch{

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

	string get_Param_name(Param param);
	Param get_Param(string name);
	bool Param_is_float(Param param);
	bool Param_is_voltage(Param param);

	typedef union param_value param_value;
	union param_value{
		int discrete;
		float real;
	};

}

class Board_spec{
	public:
		Board_spec();
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

class CAEN_Comm{
	public:
		CAEN_Comm();
		CAEN_Comm(string IP, string username, string passwd);
		~CAEN_Comm();
		map<int,map<int,string> > get_Ch_name();
		map<int,map<int,CAEN_Ch::param_value> > get_Ch_param(map<int,vector<int> > Ch_list, CAEN_Ch::Param param);
	protected:
		int handle;
		map<int,Board_spec> boards;
};