#ifndef livedisplay_h
#define livedisplay_h

#include "detector.h"
#include "tomography.h"

#include <vector>
#include <map>
#include <string>
#include <sys/inotify.h>

using std::vector;
using std::string;
using std::map;

//depreacted, currently broken
//used before FeuUdpControl -q option was available
class liveDisplay: public CosmicBench{
	public:
		liveDisplay();
		~liveDisplay();
		liveDisplay(string config_file);
		int start_inotify(string filename);
		int pause_inotify();
		int resume_inotify();
		int stop_inotify();
		unsigned int read_inotify();
		void add_file(string filename);
		void add_files(int first,int last);
		void flux_map(double z);
	protected:
		Tomography::elec_type electronic_type;
		vector<string> filenames;
		long max_event;
		int inotify_descriptor;
		int file_descriptor;
		bool inotify_started;
		string current_file;
		map<int,Tomography::det_type> det_type_by_asic;
		map<int,int> det_n_by_asic;
		bool use_srf;
		string data_file_basename;
		long max_file_size;
		map<Tomography::det_type,vector<vector<float> > > Pedestal;
		//float (*Pedestal_MG)[61];
		//float (*Pedestal_CM)[64];
};

#endif