#ifndef livedisplay_h
#define livedisplay_h

#include <vector>
#include <string>
#include <sys/inotify.h>

using std::vector;
using std::string;

class liveDisplay{
	public:
		liveDisplay();
		~liveDisplay();
		liveDisplay(int max_event_);
		int start_inotify(string filename);
		int pause_inotify();
		int resume_inotify();
		int stop_inotify();
		unsigned int read_inotify();
		void add_file(string filename);
		void flux_map(double z);
	protected:
		vector<string> filenames;
		int max_event;
		int inotify_descriptor;
		int file_descriptor;
		bool inotify_started;
		string current_file;
};

#endif