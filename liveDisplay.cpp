#define livedisplay_cpp

#include "liveDisplay.h"

#include <string>
#include <sys/inotify.h>
#include <unistd.h>

#include <TH2D.h>

using std::string;

//inotify stuff
#define MAX_EVENTS 1024
#define LEN_NAME 16
#define EVENT_SIZE  ( sizeof (struct inotify_event) )
#define BUF_LEN     ( MAX_EVENTS * ( EVENT_SIZE + LEN_NAME ))

liveDisplay::liveDisplay(){
	filenames.clear();
	max_event = -1;
	inotify_descriptor = -1;
	file_descriptor = -1;
	inotify_started = false;
	current_file = "";
}

liveDisplay::~liveDisplay(){

}

liveDisplay::liveDisplay(int max_event_){
	filenames.clear();
	max_event = max_event_;
}

void liveDisplay::add_file(string filename){
	filenames.push_back(filename);
}

int liveDisplay::start_inotify(string filename){
	inotify_descriptor = inotify_init();
	current_file = filename;
	file_descriptor = inotify_add_watch(inotify_descriptor,current_file.c_str(), IN_ALL_EVENTS);
	if(inotify_descriptor>=0 && file_descriptor>=0) inotify_started = true;
	if(inotify_descriptor<0) return inotify_descriptor;
	return file_descriptor;
}
int liveDisplay::pause_inotify(){
	if(!inotify_started) return -1;
	inotify_rm_watch(inotify_descriptor,file_descriptor);
	close(inotify_descriptor);
	return 1;
}
int liveDisplay::resume_inotify(){
	if(!inotify_started) return -1;
	inotify_descriptor = inotify_init();
	file_descriptor = inotify_add_watch(inotify_descriptor,current_file.c_str(), IN_ALL_EVENTS);
	if(inotify_descriptor<0) return inotify_descriptor;
	return file_descriptor;
}
int liveDisplay::stop_inotify(){
	if(!inotify_started) return -1;
	inotify_rm_watch(inotify_descriptor,file_descriptor);
	close(inotify_descriptor);
	inotify_started = false;
	return 1;
}

unsigned int liveDisplay::read_inotify(){
	unsigned int global_mask = 0x00000000;
	if(!inotify_started) return global_mask;
	char buffer[BUF_LEN];
	int length = read(inotify_descriptor,buffer, BUF_LEN);
	int i = 0;
	while(i<length){
		struct inotify_event * event = (struct inotify_event*) &buffer[i];
		global_mask |= event->mask;
		i+= EVENT_SIZE + event->len;
	}
	return global_mask;
}

void liveDisplay::flux_map(double z){

}