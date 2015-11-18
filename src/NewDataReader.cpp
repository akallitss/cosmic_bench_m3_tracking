#include <string>
#include <iostream>
#include <csignal>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "datareader.h"
#include "Tanalyse_W.h"
#include "event.h"
#include "tomography.h"

#include <pthread.h>
#include <queue>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::flush;
using std::setw;

using boost::property_tree::ptree;

using std::queue;

void * reader_thread(void *);
void * ped_thread(void *);
void * multicluster_thread(void *);
void * writer_thread(void *);

struct event_raw_data{
	int Nevent;
	double evttime;
	map<Tomography::det_type,vector<vector<vector<float> > > > strip_data;
};
struct event_ped_data{
	int Nevent;
	double evttime;
	map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
};
struct event_object_data{
	map<Tomography::det_type,vector<Event*> > event_data;
	int Nevent;
	double evttime;
};
bool reader_active = true;
bool ped_active = true;
bool multicluster_active = true;
bool writer_active = true;
int event_read = 0;
int event_corr = 0;
int event_demux = 0;
int event_written = 0;
pthread_attr_t reader_attr;
DataReader * blah;
int total_det;
map<Tomography::det_type,vector<vector<float> > > current_ped;
queue<event_raw_data> event_raw_data_queue;
pthread_mutex_t event_raw_data_mutex;
queue<event_ped_data> event_ped_data_queue;
pthread_mutex_t event_ped_data_mutex;
queue<event_object_data> event_objects_queue;
pthread_mutex_t event_objects_mutex;
CosmicBench * bench;
Tanalyse_W * analysisFile;
void * reader_thread(void *){
	reader_active = true;
	while((!(blah->is_end())) && Tomography::get_instance()->get_can_continue()){
		blah->process_event();
		struct event_raw_data current_data;
		current_data.Nevent = blah->get_event_n();
		current_data.evttime = blah->get_evttime();
		current_data.strip_data = blah->get_data();
		pthread_mutex_lock(&event_raw_data_mutex);
		event_raw_data_queue.push(current_data);
		pthread_mutex_unlock(&event_raw_data_mutex);
		event_read++;
		while(event_raw_data_queue.size()>5000){

		}
	}
	reader_active = false;
}
void * ped_thread(void *){
	ped_active = true;
	while(reader_active && Tomography::get_instance()->get_can_continue()){
		while(!event_raw_data_queue.empty()){
			struct event_ped_data current_ped_data;
			struct event_raw_data current_raw_data;
			pthread_mutex_lock(&event_raw_data_mutex);
			current_raw_data = event_raw_data_queue.front();
			event_raw_data_queue.pop();
			pthread_mutex_unlock(&event_raw_data_mutex);
			current_ped_data.Nevent = current_raw_data.Nevent;
			current_ped_data.evttime = current_raw_data.evttime;
			current_ped_data.strip_data = DataReader::do_ped_CMN_sub_event<double,float>(current_raw_data.strip_data,current_ped);
			pthread_mutex_lock(&event_ped_data_mutex);
			event_ped_data_queue.push(current_ped_data);
			pthread_mutex_unlock(&event_ped_data_mutex);
			event_corr++;
			while(event_ped_data_queue.size()>5000){
				
			}
		}
	}
	ped_active = false;
}
void * multicluster_thread(void *){
	multicluster_active = true;
	while(ped_active && Tomography::get_instance()->get_can_continue()){
		while(!event_ped_data_queue.empty()){
			struct event_ped_data current_ped_data;
			pthread_mutex_lock(&event_ped_data_mutex);
			current_ped_data = event_ped_data_queue.front();
			event_ped_data_queue.pop();
			pthread_mutex_unlock(&event_ped_data_mutex);
			struct event_object_data current_object_data;
			current_object_data.Nevent = current_ped_data.Nevent;
			current_object_data.evttime = current_ped_data.evttime;
			for(int i=0;i<total_det;i++){
				Detector * det = bench->get_detector(i);
				current_object_data.event_data[det->get_type()].push_back(det->build_event(current_ped_data.strip_data[det->get_type()][det->get_n_in_tree()],current_ped_data.Nevent));
				(current_object_data.event_data[det->get_type()].back())->MultiCluster();
			}
			pthread_mutex_lock(&event_objects_mutex);
			event_objects_queue.push(current_object_data);
			pthread_mutex_unlock(&event_objects_mutex);
			event_demux++;
			while(event_objects_queue.size()>5000){
				
			}
		}
	}
	multicluster_active = false;
}
void * writer_thread(void *){
	writer_active = true;
	while(multicluster_active && Tomography::get_instance()->get_can_continue()){
		while(!event_objects_queue.empty()){
			struct event_object_data current_object_data;
			pthread_mutex_lock(&event_objects_mutex);
			current_object_data = event_objects_queue.front();
			event_objects_queue.pop();
			pthread_mutex_unlock(&event_objects_mutex);
			analysisFile->fillTree(current_object_data.Nevent,current_object_data.evttime,current_object_data.event_data);
			for(map<Tomography::det_type,vector<Event*> >::iterator type_it = current_object_data.event_data.begin();type_it!=current_object_data.event_data.end();++type_it){
				for(vector<Event*>::iterator ev_it = (type_it->second).begin();ev_it!=(type_it->second).end();++ev_it){
					delete *ev_it;
				}
			}
			event_written++;
		}
	}
	writer_active = false;
}

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the Run caracs" << endl;
		return 1;
	}
	enum process_type{
		data_run,
		ped_run,
		analysis
	};
	process_type operation = data_run;
	for(int i=2;i<argc;i++){
		if(argv[i] == string("ped")){
			operation = ped_run;
		}
		else if(argv[i] == string("data")){
			operation = data_run;
		}
		else if(argv[i] == string("analyse")){
			operation = analysis;
		}
	}
	string config_file = argv[1];
	ptree config_tree;
	read_json(config_file, config_tree);
	Tomography::Init(config_tree);
	if(operation != analysis){
		blah = new DataReader(config_tree,true);
		blah->process();
		if(operation == ped_run) blah->compute_ped();
		blah->read_ped();
		blah->do_ped_sub();
		blah->do_common_noise_sub();
		if(operation == ped_run) blah->compute_RMSPed();
		delete blah;
	}
	else{
		bench = new CosmicBench(config_tree);
		analysisFile = new Tanalyse_W(config_tree.get<string>("Tree"),bench->get_det_N());
		blah = new DataReader(config_tree,false);
		blah->read_ped();
		current_ped = blah->get_Ped();
		total_det = bench->get_det_N_tot();
		//long event_nb = 0;
		//int Nevent = 0;
		//double evttime = 0;
		pthread_mutex_init(&event_raw_data_mutex, NULL);
		pthread_mutex_init(&event_ped_data_mutex, NULL);
		pthread_mutex_init(&event_objects_mutex, NULL);
		
		pthread_attr_init(&reader_attr);
		pthread_attr_setdetachstate(&reader_attr, PTHREAD_CREATE_JOINABLE);
		pthread_t reader_id;
		pthread_t ped_id;
		pthread_t multicluster_id;
		pthread_t writer_id;
		
		bool thread_launched = true;
		int result = pthread_create(&reader_id,NULL, reader_thread, NULL);
		if(result == 0 && thread_launched){
			result = pthread_create(&ped_id,NULL, ped_thread, NULL);
		}
		else{
			thread_launched = false;
			cout << "cannot create reader thread" << endl;
		}
		if(result == 0 && thread_launched){
			result = pthread_create(&multicluster_id,NULL, multicluster_thread, NULL);
		}
		else{
			thread_launched = false;
			cout << "cannot create ped correcter thread" << endl;
		}
		if(result == 0 && thread_launched){
			result = pthread_create(&writer_id,NULL, writer_thread, NULL);
		}
		else{
			thread_launched = false;
			cout << "cannot create multicluster thread" << endl;
		}
		if(result !=0 && thread_launched){
			thread_launched = false;
			cout << "cannot create writer thread" << endl;
		}
		
		pthread_attr_destroy(&reader_attr);
		if(thread_launched){
			cout << setw(15) << "event read" << "|" << setw(15) << "queue for corr" << "|" << setw(15) << "event corr" << "|" << setw(15) << "queue for demux" << "|" << setw(15) << "event demux" << "|" << setw(15) << "queue for write" << "|" << setw(15) << "event written" << endl;
			while(reader_active || ped_active || multicluster_active || writer_active){
				cout << "\r" << setw(15) << event_read << "|" << setw(15) << event_corr << "|" << setw(15) << event_demux << "|" << setw(15) << event_written << flush;
				usleep(100000);
			}
			cout << "\r" << setw(15) << event_read << "|" << setw(15) << event_corr << "|" << setw(15) << event_demux << "|" << setw(15) << event_written << endl;
			void * status;
			result = pthread_join(reader_id,&status);
			if(result !=0){
				cout << "cannot join reader thread" << endl;
			}
			result = pthread_join(ped_id,&status);
			if(result !=0){
				cout << "cannot join ped correcter thread" << endl;
			}
			result = pthread_join(multicluster_id,&status);
			if(result !=0){
				cout << "cannot join multicluster thread" << endl;
			}
			result = pthread_join(writer_id,&status);
			if(result !=0){
				cout << "cannot join writer thread" << endl;
			}
		}
		analysisFile->Write();
		analysisFile->CloseFile();
		delete blah;
		delete analysisFile;
		delete bench;
	}
	Tomography::Quit();
	return 0;
}