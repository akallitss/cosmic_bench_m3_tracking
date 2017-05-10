#include <string>
#include <iostream>
#include <csignal>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "datareader.h"
#include "Tanalyse_W.h"
#include "Tsignal_W.h"
#include "Tray.h"
#include "event.h"
#include "tomography.h"

#include "MT_tomography.h"
#include "task/ped_task.h"
#include "task/multicluster_task.h"
#include "task/write_analyse_task.h"
#include "task/write_signal_task.h"
#include "task/write_rays_task.h"
#include "task/read_elec_task.h"
#include "task/tracking_task.h"

//#include <pthread.h>
//#include <queue>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::flush;
using std::setw;

using boost::property_tree::ptree;

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the Run caracs" << endl;
		return 1;
	}
	enum process_type{
		data_run,
		ped_run,
		analysis,
		pyramids,
		pyrarays,
		mcube,
		live,
		watto,
		read
	};
	process_type operation = read;
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
		else if(argv[i] == string("live")){
			operation = live;
		}
		else if(argv[i] == string("read")){
			operation = read;
		}
		else if(argv[i] == string("pyramids")){
			operation = pyramids;
		}
		else if(argv[i] == string("pyrarays")){
			operation = pyrarays;
		}
		else if(argv[i] == string("mcube")){
			operation = mcube;
		}
		else if(argv[i] == string("watto")){
			operation = watto;
		}
	}
	string config_file = argv[1];
	ptree config_tree;
	read_json(config_file, config_tree);
	Tomography::Init(config_tree);
	if(operation == read){
		DataReader * blah = new DataReader(config_tree,true,false);
		blah->process();
		delete blah;
	}
	else if((operation == analysis) || (operation == watto)){
		CosmicBench * bench = new CosmicBench(config_tree);
		Tanalyse_W * analysisFile = new Tanalyse_W(config_tree.get<string>("Tree"),bench->get_det_N());
		DataReader * blah = NULL;
		if(operation==analysis) blah = new DataReader(config_tree,false);
		else blah = new DataReader(config_tree,"./");
		blah->read_ped();
		map<Tomography::det_type,vector<vector<float> > > current_ped = blah->get_Ped();

		Output_Task<event_data> * to_write = new Write_Analyse_Task(analysisFile);
		Input_Task * to_do = new Read_Elec_Task(blah, new Ped_Corr_Task(current_ped, new Multicluster_Task(bench,to_write)));
		vector<Thread*> threads;
		threads.push_back(new Writer_Thread(to_write));
		(threads.back())->start();
		threads.push_back(new Reader_Thread(to_do));
		(threads.back())->start();
		const unsigned short n_thread = (Tomography::get_instance()->get_thread_number() > threads.size()) ? (Tomography::get_instance()->get_thread_number() - threads.size()) : 1;
		cout << "1 | " << n_thread << " | 1" << endl;
		for(unsigned short i=0;i<n_thread;i++){
			threads.push_back(new Worker_Thread());
			(threads.back())->start();
		}
		Display_Thread::get_instance()->start_count();
		//cout << Tomography::get_instance()->init_count() << "|" << setw(7) << "tasks" << endl;
		bool has_working_thread = true;
		while(has_working_thread && Tomography::get_instance()->get_can_continue()){
			//cout << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << Task::task_left() << flush;
			has_working_thread = false;
			for(unsigned short i=0;i<threads.size();i++){
				if(threads[i]->is_working()){
					has_working_thread = true;
					break;
				}
			}
			usleep(10000);
		}
		for(unsigned short i=0;i<threads.size();i++){
			threads[i]->stop();
			delete threads[i];
		}
		Display_Thread::get_instance()->stop_count();
		//cout << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << Task::task_left() << endl;
		analysisFile->Write();
		analysisFile->CloseFile();
		delete blah;
		delete analysisFile;
		delete bench;
	}
	else if((operation == pyramids) || (operation == pyrarays) || (operation == live) || (operation == mcube)){
		CosmicBench * bench = new CosmicBench(config_tree);
		Tsignal_W * signalFile = NULL;
		Tanalyse_W * analysisFile = NULL;
		Tray * raysFile = NULL;
		DataReader * blah = new DataReader(config_tree,false,true);
		blah->read_ped();
		map<Tomography::det_type,vector<vector<float> > > current_ped = blah->get_Ped();

		Output_Task<event_data> * to_write_analyse = NULL;
		Output_Task<ray_data> * to_write_rays = NULL;
		vector<Thread*> threads;
		Typed_Task<event_data> * follow_up_analyse_task = NULL;
		if(operation==pyrarays){
			cout << "creating rays output chain" << endl;
			raysFile = new Tray(config_tree.get<string>("metadata") + "_rays.root");
			to_write_rays = new Write_Rays_Task(raysFile,bench->get_z_Up(),bench->get_z_Down(),config_tree.get<string>("metadata") + "_ampl.txt");
			threads.push_back(new Writer_Thread(to_write_rays));
			follow_up_analyse_task = new Tracking_Abs_Task(bench,to_write_rays);
		}
		Output_Task<raw_data> * to_write_signal = NULL;
		Typed_Task<raw_data> * follow_up_signal_task = NULL;
		if(operation!=live){
			cout << "creating analyse output chain" << endl;
			analysisFile = new Tanalyse_W(config_tree.get<string>("Tree"),bench->get_det_N());
			to_write_analyse = new Write_Analyse_Task(analysisFile, follow_up_analyse_task);
			threads.push_back(new Writer_Thread(to_write_analyse));
			follow_up_signal_task = new Ped_Corr_Task(current_ped,new Multicluster_Task(bench,to_write_analyse));
		}
		Input_Task * to_do = NULL;
		if(operation!=mcube){
			cout << "creating signal output chain" << endl;
			signalFile = new Tsignal_W(config_tree.get<string>("signal_file"),bench->get_det_N());
			to_write_signal = new Write_Signal_Task<raw_data>(signalFile,follow_up_signal_task);
			threads.push_back(new Writer_Thread(to_write_signal));
			to_do = new Read_Elec_Task(blah,to_write_signal);
		}
		else{
			cout << "skipping signal file" << endl;
			to_do = new Read_Elec_Task(blah,follow_up_signal_task);
		}
		threads.push_back(new Reader_Thread(to_do));
		for(vector<Thread*>::reverse_iterator thread_it=threads.rbegin();thread_it!=threads.rend();++thread_it){
			(*thread_it)->start();
		}

		const unsigned short n_thread = (Tomography::get_instance()->get_thread_number() > threads.size()) ? (Tomography::get_instance()->get_thread_number() - threads.size()) : 1;
		cout << " reader : 1 | worker : " << n_thread << " | writer : " << (threads.size() - 1) << endl;
		for(unsigned short i=0;i<n_thread;i++){
			threads.push_back(new Worker_Thread());
			(threads.back())->start();
		}
		Display_Thread::get_instance()->start_count();
		//cout << Tomography::get_instance()->init_count() << "|" << setw(7) << "tasks" << endl;
		bool has_working_thread = true;
		while(has_working_thread && Tomography::get_instance()->get_can_continue()){
			//cout << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << Task::task_left() << flush;
			has_working_thread = false;
			for(unsigned short i=0;i<threads.size();i++){
				if(threads[i]->is_working()){
					has_working_thread = true;
					break;
				}
			}
			usleep(10000);
		}
		for(unsigned short i=0;i<threads.size();i++){
			threads[i]->stop();
			delete threads[i];
		}
		Display_Thread::get_instance()->stop_count();
		//cout << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << Task::task_left() << endl;
		if(operation!=mcube){
			signalFile->Write();
			signalFile->CloseFile();
			delete signalFile;
		}
		if(operation != live){
			analysisFile->Write();
			analysisFile->CloseFile();
			delete analysisFile;
			if(operation==pyrarays){
				raysFile->Write();
				raysFile->CloseFile();
				delete raysFile;
			}
		}
		delete blah;
		delete bench;
		delete to_do;
	}
	else{
		DataReader * blah = new DataReader(config_tree,true);
		if(operation == ped_run) blah->compute_ped();
		blah->read_ped();
		blah->do_ped_sub();
		blah->do_common_noise_sub();
		if(operation == ped_run) blah->compute_RMSPed();
		delete blah;
	}
	Tomography::Quit();
	return 0;
}
