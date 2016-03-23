#include <string>
#include <iostream>
#include <csignal>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "datareader.h"
#include "Tanalyse_W.h"
#include "Tsignal_W.h"
#include "event.h"
#include "tomography.h"

#include "MT_tomography.h"
#include "task/ped_task.h"
#include "task/multicluster_task.h"
#include "task/write_analyse_task.h"
#include "task/write_signal_task.h"
#include "task/read_elec_task.h"

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
		live,
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
	}
	string config_file = argv[1];
	ptree config_tree;
	read_json(config_file, config_tree);
	Tomography::Init(config_tree);
	if((operation == live) || (operation == read)){
		DataReader * blah = new DataReader(config_tree,true,operation==live);
		blah->process();
		delete blah;
	}
	else if(operation == analysis){
		CosmicBench * bench = new CosmicBench(config_tree);
		Tanalyse_W * analysisFile = new Tanalyse_W(config_tree.get<string>("Tree"),bench->get_det_N());
		DataReader * blah = new DataReader(config_tree,false);
		blah->read_ped();
		map<Tomography::det_type,vector<vector<float> > > current_ped = blah->get_Ped();

		Output_Task<event_data> * to_write = new Write_Analyse_Task(analysisFile);
		Input_Task * to_do = new Read_Elec_Task(blah, new Ped_Corr_Task(current_ped, new Multicluster_Task(bench,to_write)));
		vector<Thread*> threads;
		threads.push_back(new Writer_Thread(to_write));
		(threads.back())->start();
		threads.push_back(new Reader_Thread(to_do));
		(threads.back())->start();
		unsigned short n_thread = Tomography::get_instance()->get_thread_number() - threads.size();
		if(n_thread<1) n_thread = 1;
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
	else if(operation == pyramids){
		CosmicBench * bench = new CosmicBench(config_tree);
		Tanalyse_W * analysisFile = new Tanalyse_W(config_tree.get<string>("Tree"),bench->get_det_N());
		DataReader * blah = new DataReader(config_tree,false);
		blah->read_ped();
		map<Tomography::det_type,vector<vector<float> > > current_ped = blah->get_Ped();
		Tsignal_W * signalFile = new Tsignal_W(config_tree.get<string>("signal_file"),bench->get_det_N());
		Output_Task<event_data> * to_write_analyse = new Write_Analyse_Task(analysisFile);
		Output_Task<raw_data> * to_write_signal = new Write_Signal_Task<raw_data>(signalFile,new Ped_Corr_Task(current_ped,new Multicluster_Task(bench,to_write_analyse)));
		Input_Task * to_do = new Read_Elec_Task(blah,to_write_signal);
		vector<Thread*> threads;
		threads.push_back(new Reader_Thread(to_do));
		(threads.back())->start();
		threads.push_back(new Writer_Thread(to_write_signal));
		(threads.back())->start();
		threads.push_back(new Writer_Thread(to_write_analyse));
		(threads.back())->start();
		unsigned short n_thread = Tomography::get_instance()->get_thread_number() - threads.size();
		if(n_thread<1) n_thread = 1;
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
		delete signalFile;
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
