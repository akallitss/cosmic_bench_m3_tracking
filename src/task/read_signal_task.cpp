#define read_signal_task_cpp

#include "task/read_signal_task.h"

#include "Tsignal_R.h"

Read_Signal_Task::Read_Signal_Task(long max_event_, Tsignal_R * reader_, Tomography::signal_type type_): Input_Task(max_event_){
	reader = reader_;
	type = type_;
}
Read_Signal_Task::Read_Signal_Task(long max_event_, Tsignal_R * reader_, Tomography::signal_type type_,Task * next_task_): Input_Task(max_event_, next_task_){
	reader = reader_;
	type = type_;
}
Read_Signal_Task::~Read_Signal_Task(){

}
bool Read_Signal_Task::do_task(){
	if(type == Tomography::raw){
		struct raw_data current_data;
		pthread_mutex_lock(&IO_mutex);
		bool has_read = reader->GetNext();
		if(has_read){
			current_data.Nevent = reader->Nevent;
			current_data.evttime = reader->evttime;
			for(map<Tomography::det_type,unsigned short>::iterator type_it=(reader->det_N).begin();type_it!=(reader->det_N).end();++type_it){
				current_data.strip_data[type_it->first] = vector<vector<vector<float> > >(type_it->second,vector<vector<float> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<float>(Tomography::get_instance()->get_Nsample(),0)));
				for(unsigned short i=0;i<type_it->second;i++){
					current_data.strip_data[type_it->first][i] = reader->get_ampl_raw<float>(type_it->first,i);
				}
			}
			pthread_mutex_unlock(&IO_mutex);
		}
		else pthread_mutex_unlock(&IO_mutex);
		Tomography::get_instance()->push_next_raw_data(current_data);
		return has_read;
	}
	else{
		struct ped_data current_data;
		pthread_mutex_lock(&IO_mutex);
		bool has_read = reader->GetNext();
		if(has_read){
			current_data.Nevent = reader->Nevent;
			current_data.evttime = reader->evttime;
			for(map<Tomography::det_type,unsigned short>::iterator type_it=(reader->det_N).begin();type_it!=(reader->det_N).end();++type_it){
				current_data.strip_data[type_it->first] = vector<vector<vector<double> > >(type_it->second,vector<vector<double> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
				for(unsigned short i=0;i<type_it->second;i++){
					if(type == Tomography::ped) current_data.strip_data[type_it->first][i] = reader->get_ampl_ped<double>(type_it->first,i);
					else if(type == Tomography::corr) current_data.strip_data[type_it->first][i] = reader->get_ampl<double>(type_it->first,i);
				}
			}
			pthread_mutex_unlock(&IO_mutex);
		}
		else pthread_mutex_unlock(&IO_mutex);
		if(type == Tomography::ped) Tomography::get_instance()->push_next_ped_data(current_data);
		else if(type == Tomography::corr) Tomography::get_instance()->push_next_corr_data(current_data);
		return has_read;
	}
}
bool Read_Signal_Task::can_exec(){
	if(reader==NULL) return false;
	if(type==Tomography::unknown_signal) return false;
	if(max_event<0) return true;
	return (reader->fChain->GetEntriesFast() < max_event);
}