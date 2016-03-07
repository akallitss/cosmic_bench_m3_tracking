#define read_signal_task_cpp

#include "task/read_signal_task.h"

#include "Tsignal_R.h"
/*
Read_Signal_Task::Read_Signal_Task(long max_event_, Tsignal_R * reader_): Input_Task(max_event_){
	reader = reader_;
	next_task = NULL;
}
*/
template<typename T>
Read_Signal_Task<T>::Read_Signal_Task(long max_event_, Tsignal_R * reader_,Typed_Task<T> * next_task_): Input_Task(max_event_){
	reader = reader_;
	next_task = next_task_;
}

template<typename T>
Read_Signal_Task<T>::~Read_Signal_Task(){

}

template<>
bool Read_Signal_Task<raw_data>::do_task(){
	raw_data * current_data = new raw_data();
	pthread_mutex_lock(&IO_mutex);
	bool has_read = reader->GetNext();
	if(has_read){
		current_data->Nevent = reader->Nevent;
		current_data->evttime = reader->evttime;
		for(map<Tomography::det_type,unsigned short>::iterator type_it=(reader->det_N).begin();type_it!=(reader->det_N).end();++type_it){
			(current_data->strip_data)[type_it->first] = vector<vector<vector<float> > >(type_it->second,vector<vector<float> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<float>(Tomography::get_instance()->get_Nsample(),0)));
			for(unsigned short i=0;i<type_it->second;i++){
				(current_data->strip_data)[type_it->first][i] = reader->get_ampl_raw<float>(type_it->first,i);
			}
		}
		pthread_mutex_unlock(&IO_mutex);
		next_task->push_next_data(current_data);
	}
	else{
		pthread_mutex_unlock(&IO_mutex);
		delete current_data;
	}
	return has_read;
}
template<>
bool Read_Signal_Task<ped_data>::do_task(){
	ped_data * current_data = new ped_data();
	pthread_mutex_lock(&IO_mutex);
	bool has_read = reader->GetNext();
	if(has_read){
		current_data->Nevent = reader->Nevent;
		current_data->evttime = reader->evttime;
		for(map<Tomography::det_type,unsigned short>::iterator type_it=(reader->det_N).begin();type_it!=(reader->det_N).end();++type_it){
			(current_data->strip_data)[type_it->first] = vector<vector<vector<double> > >(type_it->second,vector<vector<double> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
			for(unsigned short i=0;i<type_it->second;i++){
				(current_data->strip_data)[type_it->first][i] = reader->get_ampl_ped<double>(type_it->first,i);
			}
		}
		pthread_mutex_unlock(&IO_mutex);
		next_task->push_next_data(current_data);
	}
	else{
		pthread_mutex_unlock(&IO_mutex);
		delete current_data;
	}
	return has_read;
}
template<>
bool Read_Signal_Task<corr_data>::do_task(){
	corr_data * current_data = new corr_data();
	pthread_mutex_lock(&IO_mutex);
	bool has_read = reader->GetNext();
	if(has_read){
		current_data->Nevent = reader->Nevent;
		current_data->evttime = reader->evttime;
		for(map<Tomography::det_type,unsigned short>::iterator type_it=(reader->det_N).begin();type_it!=(reader->det_N).end();++type_it){
			(current_data->strip_data)[type_it->first] = vector<vector<vector<double> > >(type_it->second,vector<vector<double> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
			for(unsigned short i=0;i<type_it->second;i++){
				(current_data->strip_data)[type_it->first][i] = reader->get_ampl<double>(type_it->first,i);
			}
		}
		pthread_mutex_unlock(&IO_mutex);
		next_task->push_next_data(current_data);
	}
	else{
		pthread_mutex_unlock(&IO_mutex);
		delete current_data;
	}
	return has_read;
}
template<typename T>
bool Read_Signal_Task<T>::can_exec() const{
	if(reader==NULL) return false;
	if(max_event<0) return true;
	return (reader->fChain->GetEntriesFast() < max_event);
}
template<typename T>
void Read_Signal_Task<T>::update_task_list() const{
	add_task(next_task);
}

template class Read_Signal_Task<raw_data>;
template class Read_Signal_Task<ped_data>;
template class Read_Signal_Task<corr_data>;