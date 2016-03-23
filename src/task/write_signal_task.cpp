#define write_signal_task_cpp

#include "task/write_signal_task.h"

#include "Tsignal_W.h"

template<typename T>
Write_Signal_Task<T>::Write_Signal_Task(Tsignal_W * writer_): Output_Task<T>(){
	writer = writer_;
}

template<typename T>
Write_Signal_Task<T>::Write_Signal_Task(Tsignal_W * writer_, Typed_Task<T> * next_task_): Output_Task<T>(next_task_){
	writer = writer_;
}

template<typename T>
Write_Signal_Task<T>::~Write_Signal_Task(){

}
template<>
bool Write_Signal_Task<raw_data>::do_task(){
	raw_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree_raw(current_data->Nevent, current_data->evttime, current_data->strip_data);
	pthread_mutex_unlock(&IO_mutex);
	if(next_task==NULL) delete current_data;
	else next_task->push_next_data(current_data);
	return true;
}
template<>
bool Write_Signal_Task<ped_data>::do_task(){
	ped_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree_ped(current_data->strip_data);
	pthread_mutex_unlock(&IO_mutex);
	if(next_task==NULL) delete current_data;
	else next_task->push_next_data(current_data);
	return true;
}
template<>
bool Write_Signal_Task<corr_data>::do_task(){
	corr_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree_corr(current_data->strip_data);
	pthread_mutex_unlock(&IO_mutex);
	if(next_task==NULL) delete current_data;
	else next_task->push_next_data(current_data);
	return true;
}
template<typename T>
bool Write_Signal_Task<T>::can_exec() const{
	return ((writer!=NULL) && !(this->is_queue_empty()));
}
template class Write_Signal_Task<raw_data>;
template class Write_Signal_Task<ped_data>;
template class Write_Signal_Task<corr_data>;
