#define write_signal_task_cpp

#include "task/write_signal_task.h"

#include "Tsignal_W.h"

template<typename T>
Write_Signal_Task<T>::Write_Signal_Task(Tsignal_W * writer_): Output_Task<T>(){
	writer = writer_;
	type = type_;
}

/*
template<typename T>
Write_Signal_Task<T>::Write_Signal_Task(Tsignal_W * writer_, Task * next_task_): Output_Task<T>(){
	writer = writer_;
	type = type_;
}
*/
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
	delete current_data;
	return true;
}
template<>
bool Write_Signal_Task<ped_data>::do_task(){
	ped_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree_ped(current_data->strip_data);
	pthread_mutex_unlock(&IO_mutex);
	delete current_data;
	return true;
}
template<>
bool Write_Signal_Task<corr_data>::do_task(){
	corr_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree_corr(current_data->strip_data);
	pthread_mutex_unlock(&IO_mutex);
	delete current_data;
	return true;
}
template<typename T>
bool Write_Signal_Task<T>::can_exec() const{
	return ((writer!=NULL) && !is_queue_empty());
}
template<typename T>
void Write_Signal_Task<T>::update_task_list() const{
	//add_task(next_task);
}