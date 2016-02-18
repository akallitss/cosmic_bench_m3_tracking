#define read_elec_task_cpp

#include "task/read_elec_task.h"

#include "datareader.h"
/*
Read_Elec_Task::Read_Elec_Task(DataReader * reader_): Input_Task(reader_->get_max_event()){
	reader = reader_;
	next_task = NULL;
}
*/
Read_Elec_Task::Read_Elec_Task(DataReader * reader_, Typed_Task<raw_data> * next_task_): Input_Task(reader_->get_max_event()){
	reader = reader_;
	next_task = next_task_;
}
Read_Elec_Task::~Read_Elec_Task(){

}
bool Read_Elec_Task::do_task(){
	raw_data * current_data = new raw_data();
	pthread_mutex_lock(&IO_mutex);
	bool has_read = !(reader->is_end());
	if(has_read){
		reader->process_event();
		current_data->Nevent = reader->get_event_n();
		current_data->evttime = reader->get_evttime();
		current_data->strip_data = reader->get_data();
		pthread_mutex_unlock(&IO_mutex);
		next_task->push_next_data(current_data);
	}
	else{
		pthread_mutex_unlock(&IO_mutex);
		delete current_data;
	}
	return has_read;
}
bool Read_Elec_Task::can_exec() const{
	if(reader==NULL) return false;
	if(max_event<0) return true;
	return (reader->get_event_n() < max_event);
}
void Read_Elec_Task::update_task_list() const{
	add_task(next_task);
}