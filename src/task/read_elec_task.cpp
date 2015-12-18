#define read_elec_task_cpp

#include "task/read_elec_task.h"

#include "datareader.h"

Read_Elec_Task::Read_Elec_Task(DataReader * reader_): Input_Task(reader_->get_max_event()){
	reader = reader_;
}
Read_Elec_Task::Read_Elec_Task(DataReader * reader_, Task * next_task_): Input_Task(reader_->get_max_event(), next_task_){
	reader = reader_;
}
Read_Elec_Task::~Read_Elec_Task(){

}
bool Read_Elec_Task::do_task(){
	struct raw_data current_data;
	pthread_mutex_lock(&IO_mutex);
	bool has_read = !(reader->is_end());
	if(has_read){
		reader->process_event();
		current_data.Nevent = reader->get_event_n();
		current_data.evttime = reader->get_evttime();
		current_data.strip_data = reader->get_data();
		pthread_mutex_unlock(&IO_mutex);
		Tomography::get_instance()->push_next_raw_data(current_data);
	}
	else pthread_mutex_unlock(&IO_mutex);
	return has_read;
}
bool Read_Elec_Task::can_exec(){
	if(reader==NULL) return false;
	if(max_event<0) return true;
	return (reader->get_event_n() < max_event);
}