#define read_live_task_cpp

#include <sys/msg.h>

#include "task/read_live_task.h"

#include "ElecReader.h"

Read_Live_Task::Read_Live_Task(int queue_id_): Input_Task(-1){
	queue_id = queue_id_;
	status = 0;
}
Read_Live_Task::~Read_Live_Task(){

}
bool Read_Live_Task::do_task(){
	data_message * current_data = new data_message();
	if(msgrcv(queue_id,current_data,sizeof(current_data->data),0,0) == -1) return false;
	if(current_data->mtype == 2){
		pthread_mutex_lock(&IO_mutex);
		data_queue.push(current_data);
		pthread_mutex_unlock(&IO_mutex);
		return true;
	}
	else if(current_data->mtype == 1){
		status |= (current_data->data)[0];
		delete current_data;
		return true;
	}
	delete current_data;
	return false;
}
bool Read_Live_Task::can_exec() const{
	return (!(status & 0x8000));
}
void Read_Live_Task::update_task_list() const{

}
data_message * Read_Live_Task::get_next_data(){
	data_message * current_data;
	pthread_mutex_lock(&IO_mutex);
	current_data = data_queue.front();
	data_queue.pop();
	pthread_mutex_unlock(&IO_mutex);
	return current_data;
}
bool Read_Live_Task::has_new_data() const{
	return !(data_queue.empty());
}
int Read_Live_Task::get_status() const{
	return status;
}