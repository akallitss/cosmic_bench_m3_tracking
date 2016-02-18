#define write_analyse_task_cpp

#include "task/write_analyse_task.h"

#include "Tanalyse_W.h"
#include "event.h"

Write_Analyse_Task::Write_Analyse_Task(Tanalyse_W * writer_): Output_Task<event_data>(){
	writer = writer_;
}
/*
Write_Analyse_Task::Write_Analyse_Task(Tanalyse_W * writer_, Typed_Task<event_data>() * next_task_): Output_Task<event_data>(){
	writer = writer_;
}
*/
Write_Analyse_Task::~Write_Analyse_Task(){

}
bool Write_Analyse_Task::do_task(){
	event_data * current_data = get_next_data();
	if(current_data->Nevent<0) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree(current_data->Nevent, current_data->evttime, current_data->det_data);
	pthread_mutex_unlock(&IO_mutex);
	delete current_data;
	return true;
}
bool Write_Analyse_Task::can_exec() const{
	return (writer!=NULL && (!is_queue_empty()));
}
void Write_Analyse_Task::update_task_list() const{
	//add_task(next_task);
}