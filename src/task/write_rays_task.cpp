#define write_rays_task_cpp

#include "task/write_rays_task.h"

#include "Tray.h"
#include "event.h"

Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_): Output_Task<ray_data>(){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
}

Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_, Typed_Task<ray_data> * next_task_): Output_Task<ray_data>(next_task_){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
}

Write_Rays_Task::~Write_Rays_Task(){

}
bool Write_Rays_Task::do_task(){
	ray_data * current_data = get_next_data();
	if(current_data->CBevent == NULL){
		delete current_data;
		return false;
	}
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree((current_data->CBevent)->get_evn(), (current_data->CBevent)->get_evttime(), current_data->rays, z_up, z_down);
	if((data_treated%1000) == 0) writer->Write();
	pthread_mutex_unlock(&IO_mutex);
	if(next_task==NULL) delete current_data;
	else next_task->push_next_data(current_data);
	return true;
}
bool Write_Rays_Task::can_exec() const{
	return (writer!=NULL && (!is_queue_empty()));
}
