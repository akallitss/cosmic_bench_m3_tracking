#define write_rays_task_cpp

#include "task/write_rays_task.h"

#include "Tray.h"

Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_): Output_Task<ray_data>(){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
}
/*
Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_, Typed_Task<ray_data> * next_task_): Output_Task<ray_data>(){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
}
*/
Write_Rays_Task::~Write_Rays_Task(){

}
bool Write_Rays_Task::do_task(){
	ray_data * current_data = get_next_data();
	if(current_data->CBevent == NULL) return false;
	pthread_mutex_lock(&IO_mutex);
	writer->fillTree((current_data->CBevent)->get_evn(), (current_data->CBevent)->get_evttime(), current_data->rays, z_up, z_down);
	pthread_mutex_unlock(&IO_mutex);
	delete current_data;
	return true;
}
bool Write_Rays_Task::can_exec() const{
	return (writer!=NULL && (!is_queue_empty()));
}
void Write_Rays_Task::update_task_list() const{
	//add_task(next_task);
}