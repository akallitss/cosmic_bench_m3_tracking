#define tracking_task_cpp

#include "task/tracking_task.h"

#include "detector.h"
#include "event.h"
/*
Tracking_Abs_Task::Tracking_Abs_Task(const CosmicBench * const detectors_): Typed_Task<event_data>(){
	detectors = detectors_;
	next_task = NULL;
}
*/
Tracking_Abs_Task::Tracking_Abs_Task(const CosmicBench * const detectors_, Typed_Task<ray_data> * next_task_): Typed_Task<event_data>(){
	detectors = detectors_;
	next_task = next_task_;
}
Tracking_Abs_Task::~Tracking_Abs_Task(){

}
bool Tracking_Abs_Task::do_task(){
	event_data * current_data = get_next_data();
	if(current_data->Nevent < 0) return false;
	vector<Event*> all_event;
	for(map<Tomography::det_type,vector<Event*> >::iterator type_it = (current_data->det_data).begin();type_it!=(current_data->det_data).end();++type_it){
		all_event.insert(all_event.end(),(type_it->second).begin(),(type_it->second).end());
	}
	ray_data * tracked_data = new ray_data();
	tracked_data->CBevent = new CosmicBenchEvent(detectors,all_event);
	tracked_data->rays = (tracked_data->CBevent)->get_absorption_rays();
	next_task->push_next_data(tracked_data);
	delete current_data;
	return true;
}
bool Tracking_Abs_Task::can_exec() const{
	return (!is_queue_empty());
}
void Tracking_Abs_Task::update_task_list() const{
	add_task(next_task);
}
/*
Tracking_Dev_Task::Tracking_Dev_Task(const CosmicBench * const detectors_): Typed_Task<event_data>(){
	detectors = detectors_;
	next_task = NULL;
}
*/
Tracking_Dev_Task::Tracking_Dev_Task(const CosmicBench * const detectors_, Typed_Task<deviation_data> * next_task_): Typed_Task<event_data>(){
	detectors = detectors_;
	next_task = next_task_;
}
Tracking_Dev_Task::~Tracking_Dev_Task(){

}
bool Tracking_Dev_Task::do_task(){
	event_data * current_data = get_next_data();
	if(current_data->Nevent < 0) return false;
	vector<Event*> all_event;
	for(map<Tomography::det_type,vector<Event*> >::iterator type_it = (current_data->det_data).begin();type_it!=(current_data->det_data).end();++type_it){
		all_event.insert(all_event.end(),(type_it->second).begin(),(type_it->second).end());
	}
	deviation_data * tracked_data = new deviation_data();
	tracked_data->CBevent = new CosmicBenchEvent(detectors,all_event);
	(tracked_data->CBevent)->createPairs();
	tracked_data->rays = (tracked_data->CBevent)->get_rayPairs();
	next_task->push_next_data(tracked_data);
	delete current_data;
	return true;
}
bool Tracking_Dev_Task::can_exec() const{
	return (!is_queue_empty());
}
void Tracking_Dev_Task::update_task_list() const{
	add_task(next_task);
}