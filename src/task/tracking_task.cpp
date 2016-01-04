#define tracking_task_cpp

#include "task/tracking_task.h"

#include "detector.h"
#include "event.h"

Tracking_Abs_Task::Tracking_Abs_Task(const CosmicBench * const detectors_): Task(){
	detectors = detectors_;
}
Tracking_Abs_Task::Tracking_Abs_Task(const CosmicBench * const detectors_, Task * next_task_): Task(next_task_){
	detectors = detectors_;
}
Tracking_Abs_Task::~Tracking_Abs_Task(){

}
bool Tracking_Abs_Task::do_task(){
	struct event_data current_data = Tomography::get_instance()->get_next_event_data();
	if(current_data.Nevent < 0) return false;
	vector<Event*> all_event;
	for(map<Tomography::det_type,vector<Event*> >::iterator type_it = current_data.det_data.begin();type_it!=current_data.det_data.end();++type_it){
		all_event.insert(all_event.end(),(type_it->second).begin(),(type_it->second).end());
	}
	struct ray_data tracked_data;
	tracked_data.CBevent = new CosmicBenchEvent(detectors,all_event);
	for(vector<Event*>::iterator ev_it = all_event.begin();ev_it!=all_event.end();++ev_it){
		delete *ev_it;
	}
	tracked_data.rays = (tracked_data.CBevent)->get_absorption_rays();
	Tomography::get_instance()->push_next_ray_data(tracked_data);
	return true;
}
bool Tracking_Abs_Task::can_exec(){
	return (!(Tomography::get_instance()->is_event_data_empty()));
}
void Tracking_Abs_Task::update_task_list(){
	add_task(next_task);
}

Tracking_Dev_Task::Tracking_Dev_Task(const CosmicBench * const detectors_): Task(){
	detectors = detectors_;
}
Tracking_Dev_Task::Tracking_Dev_Task(const CosmicBench * const detectors_, Task * next_task_): Task(next_task_){
	detectors = detectors_;
}
Tracking_Dev_Task::~Tracking_Dev_Task(){

}
bool Tracking_Dev_Task::do_task(){
	struct event_data current_data = Tomography::get_instance()->get_next_event_data();
	if(current_data.Nevent < 0) return false;
	vector<Event*> all_event;
	for(map<Tomography::det_type,vector<Event*> >::iterator type_it = current_data.det_data.begin();type_it!=current_data.det_data.end();++type_it){
		all_event.insert(all_event.end(),(type_it->second).begin(),(type_it->second).end());
	}
	struct deviation_data tracked_data;
	tracked_data.CBevent = new CosmicBenchEvent(detectors,all_event);
	for(vector<Event*>::iterator ev_it = all_event.begin();ev_it!=all_event.end();++ev_it){
		delete *ev_it;
	}
	(tracked_data.CBevent)->createPairs();
	tracked_data.rays = (tracked_data.CBevent)->get_rayPairs();
	Tomography::get_instance()->push_next_deviation_data(tracked_data);
	return true;
}
bool Tracking_Dev_Task::can_exec(){
	return (!(Tomography::get_instance()->is_event_data_empty()));
}
void Tracking_Dev_Task::update_task_list(){
	add_task(next_task);
}