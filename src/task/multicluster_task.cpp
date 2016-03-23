#define multicluster_task_cpp

#include "task/multicluster_task.h"

#include "detector.h"
#include "event.h"

/*
Multicluster_Task::Multicluster_Task(const CosmicBench * const detectors_): Typed_Task<ped_data>(){
	detectors = detectors_;
	next_task = NULL;
}
*/
Multicluster_Task::Multicluster_Task(const CosmicBench * const detectors_, Typed_Task<event_data> * next_task_): Typed_Task<corr_data>(){
	detectors = detectors_;
	next_task = next_task_;
}
Multicluster_Task::~Multicluster_Task(){
	delete next_task;
}
bool Multicluster_Task::do_task(){
	corr_data * current_data = get_next_data();
	if(current_data->Nevent < 0){
		delete current_data;
		return false;
	}
	event_data * demux_data = new event_data();
	demux_data->Nevent = current_data->Nevent;
	demux_data->evttime = current_data->evttime;
	for(int i=0;i<detectors->get_det_N_tot();i++){
		Detector * det = detectors->get_detector(i);
		demux_data->det_data[det->get_type()].push_back(det->build_event(current_data->strip_data[det->get_type()][det->get_n_in_tree()],demux_data->Nevent,demux_data->evttime));
		(demux_data->det_data[det->get_type()].back())->MultiCluster();
	}
	next_task->push_next_data(demux_data);
	delete current_data;
	return true;
}
void Multicluster_Task::update_task_list() const{
	add_task(next_task);
}
bool Multicluster_Task::is_queueable() const{
	return true;
}
