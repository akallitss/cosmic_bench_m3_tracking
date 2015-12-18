#define read_analyse_task_cpp

#include "task/read_analyse_task.h"

#include "Tanalyse_R.h"
#include "detector.h"

Read_Analyse_Task::Read_Analyse_Task(long max_event_, Tanalyse_R * reader_, const CosmicBench * const detectors_): Input_Task(max_event_){
	detectors = detectors_;
	reader = reader_;
}
Read_Analyse_Task::Read_Analyse_Task(long max_event_, Tanalyse_R * reader_, const CosmicBench * const detectors_, Task * next_task_): Input_Task(max_event_, next_task_){
	detectors = detectors_;
	reader = reader_;
}
Read_Analyse_Task::~Read_Analyse_Task(){

}
bool Read_Analyse_Task::do_task(){
	int det_N = detectors->get_det_N_tot();
	struct event_data current_data;
	pthread_mutex_lock(&IO_mutex);
	bool has_read = reader->GetNext();
	if(has_read){
		current_data.Nevent = reader->evn;
		current_data.evttime = reader->evttime;
		for(int i=0;i<det_N;i++){
			//events.push_back(detectors->get_detector(i)->build_event(treeObject));
			Detector * current_det = detectors->get_detector(i);
			current_data.det_data[current_det->get_type()].push_back(current_det->build_event(reader));
		}
		pthread_mutex_unlock(&IO_mutex);
		Tomography::get_instance()->push_next_event_data(current_data);
	}
	else pthread_mutex_unlock(&IO_mutex);
	return has_read;
}
bool Read_Analyse_Task::can_exec(){
	if(reader==NULL) return false;
	if(detectors==NULL) return false;
	if(max_event<0) return true;
	return (reader->fChain->GetEntriesFast() < max_event);
}