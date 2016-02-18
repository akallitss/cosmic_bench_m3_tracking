#ifndef read_analyse_task_h
#define read_analyse_task_h

#include "MT_tomography.h"

class Tanalyse_R;
class CosmicBench;

class Read_Analyse_Task: public Input_Task{
	public:
		//Read_Analyse_Task(long max_event_, Tanalyse_R * reader_, const CosmicBench * const detectors_);
		Read_Analyse_Task(long max_event_, Tanalyse_R * reader_, const CosmicBench * const detectors_, Typed_Task<event_data> * next_task_);
		~Read_Analyse_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
	protected:
		const CosmicBench * detectors;
		Tanalyse_R * reader;
		Typed_Task<event_data> * next_task;
};

#endif