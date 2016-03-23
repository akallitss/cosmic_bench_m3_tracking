#ifndef multicluster_task_h
#define multicluster_task_h

#include "MT_tomography.h"

class CosmicBench;

class Multicluster_Task: public Typed_Task<corr_data>{
	public:
		//Multicluster_Task(const CosmicBench * const detectors_);
		Multicluster_Task(const CosmicBench * const detectors_, Typed_Task<event_data> * next_task_);
		~Multicluster_Task();
		bool do_task();
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		const CosmicBench * detectors;
		Typed_Task<event_data> * next_task;
};

#endif
