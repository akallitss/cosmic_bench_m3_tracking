#ifndef tracking_task_h
#define tracking_task_h

#include "MT_tomography.h"

class CosmicBench;

class Tracking_Abs_Task: public Typed_Task<event_data>{
	public:
		//Tracking_Abs_Task(const CosmicBench * const detectors_);
		Tracking_Abs_Task(const CosmicBench * const detectors_, Typed_Task<ray_data> * next_task_);
		~Tracking_Abs_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		const CosmicBench * detectors;
		Typed_Task<ray_data> * next_task;
};
class Tracking_Dev_Task: public Typed_Task<event_data>{
	public:
		//Tracking_Dev_Task(const CosmicBench * const detectors_);
		Tracking_Dev_Task(const CosmicBench * const detectors_, Typed_Task<deviation_data> * next_task_);
		~Tracking_Dev_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		const CosmicBench * detectors;
		Typed_Task<deviation_data> * next_task;
};

#endif