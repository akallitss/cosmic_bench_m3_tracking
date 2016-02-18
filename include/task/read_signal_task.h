#ifndef read_signal_task_h
#define read_signal_task_h

#include "MT_tomography.h"
#include "tomography.h"

class Tsignal_R;

template<typename T=void>
class Read_Signal_Task: public Input_Task{
	public:
		//Read_Signal_Task(long max_event_, Tsignal_R * reader_);
		Read_Signal_Task(long max_event_, Tsignal_R * reader_, Typed_Task<T> * next_task_);
		~Read_Signal_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
	protected:
		Tsignal_R * reader;
		Typed_Task<T> * next_task;
};

#endif