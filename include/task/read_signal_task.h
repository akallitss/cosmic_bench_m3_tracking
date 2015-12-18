#ifndef read_signal_task_h
#define read_signal_task_h

#include "MT_tomography.h"
#include "tomography.h"

class Tsignal_R;

class Read_Signal_Task: public Input_Task{
	public:
		Read_Signal_Task(long max_event_, Tsignal_R * reader_, Tomography::signal_type type_);
		Read_Signal_Task(long max_event_, Tsignal_R * reader_, Tomography::signal_type type_,Task * next_task_);
		~Read_Signal_Task();
		bool do_task();
		bool can_exec();
	protected:
		Tomography::signal_type type;
		Tsignal_R * reader;
};

#endif