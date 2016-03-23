#ifndef write_signal_task_h
#define write_signal_task_h

#include "MT_tomography.h"
#include "tomography.h"

class Tsignal_W;

template<typename T=void>
class Write_Signal_Task: public Output_Task<T>{
	public:
		Write_Signal_Task(Tsignal_W * writer_);
		Write_Signal_Task(Tsignal_W * writer_, Typed_Task<T> * next_task_);
		~Write_Signal_Task();
		bool do_task();
		bool can_exec() const;
	protected:
		Tsignal_W * writer;
};

#endif