#ifndef write_analyse_task_h
#define write_analyse_task_h

#include "MT_tomography.h"

class Tanalyse_W;

class Write_Analyse_Task: public Output_Task<event_data>{
	public:
		Write_Analyse_Task(Tanalyse_W * writer_);
		Write_Analyse_Task(Tanalyse_W * writer_, Typed_Task<event_data> * next_task_);
		~Write_Analyse_Task();
		bool do_task();
		bool can_exec() const;
	protected:
		Tanalyse_W * writer;
};

#endif