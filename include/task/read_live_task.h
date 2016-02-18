#ifndef read_live_task_h
#define read_live_task_h

#include "MT_tomography.h"

class DataReader;
class data_message;

class Read_Live_Task: public Input_Task{
	public:
		Read_Live_Task(int queue_id_);
		~Read_Live_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
		data_message * get_next_data();
		bool has_new_data() const;
		int get_status() const;
	protected:
		int queue_id;
		queue<data_message*> data_queue;
		int status;
};

#endif