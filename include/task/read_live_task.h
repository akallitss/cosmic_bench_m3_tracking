#ifndef read_live_task_h
#define read_live_task_h

#include "MT_tomography.h"

#include <string>
using std::string;

class DataReader;
class data_message;

class Read_Live_Task: public Input_Task{
	public:
		Read_Live_Task(string pipe_name);
		~Read_Live_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
		data_message * get_next_data();
		data_message * wait_new_data();
		bool has_new_data() const;
		int get_status() const;
		bool is_saturated() const;
		string init_count() const;
		string print_count() const;
	protected:
		pthread_cond_t queue_cond;
		void * pipe_ptr;
		int queue_id;
		queue<data_message*> data_queue;
		unsigned long data_count;
		int status;
};

#endif
