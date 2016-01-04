#ifndef tomography_MT_h
#define tomography_MT_h

#include "tomography.h"

#include <pthread.h>

#include <queue>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using std::queue;
using std::ostringstream;
using std::string;
using std::ofstream;
using std::vector;

class TCanvas;
class TObject;

class Task{
	public:
		Task();
		Task(Task * next_task_);
		virtual ~Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() = 0;
		virtual void update_task_list() = 0;
		static Task * get_next_task();
		static void add_task(Task * new_task);
		static unsigned int task_left();
	protected:
		static queue<Task*> task_queue;
		static pthread_mutex_t queue_mutex;
		static bool is_init;
		Task * next_task;
};

class IO_Task: public Task{
	public:
		IO_Task();
		IO_Task(Task * next_task_);
		~IO_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() = 0;
		virtual void update_task_list() = 0;
	protected:
		pthread_mutex_t IO_mutex;

};

class Input_Task{
	public:
		Input_Task(long max_event_);
		Input_Task(long max_event_, Task * next_task_);
		~Input_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() = 0;
		void update_task_list();
	protected:
		pthread_mutex_t IO_mutex;
		Task * next_task;
		long max_event;
};

class Thread{
	public:
		Thread();
		virtual ~Thread();
		int start();
		int stop();
		pthread_t getThreadId() const;
		virtual bool is_working() const = 0;
	private:
		virtual void * run() = 0;
		static void * runThread(void * arg);
		virtual void pre_stop();
		pthread_attr_t attr;
		pthread_t id;
		int running;
		int detached;
};

class Worker_Thread: public Thread{
	public:
		Worker_Thread();
		~Worker_Thread();
		bool is_working() const;
	protected:
		void * run();
		void pre_stop();
		bool working;
		Task * current_task;
};

class Reader_Thread: public Thread{
	public:
		Reader_Thread(Input_Task * current_task_);
		~Reader_Thread();
		bool is_working() const;
	protected:
		void * run();
		void pre_stop();
		bool working;
		Input_Task * current_task;
};

class Display_Thread: public Thread, public ostringstream{
	public:
		Display_Thread();
		Display_Thread(string log_file_name);
		~Display_Thread();
		void set_log_file(string log_file_name);
		bool is_working() const;
		void register_canvas(TCanvas * new_canvas, unsigned short canvas_div_n = 0);
		void register_plot(TObject * new_plot, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		//void register_div_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		//void register_sub_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
	protected:
		struct plot_info{
			TObject * plot;
			unsigned short div;
			string draw_opt;
		};
		/*
		struct divided_hist{
			unsigned short div;
			string draw_opt;
			TH1 * hist_a;
			TH1 * hist_b;
		};
		struct substracted_hist{
			unsigned short div;
			string draw_opt;
			TH1 * hist_a;
			TH1 * hist_b;
		};
		*/
		struct canvas_info{
			TCanvas * addr;
			unsigned short div_n;
			vector<plot_info> plots;
			//vector<divided_hist> div_plots;
			//vector<substracted_hist> sub_plots;
		};
		void * run();
		void pre_stop();
		bool working;
		ofstream log_file;
		vector<canvas_info> canvas_list;
		void display_text();
		void display_canvas();

};

#endif