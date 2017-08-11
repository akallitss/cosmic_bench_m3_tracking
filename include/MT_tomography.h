#ifndef tomography_MT_h
#define tomography_MT_h

#include "tomography.h"

#include <pthread.h>

#include <queue>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using std::queue;
using std::ostringstream;
using std::string;
using std::ofstream;
using std::vector;
using std::map;

class TCanvas;
class TObject;

//storage class for raw amplitudes of a full event
class raw_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<float> > > > strip_data;
		raw_data();
		~raw_data();
};
//storage class for ped subtracted amplitudes of a full event
class ped_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
		ped_data();
		~ped_data();
};
//storage class for ped and common noise subtracted amplitudes of a full event
class corr_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
		corr_data();
		~corr_data();
};
//storage class for built event
class event_data{
	public:
		map<Tomography::det_type,vector<Event*> > det_data;
		int Nevent;
		double evttime;
		event_data();
		~event_data();
};
//storage class for straight muon tracks
class ray_data{
	public:
		CosmicBenchEvent * CBevent;
		vector<Ray> rays;
		ray_data();
		~ray_data();
};
//storage class for scattered muon trakcs
class deviation_data{
	public:
		CosmicBenchEvent * CBevent;
		vector<RayPair> rays;
		deviation_data();
		~deviation_data();
};
//base abstract class for what a thread can do as a task
class Task{
	public:
		Task();
		//Task(Task * next_task_);
		virtual ~Task();
		//method to make the thread actually do the job
		virtual bool do_task() = 0;
		//check if the job can be done
		virtual bool can_exec() const = 0;
		//after doing the job, update the task list about further job to do
		virtual void update_task_list() const = 0;
		//retrieve the following job to do
		static Task * get_next_task();
		//add a job to do
		static void add_task(Task * new_task);
		//retrieve the number of job left to be done
		static unsigned int task_left();
		//retrieve count header
		virtual string init_count() const = 0;
		//retrieve actual task count (done and left to be done)
		virtual string print_count() const = 0;
		//check if it is possible to add this task to the task list
		virtual bool is_queueable() const = 0;
	protected:
		static queue<Task*> task_queue;
		static pthread_mutex_t queue_mutex;
		static queue<Task*> queue_init();
		//Task * next_task;
};

//base class for tasks that need T type of information to be done
template<typename T>
class Typed_Task: public Task{
	public:
		Typed_Task();
		//Typed_Task(Task * next_task_);
		virtual ~Typed_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() const;
		virtual void update_task_list() const = 0;
		virtual string init_count() const;
		string print_count() const;
		//store data to be processed
		void push_next_data(T * next_data);
		//retrieve data queue size
		unsigned int get_queue_size() const;
		virtual bool is_queueable() const = 0;
	protected:
		//retrieve the data to process it
		T * get_next_data();
		//check if there is data to be processed
		bool is_queue_empty() const;
		queue<T*> data_queue;
		unsigned long data_treated;
		pthread_mutex_t data_queue_mutex;

};
//task class that store processed data to be used later
template<typename T>
class Buffer_Task: public Typed_Task<T>{
	public:
		Buffer_Task();
		~Buffer_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
		T * fetch_data();
		bool can_fetch_data() const;
		bool is_queueable() const;
};

//base class for tasks that output T type of information in files or elsewhere
template<typename T>
class Output_Task: public Typed_Task<T>{
	public:
		Output_Task();
		Output_Task(Typed_Task<T> * next_task_);
		~Output_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() const = 0;
		void update_task_list() const;
		string init_count() const;
		bool is_queueable() const;
	protected:
		pthread_mutex_t IO_mutex;
		Typed_Task<T> * next_task;
};

//base class for tasks that fetch information in files, message queues or elsewhere
class Input_Task: public Task{
	public:
		Input_Task(long max_event_);
		//Input_Task(long max_event_, Task * next_task_);
		~Input_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() const = 0;
		void update_task_list() const = 0;
		virtual bool is_saturated() const = 0;
		virtual string init_count() const;
		virtual string print_count() const;
		bool is_queueable() const;
	protected:
		pthread_mutex_t IO_mutex;
		//Task * next_task;
		long max_event;
};

//base abstract class for thread that will execute tasks
class Thread{
	public:
		Thread();
		virtual ~Thread();
		//spawn the thread
		int start();
		//force to stop (kill) the thread
		int stop();
		//retrieve thread id
		pthread_t getThreadId() const;
		//check if thread is running and doing tasks
		virtual bool is_working() const = 0;
		//retrieve thread name
		virtual string get_name() const = 0;
	private:
		//method executed by the thread
		virtual void * run() = 0;
		//pthread compliant function to spawn the thread (call run())
		static void * runThread(void * arg);
		//called just before the thread kill in order to clean up
		virtual void pre_stop();
		pthread_attr_t attr;
		pthread_t id;
		int running;
		int detached;
};

//class for thread that will execute non-I/O jobs
class Worker_Thread: public Thread{
	public:
		Worker_Thread();
		~Worker_Thread();
		bool is_working() const;
		string get_name() const;
	protected:
		void * run();
		void pre_stop();
		bool working;
		Task * current_task;
};

//class for thread that will execute Input_Task
class Reader_Thread: public Thread{
	public:
		//constructor with the input task to do
		Reader_Thread(Input_Task * current_task_);
		~Reader_Thread();
		bool is_working() const;
		//check if there is still working readers that will output further data to be processed
		static unsigned short has_working_readers();
		string get_name() const;
	protected:
		static unsigned short working_readers;
		static unsigned short working_readers_init();
		static pthread_mutex_t number_mutex;
		void * run();
		void pre_stop();
		bool working;
		Input_Task * current_task;
};

//class for thread that will execute Output_Task
class Writer_Thread: public Thread{
	public:
		//constructor with the output task to do
		template<typename T> Writer_Thread(Output_Task<T> * current_task_);
		~Writer_Thread();
		bool is_working() const;
		string get_name() const;
	protected:
		void * run();
		void pre_stop();
		bool working;
		Task * current_task;
};

//specific singleton class to display thread informations during the process
class Display_Thread: public Thread, public ostringstream{
	public:
		//set log file to save all the text output
		void set_log_file(string log_file_name);
		bool is_working() const;
		//register canvas to display (currently broken)
		void register_canvas(TCanvas * new_canvas, int canvas_div_n = 0);
		//register plot to display (TGraph,TH,TProfile,...) in the given canvas
		void register_plot(TNamed * new_plot, string canvas_name, string draw_opt = "", int canvas_div = 0);
		//void register_div_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		//void register_sub_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		//start to regularly output task information
		void start_count();
		//stop to output task information
		void stop_count();
		//register a task to be able to output its info
		void register_task(Task * some_task);
		//unregister a task to stop to output its info
		void unregister_task(Task * some_task);
		//get the singleton instance to feed additionnal text to output
		static Display_Thread * get_instance();
		//quit and stop the display thread
		static void Quit();
		//actually display the information text
		void display_text();
		//display and update the registered canvas with their plots
		void display_canvas();
		string get_name() const;
	protected:
		static Display_Thread * singleton_instance;
		Display_Thread();
		Display_Thread(string log_file_name);
		~Display_Thread();
		struct plot_info{
			TNamed * plot_orig;
			TNamed * plot_clone;
			int div;
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
			int div_n;
			vector<plot_info> plots;
			//vector<divided_hist> div_plots;
			//vector<substracted_hist> sub_plots;
		};
		void * run();
		void pre_stop();
		bool working;
		bool is_counting;
		ofstream log_file;
		vector<canvas_info> canvas_list;
		vector<Task*> registered_task;

};

#endif
