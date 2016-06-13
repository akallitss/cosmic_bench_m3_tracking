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

class raw_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<float> > > > strip_data;
		raw_data();
		~raw_data();
};
class ped_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
		ped_data();
		~ped_data();
};
class corr_data{
	public:
		int Nevent;
		double evttime;
		map<Tomography::det_type,vector<vector<vector<double> > > > strip_data;
		corr_data();
		~corr_data();
};
class event_data{
	public:
		map<Tomography::det_type,vector<Event*> > det_data;
		int Nevent;
		double evttime;
		event_data();
		~event_data();
};
class ray_data{
	public:
		CosmicBenchEvent * CBevent;
		vector<Ray> rays;
		ray_data();
		~ray_data();
};
class deviation_data{
	public:
		CosmicBenchEvent * CBevent;
		vector<RayPair> rays;
		deviation_data();
		~deviation_data();
};

class Task{
	public:
		Task();
		//Task(Task * next_task_);
		virtual ~Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() const = 0;
		virtual void update_task_list() const = 0;
		static Task * get_next_task();
		static void add_task(Task * new_task);
		static unsigned int task_left();
		virtual string init_count() const = 0;
		virtual string print_count() const = 0;
		virtual bool is_queueable() const = 0;
	protected:
		static queue<Task*> task_queue;
		static pthread_mutex_t queue_mutex;
		static queue<Task*> queue_init();
		//Task * next_task;
};

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
		void push_next_data(T * next_data);
		unsigned int get_queue_size() const;
		virtual bool is_queueable() const = 0;
	protected:
		T * get_next_data();
		bool is_queue_empty() const;
		queue<T*> data_queue;
		unsigned long data_treated;
		pthread_mutex_t data_queue_mutex;

};
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

class Input_Task: public Task{
	public:
		Input_Task(long max_event_);
		//Input_Task(long max_event_, Task * next_task_);
		~Input_Task();
		virtual bool do_task() = 0;
		virtual bool can_exec() const = 0;
		void update_task_list() const = 0;
		virtual bool is_saturated() const = 0;
		string init_count() const;
		string print_count() const;
		bool is_queueable() const;
	protected:
		pthread_mutex_t IO_mutex;
		//Task * next_task;
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
		static unsigned short has_working_readers();
	protected:
		static unsigned short working_readers;
		static unsigned short working_readers_init();
		static pthread_mutex_t number_mutex;
		void * run();
		void pre_stop();
		bool working;
		Input_Task * current_task;
};

class Writer_Thread: public Thread{
	public:
		template<typename T> Writer_Thread(Output_Task<T> * current_task_);
		~Writer_Thread();
		bool is_working() const;
	protected:
		void * run();
		void pre_stop();
		bool working;
		Task * current_task;
};

class Display_Thread: public Thread, public ostringstream{
	public:
		void set_log_file(string log_file_name);
		bool is_working() const;
		void register_canvas(TCanvas * new_canvas, int canvas_div_n = 0);
		void register_plot(TNamed * new_plot, string canvas_name, string draw_opt = "", int canvas_div = 0);
		//void register_div_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		//void register_sub_hist(TH1 * new_plot_a, TH1 * new_plot_b, string canvas_name, string draw_opt = "", unsigned short canvas_div = 0);
		void start_count();
		void stop_count();
		void register_task(Task * some_task);
		void unregister_task(Task * some_task);
		static Display_Thread * get_instance();
		static void Quit();
		void display_text();
		void display_canvas();
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
