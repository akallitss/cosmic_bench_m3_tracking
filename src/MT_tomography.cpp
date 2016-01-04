#define tomography_MT_cpp

#include "MT_tomography.h"

#include "datareader.h"
#include "detector.h"
#include "event.h"
#include "Tsignal_R.h"
#include "Tanalyse_R.h"
#include "Tray.h"
#include "Tsignal_W.h"
#include "Tanalyse_W.h"

#include <TObject.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;
using std::flush;

Task::Task(){
	next_task = NULL;
}
Task::Task(Task * next_task_){
	next_task = next_task_;
}
Task::~Task(){

}
Task * Task::get_next_task(){
	Task * priority_task;
	pthread_mutex_lock(&queue_mutex);
	if(task_queue.empty()) priority_task = NULL;
	else{
		priority_task = task_queue.front();
		task_queue.pop();
	}
	pthread_mutex_unlock(&queue_mutex);
	return priority_task;
}
bool Task::is_init = false;
pthread_mutex_t Task::queue_mutex;
queue<Task*> Task::task_queue;

void Task::add_task(Task * new_task){
	if(!is_init){
		pthread_mutex_init(&queue_mutex, NULL);
		is_init = true;
	}
	pthread_mutex_lock(&queue_mutex);
	task_queue.push(new_task);
	pthread_mutex_unlock(&queue_mutex);
}
unsigned int Task::task_left(){
	return task_queue.size();
}
IO_Task::IO_Task(): Task(){
	pthread_mutex_init(&IO_mutex, NULL);
}
IO_Task::IO_Task(Task * next_task_): Task(next_task_){
	pthread_mutex_init(&IO_mutex, NULL);
}
IO_Task::~IO_Task(){
	pthread_mutex_destroy(&IO_mutex);
}

Input_Task::Input_Task(long max_event_){
	pthread_mutex_init(&IO_mutex, NULL);
	next_task = NULL;
	max_event = max_event_;
}
Input_Task::Input_Task(long max_event_, Task * next_task_){
	pthread_mutex_init(&IO_mutex, NULL);
	next_task = next_task_;
	max_event = max_event_;
}
Input_Task::~Input_Task(){
	pthread_mutex_destroy(&IO_mutex);
}
void Input_Task::update_task_list(){
	Task::add_task(next_task);
}

Thread::Thread(){
	id = 0;
	running = 0;
	detached = 0;
}
Thread::~Thread(){
	if(running ==1 && detached == 0) pthread_detach(id);
	if(running == 1) pthread_cancel(id);
}
void * Thread::runThread(void * arg){
	return ((Thread*) arg)->run();
}
int Thread::start(){
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	if(running == 0){
		int result = pthread_create(&id,NULL,runThread,this);
		if(result == 0) running = 1;
		pthread_attr_destroy(&attr);
		return result;
	}
	else return -1;
}
void Thread::pre_stop(){

}
int Thread::stop(){
	if(running == 0){
		return -1;
	}
	pre_stop();
	void * status;
	int result = pthread_join(id,&status);
	if(result == 0) running = 0;
	return result;
}
pthread_t Thread::getThreadId() const{
	return id;
}

Worker_Thread::Worker_Thread(): Thread(){
	working = false;
	current_task = NULL;
}
Worker_Thread::~Worker_Thread(){
	working = false;
	current_task = NULL;
}
void * Worker_Thread::run(){
	working = true;
	unsigned int wait_time = 0;
	while(working){
		current_task = Task::get_next_task();
		if(current_task!=NULL){
			wait_time = 0;
			while(!(current_task->can_exec())) usleep(100);
			if(current_task->do_task()) current_task->update_task_list();
		}
		else{
			usleep(100);
			wait_time++;
			if(wait_time>50000) working = false;
		}
	}
	return 0;
}
void Worker_Thread::pre_stop(){
	working = false;
}
bool Worker_Thread::is_working() const{
	return working;
}

Reader_Thread::Reader_Thread(Input_Task * current_task_): Thread(){
	working = false;
	current_task = current_task_;
}
Reader_Thread::~Reader_Thread(){
	working = false;
	current_task = NULL;
}
bool Reader_Thread::is_working() const{
	return working;
}
void * Reader_Thread::run(){
	working = current_task!=NULL;
	unsigned int wait_time = 0;
	while(working){
		if(current_task->can_exec()){
			wait_time = 0;
			while(!(current_task->do_task())){
				usleep(1000);
				wait_time++;
				if(wait_time>10){
					working = false;
					break;
				}
			}
			if(working) current_task->update_task_list();
		}
	}
	return 0;
}
void Reader_Thread::pre_stop(){
	working = false;
}

Display_Thread::Display_Thread(): Thread(){
	start();
}
Display_Thread::Display_Thread(string log_file_name): Thread(), ostringstream(){
	log_file.open(log_file_name.c_str());
	start();
}
Display_Thread::~Display_Thread(){
	stop();
	if(log_file.is_open() && log_file.good()){
		log_file.flush();
		log_file.close();
	}
	cout << endl;
}
void Display_Thread::set_log_file(string log_file_name){
	if(log_file.is_open() && log_file.good()){
		log_file.flush();
		log_file.close();
	}
	log_file.open(log_file_name.c_str());
}
bool Display_Thread::is_working() const{
	return working;
}
void Display_Thread::register_canvas(TCanvas * new_canvas, unsigned short canvas_div_n){
	for(vector<canvas_info>::iterator it=canvas_list.begin();it!=canvas_list.end();++it){
		if(new_canvas->GetName() == it->addr->GetName()){
			cout << "canvas \"" << new_canvas->GetName() << "\" already registered !" << endl;
			return;
		}
	}
	struct canvas_info next_canvas;
	next_canvas.addr = new_canvas;
	next_canvas.div_n = canvas_div_n;
	canvas_list.push_back(next_canvas);
}
void Display_Thread::register_plot(TObject * new_plot, string canvas_name, string draw_opt, unsigned short canvas_div){
	for(vector<canvas_info>::iterator it=canvas_list.begin();it!=canvas_list.end();++it){
		if((canvas_name == it->addr->GetName()) && (canvas_div <= it->div_n)){
			struct plot_info next_plot;
			next_plot.plot = new_plot;
			next_plot.div = canvas_div;
			next_plot.draw_opt = draw_opt;
			(it->plots).push_back(next_plot);
			return;
		}
	}
}
void Display_Thread::pre_stop(){
	working = false;
	display_text();
	display_canvas();
}
void * Display_Thread::run(){
	unsigned int delay = 0;
	while(working){
		display_text();
		if(delay>300){
			display_canvas();
		}
		usleep(10000);
		delay++;
	}
	return 0;
}
void Display_Thread::display_text(){
	string buffer = str();
	str("");
	cout << buffer;
	cout.flush();
	if(log_file.is_open() && log_file.good()){
		log_file << buffer;
		log_file.flush();
	}
}
void Display_Thread::display_canvas(){
	for(vector<canvas_info>::iterator it=canvas_list.begin();it!=canvas_list.end();++it){
		for(vector<plot_info>::iterator jt = (it->plots).begin();jt!=(it->plots).end();++jt){
			it->addr->cd(jt->div);
			jt->plot->DrawClone((jt->draw_opt).c_str());
		}
		it->addr->Modified();
		it->addr->Update();
	}
}