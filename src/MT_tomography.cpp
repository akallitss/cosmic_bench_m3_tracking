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
#include <iomanip>

using std::vector;
using std::cout;
using std::endl;
using std::flush;
using std::right;
using std::left;
using std::setw;

raw_data::raw_data(){
	Nevent = -1;
	evttime = 0;
}
raw_data::~raw_data(){

}
ped_data::ped_data(){
	Nevent = -1;
	evttime = 0;
}
ped_data::~ped_data(){

}
corr_data::corr_data(){
	Nevent = -1;
	evttime = 0;
}
corr_data::~corr_data(){

}
event_data::event_data(){
	Nevent = -1;
	evttime = 0;
}
event_data::~event_data(){
	for(map<Tomography::det_type,vector<Event*> >::iterator type_it = det_data.begin();type_it!=det_data.end();++type_it){
		while((type_it->second).size()>0){
			delete (type_it->second).back();
			(type_it->second).pop_back();
		}
	}
}
ray_data::ray_data(){
	CBevent = NULL;
}
ray_data::~ray_data(){
	if(CBevent!=NULL) delete CBevent;
}
deviation_data::deviation_data(){
	CBevent = NULL;
}
deviation_data::~deviation_data(){
	if(CBevent!=NULL) delete CBevent;
}

Task::Task(){

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
	if(new_task->is_queueable()){
		pthread_mutex_lock(&queue_mutex);
		task_queue.push(new_task);
		pthread_mutex_unlock(&queue_mutex);
	}
}
unsigned int Task::task_left(){
	return task_queue.size();
}

template<typename T>
Typed_Task<T>::Typed_Task(){
	pthread_mutex_init(&data_queue_mutex, NULL);
	data_treated = 0;
	Display_Thread::get_instance()->register_task(this);
}
template<typename T>
Typed_Task<T>::~Typed_Task(){
	pthread_mutex_lock(&data_queue_mutex);
	while(!(data_queue.empty())){
		T * next_data = data_queue.front();
		data_queue.pop();
		delete next_data;
	}
	pthread_mutex_unlock(&data_queue_mutex);
	pthread_mutex_destroy(&data_queue_mutex);
	Display_Thread::get_instance()->unregister_task(this);
}
template<typename T>
void Typed_Task<T>::push_next_data(T * next_data){
	pthread_mutex_lock(&data_queue_mutex);
	data_queue.push(next_data);
	data_treated++;
	pthread_mutex_unlock(&data_queue_mutex);
}
template<typename T>
T * Typed_Task<T>::get_next_data(){
	T * next_data;
	pthread_mutex_lock(&data_queue_mutex);
	next_data = data_queue.front();
	data_queue.pop();
	pthread_mutex_unlock(&data_queue_mutex);
	return next_data;
}
template<typename T>
bool Typed_Task<T>::is_queue_empty() const{
	return data_queue.empty();
}
template<typename T>
bool Typed_Task<T>::can_exec() const{
	return !(data_queue.empty());
}
template<typename T>
unsigned int Typed_Task<T>::get_queue_size() const{
	return data_queue.size();
}
template<typename T>
string Typed_Task<T>::print_count() const{
	ostringstream outstring;
	outstring << right << setw(10) << data_treated << " - " << left << setw(6) << data_queue.size();
	return outstring.str();
}
template<>
string Typed_Task<raw_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "raw evt";
	return outstring.str();
}
template<>
string Typed_Task<ped_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "ped evt";
	return outstring.str();
}
template<>
string Typed_Task<corr_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "corr evt";
	return outstring.str();
}
template<>
string Typed_Task<event_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "demux evt";
	return outstring.str();
}
template<>
string Typed_Task<ray_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "abs tracked evt";
	return outstring.str();
}
template<>
string Typed_Task<deviation_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "dev tracked evt";
	return outstring.str();
}

template class Typed_Task<raw_data>;
template class Typed_Task<ped_data>;
template class Typed_Task<corr_data>;
template class Typed_Task<event_data>;
template class Typed_Task<ray_data>;
template class Typed_Task<deviation_data>;

template<typename T>
Buffer_Task<T>::Buffer_Task(): Typed_Task<T>(){

}
template<typename T>
Buffer_Task<T>::~Buffer_Task(){

}
template<typename T>
bool Buffer_Task<T>::do_task(){
	return true;
}
template<typename T>
bool Buffer_Task<T>::can_exec() const{
	return true;
}
template<typename T>
void Buffer_Task<T>::update_task_list() const{

}
template<typename T>
T * Buffer_Task<T>::fetch_data(){
	return this->get_next_data();
}
template<typename T>
bool Buffer_Task<T>::can_fetch_data() const{
	return !(this->is_queue_empty());
}
template<typename T>
bool Buffer_Task<T>::is_queueable() const{
	return false;
}

template class Buffer_Task<raw_data>;
template class Buffer_Task<ped_data>;
template class Buffer_Task<corr_data>;
template class Buffer_Task<event_data>;
template class Buffer_Task<ray_data>;
template class Buffer_Task<deviation_data>;

template<typename T>
Output_Task<T>::Output_Task(): Typed_Task<T>(){
	pthread_mutex_init(&IO_mutex, NULL);
	next_task = NULL;
}
template<typename T>
Output_Task<T>::Output_Task(Typed_Task<T> * next_task_): Typed_Task<T>(){
	pthread_mutex_init(&IO_mutex, NULL);
	next_task = next_task_;
}
template<typename T>
Output_Task<T>::~Output_Task(){
	pthread_mutex_destroy(&IO_mutex);
	if(next_task != NULL) delete next_task;
}
template<>
string Output_Task<raw_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "raw evt (w)";
	return outstring.str();
}
template<>
string Output_Task<ped_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "ped evt (w)";
	return outstring.str();
}
template<>
string Output_Task<corr_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "corr evt (w)";
	return outstring.str();
}
template<>
string Output_Task<event_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "demux evt (w)";
	return outstring.str();
}
template<>
string Output_Task<ray_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "abs tracked evt (w)";
	return outstring.str();
}
template<>
string Output_Task<deviation_data>::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "dev tracked evt (w)";
	return outstring.str();
}
template<typename T>
bool Output_Task<T>::is_queueable() const{
	return false;
}
template<typename T>
void Output_Task<T>::update_task_list() const{
	if(next_task!=NULL) Task::add_task(next_task);
}

template class Output_Task<raw_data>;
template class Output_Task<ped_data>;
template class Output_Task<corr_data>;
template class Output_Task<event_data>;
template class Output_Task<ray_data>;
template class Output_Task<deviation_data>;

Input_Task::Input_Task(long max_event_): Task(){
	pthread_mutex_init(&IO_mutex, NULL);
	max_event = max_event_;
}
Input_Task::~Input_Task(){
	pthread_mutex_destroy(&IO_mutex);
}
string Input_Task::init_count() const{
	return "";
}
string Input_Task::print_count() const{
	return "";
}
bool Input_Task::is_queueable() const{
	return false;
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
		int result = pthread_create(&id,&attr,runThread,this);
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
			if(wait_time>50000){
				working = false;
				*(Display_Thread::get_instance()) << "worker thread timeout" << endl;
			}
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
		if(current_task->is_saturated()){
			usleep(1000);
			continue;
		}
		if(current_task->can_exec()){
			wait_time = 0;
			while(!(current_task->do_task())){
				usleep(1000);
				wait_time++;
				if(wait_time>100){
					working = false;
					*(Display_Thread::get_instance()) << "reader thread execution timeout" << endl;
					break;
				}
			}
			if(working) current_task->update_task_list();
		}
		else{
			usleep(1000);
			wait_time++;
			if(wait_time>10){
				working = false;
				*(Display_Thread::get_instance()) << "reader thread check timeout" << endl;
				break;
			}
		}
	}
	return 0;
}
void Reader_Thread::pre_stop(){
	working = false;
}

template<typename T>
Writer_Thread::Writer_Thread(Output_Task<T> * current_task_): Thread(){
	working = false;
	current_task = current_task_;
}
template Writer_Thread::Writer_Thread(Output_Task<event_data> * current_task_);
template Writer_Thread::Writer_Thread(Output_Task<raw_data> * current_task_);
template Writer_Thread::Writer_Thread(Output_Task<ped_data> * current_task_);
template Writer_Thread::Writer_Thread(Output_Task<corr_data> * current_task_);
template Writer_Thread::Writer_Thread(Output_Task<ray_data> * current_task_);
template Writer_Thread::Writer_Thread(Output_Task<deviation_data> * current_task_);

Writer_Thread::~Writer_Thread(){
	working = false;
	current_task = NULL;
}
bool Writer_Thread::is_working() const{
	return working;
}
void * Writer_Thread::run(){
	working = current_task!=NULL;
	unsigned int wait_time = 0;
	while(working){
		if(current_task->can_exec()){
			wait_time = 0;
			while(!(current_task->do_task())){
				usleep(1000);
				wait_time++;
				if(wait_time>100){
					working = false;
					*(Display_Thread::get_instance()) << "writer thread execution timeout" << endl;
					break;
				}
			}
			if(working) current_task->update_task_list();
		}
		else{
			usleep(1000);
			wait_time++;
			if(wait_time>10){
				working = false;
				*(Display_Thread::get_instance()) << "writer thread check timeout" << endl;
				break;
			}
		}
	}
	return 0;
}
void Writer_Thread::pre_stop(){
	working = false;
}

Display_Thread * Display_Thread::singleton_instance = 0;
Display_Thread * Display_Thread::get_instance(){
	if(!singleton_instance){
		singleton_instance = new Display_Thread();
		singleton_instance->start();
	}
	return singleton_instance;
}
void Display_Thread::Quit(){
	if(singleton_instance) delete singleton_instance;
	singleton_instance = 0;
}

Display_Thread::Display_Thread(): Thread(){
	is_counting = false;
}
Display_Thread::Display_Thread(string log_file_name): Thread(), ostringstream(){
	log_file.open(log_file_name.c_str());
	is_counting = false;
}
Display_Thread::~Display_Thread(){
	stop();
	if(log_file.is_open() && log_file.good()){
		log_file.flush();
		log_file.close();
	}
	cout << endl;
	singleton_instance = 0;
}
void Display_Thread::set_log_file(string log_file_name){
	if(log_file.is_open()){
		if(log_file.good()) log_file.flush();
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
	is_counting = false;
}
void * Display_Thread::run(){
	working = true;
	unsigned int delay = 0;
	while(working){
		display_text();
		if(delay>300){
			display_canvas();
			delay = 0;
		}
		usleep(10000);
		delay++;
	}
	return 0;
}
void Display_Thread::display_text(){
	ostringstream temp;
	string buffer = str();
	str("");
	if(is_counting && (!registered_task.empty())){
		if(buffer.size()>0) temp << "\n" << buffer << "\n";
		temp << "\r" << setw(19) << Task::task_left();
		vector<Task*>::reverse_iterator task_it=registered_task.rbegin();
		while(task_it!=registered_task.rend()){
			temp << "|" << (*task_it)->print_count();
			++task_it;
		}
	}
	else temp << buffer;
	cout << temp.str();
	cout.flush();
	if(log_file.is_open() && log_file.good()){
		log_file << temp.str();
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
void Display_Thread::start_count(){
	if(registered_task.empty()) return;
	*this << "\r" << setw(19) << "task count";
	vector<Task*>::reverse_iterator task_it=registered_task.rbegin();
	while(task_it!=registered_task.rend()){
		*this << "|" << (*task_it)->init_count();
		++task_it;
	}
	is_counting = true;
}
void Display_Thread::stop_count(){
	is_counting = false;
}
void Display_Thread::register_task(Task * some_task){
	registered_task.push_back(some_task);
}
void Display_Thread::unregister_task(Task * some_task){
	bool constant = true;
	vector<Task*>::iterator it = registered_task.begin();
	while(it!=registered_task.end()){
		if(some_task == *it){
			constant = false;
			it = registered_task.erase(it);
		}
		else{
			++it;
		}
	}
	if(constant){
		*this << "Task was not registered, counld not unregister" << endl;
	}
}
