#define read_live_task_cpp

//#include <sys/msg.h>
extern "C" {
#include "Pipes.h"
}
#include "task/read_live_task.h"

#include "ElecReader.h"

#include <sys/time.h>

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::left;
using std::right;
using std::setw;
/*
using std::hex;
using std::dec;
using std::noshowbase;
using std::showbase;
*/

Read_Live_Task::Read_Live_Task(string pipe_name): Input_Task(-1){
	pthread_cond_init(&queue_cond,NULL);
	pipe_ptr = 0;
	queue_id = Pipe_Create(&pipe_ptr,const_cast<char*>(pipe_name.c_str()), DEF_Pipe_Type_Rd);
	if(queue_id < 0){
		cout << "failed to create message queue with name : " << pipe_name << endl;
	}
	else{
		cout << "created message queue (id : " << queue_id << ") with name " << pipe_name << endl;
	}
	status = 0;
	data_count = 0;
	Display_Thread::get_instance()->register_task(this);
}
Read_Live_Task::~Read_Live_Task(){
	pthread_cond_destroy(&queue_cond);
	if(pipe_ptr) Pipe_Delete(&pipe_ptr);
	pipe_ptr = 0;
	Display_Thread::get_instance()->unregister_task(this);
}
bool Read_Live_Task::do_task(){
	data_message * current_data = new data_message();
	//cout << "###" <<  current_data << " -> " << current_data+1 << " | " << ((unsigned long)current_data)+sizeof(*current_data) << endl;
	//cout << "trying to read message..." << endl;
	//if(msgrcv(queue_id,current_data,sizeof(current_data->data),0,0) == -1) return false;
	int rd_size;
	if(Pipe_Read(pipe_ptr,(char*)current_data, &rd_size) < 0){
		cout << "error while reading pipe" << endl;
	}
	//cout << "read message of type " << current_data->mtype << endl;
	//cout << "###" <<  current_data << " -> " << current_data+1 << " | " << ((unsigned long)current_data)+sizeof(*current_data) << " | " << ((unsigned long)current_data)+rd_size << endl;
	if(current_data->mtype == 1){
		//cout << "receiving message" << endl;
		pthread_mutex_lock(&IO_mutex);
		data_queue.push(current_data);
		pthread_cond_signal(&queue_cond);
		pthread_mutex_unlock(&IO_mutex);
		data_count++;
		//cout << "message recieved" << endl;
		return true;
	}
	else{
		//cout << "recieved control message : " << hex << showbase << current_data->mtype << dec << noshowbase << endl;
		status |= (current_data->mtype);
		delete current_data;
		return true;
	}
	delete current_data;
	return false;
}
bool Read_Live_Task::can_exec() const{
	//cout << "current status : " << hex << showbase << status << dec << noshowbase << endl;
	//if(status & 0x8000) cout << "recieved quit message" << endl;
	return (!(status & 0x8000));
}
void Read_Live_Task::update_task_list() const{

}
data_message * Read_Live_Task::wait_new_data(){
	//cout << "---" << endl;
	struct timespec wait_time;
	wait_time.tv_nsec=0;
	struct timeval now;
	data_message * current_data;
	//cout << "###" << endl;
	pthread_mutex_lock(&IO_mutex);
	//cout << "requesting message" << endl;
	while(Tomography::get_instance()->get_can_continue() && data_queue.empty() && can_exec()){
		//cout << "waiting message...";
		gettimeofday(&now,NULL);
		wait_time.tv_sec = now.tv_sec+1;
		pthread_cond_timedwait(&queue_cond,&IO_mutex,&wait_time);
		//cout << "!" << endl;
	}
	if(data_queue.empty()) current_data = new data_message();
	else{
		current_data = data_queue.front();
		data_queue.pop();
	}
	pthread_mutex_unlock(&IO_mutex);
	return current_data;
}
data_message * Read_Live_Task::get_next_data(){
	data_message * current_data;
	pthread_mutex_lock(&IO_mutex);
	current_data = data_queue.front();
	data_queue.pop();
	pthread_mutex_unlock(&IO_mutex);
	return current_data;
}
bool Read_Live_Task::has_new_data() const{
	//cout << "size : " << data_queue.size() << endl;
	return !(data_queue.empty());
}
int Read_Live_Task::get_status() const{
	return status;
}
bool Read_Live_Task::is_saturated() const{
	return false;
}

string Read_Live_Task::init_count() const{
	ostringstream outstring;
	outstring << left << setw(19) << "messages read";
	return outstring.str();
}
string Read_Live_Task::print_count() const{
	ostringstream outstring;
	outstring << right << setw(10) << data_count << " - " << left << setw(6) << data_queue.size();
	return outstring.str();
}