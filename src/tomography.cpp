#define tomography_cpp
#include "tomography.h"
#include "detector.h"

#include <sstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <utility>
#include <TROOT.h>
#include <TCanvas.h>

#include <boost/foreach.hpp>
#include <boost/property_tree/json_parser.hpp>

using std::ostringstream;
using std::setw;
using std::left;
using std::right;
using std::setfill;
using std::string;
using std::pair;

/*
bool Tomography::get_instance()->get_can_continue() = true;
bool Tomography::get_instance()->get_is_batch() = gROOT->IsBatch();
bool Tomography::get_instance()->get_live_graphic_display() = !Tomography::get_instance()->get_is_batch();
*/
Tomography* Tomography::singleton_instance = 0;
const string Tomography::DreamExt = "fdf";
const string Tomography::FeminosExt = "aqs";

ostream& operator<<(ostream& os, const Tomography::det_type& det){
	switch(det){
		case Tomography::CM : os << "CM"; break;
		case Tomography::MG : os << "MG"; break;
		case Tomography::MGv2 : os << "MGv2"; break;
		default : os << "unknown det";
	}
	return os;
}
ostream& operator<<(ostream& os, const Tomography::strip_type& strip){
	switch(strip){
		case Tomography::Wide : os << "Wide"; break;
		case Tomography::Thin : os << "Thin"; break;
		case Tomography::Demux : os << "Demux"; break;
		default : os << "unknown strip type";
	}
	return os;
}
ostream& operator<<(ostream& os, const Tomography::elec_type& elec){
	switch(elec){
		case Tomography::Dream : os << "Dream"; break;
		case Tomography::Feminos : os << "Feminos"; break;
		default : os << "unknown electronic type";
	}
	return os;
}

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_){
	if(map_.size()<1){
		os << "[]";
		return os;
	}
	os << "[ ";
	typename map<T,R>::const_iterator it=map_.begin();
	typename map<T,R>::const_iterator jt = map_.end();
	jt--;
	while(it!=jt){
		os << it->first << " -> " << it->second << " ; ";
		++it;
	}
	os << it->first << " -> " << it->second << " ]";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_){
	if(vec_.empty()){
		os << "{}";
		return os;
	}
	os << "{ ";
	typename vector<T>::const_iterator it = vec_.begin();
	typename vector<T>::const_iterator jt = vec_.end();
	jt--;
	while(it!=jt){
		os << *it << " ; ";
		++it;
	}
	os << *it << " }";
	return os;
}

template ostream& operator<<(ostream& os, const map<int,double>& map_);
template ostream& operator<<(ostream& os, const map<double,int>& map_);
template ostream& operator<<(ostream& os, const map<int,int>& map_);
template ostream& operator<<(ostream& os, const map<bool,map<int,int> >& map_);

Tomography::elec_type Tomography::str_to_elec(string str){
	elec_type return_value = unknown_elec;
	if(str == "dream") return_value = Dream;
	else if(str == "feminos") return_value = Feminos;
	return return_value;
}

static map<const Tomography::det_type,const Detector* const> Static_Detector_build(){
	map<const Tomography::det_type,const Detector* const> return_map;
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::CM,new CM_Detector()));
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::MG,new MG_Detector()));
	return_map.insert(pair<const Tomography::det_type,const Detector* const>(Tomography::MGv2,new MGv2_Detector()));
	return return_map;
}

map<const Tomography::det_type,const Detector* const> Tomography::Static_Detector = Static_Detector_build();
Display_Thread Tomography::TomoCout;// = Display_Thread();

Tomography::Tomography(){
	cout << "Warning ! Dummy Tomography instanciated" << endl;
	root_interpreter = 0;
	config_tree = ptree();
	Nsample = -1;
	XY_size = 0;
	SampleMin = -1;
	SampleMax = -1;
	sigma = 0;
	TOTCut = -1;
	chisquare_threshold = 0;
	live_graphic_display = false;
	is_batch = true;
	thread_number = 0;
	pthread_mutex_init(&raw_data_mutex,NULL);
	pthread_mutex_init(&ped_data_mutex,NULL);
	pthread_mutex_init(&corr_data_mutex,NULL);
	pthread_mutex_init(&event_mutex,NULL);
	pthread_mutex_init(&ray_mutex,NULL);
	pthread_mutex_init(&deviation_mutex,NULL);
	thread_number = 0;
	raw_data_treated = 0;
	ped_data_treated = 0;
	corr_data_treated = 0;
	event_treated = 0;
	ray_treated = 0;
	deviation_treated = 0;
	sigIntHandler.sa_handler = signal_handler;
	gROOT->SetBatch(true);
	can_continue = true;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);
}
Tomography::Tomography(ptree config_tree_){
	cout << "Loading Tomography !" << endl;
	config_tree = config_tree_;
	Nsample = config_tree.get<int>("Nsample");
	if(Nsample>Max_Nsample){
		cout << "Warning ! Nsample (" << Nsample << ") superior to max Nsample allowed (" << Max_Nsample << ")" << endl;
		Nsample = Max_Nsample;
	}
	XY_size = config_tree.get<double>("XY_size");
	SampleMin = config_tree.get<int>("SampleMin");
	SampleMax = config_tree.get<int>("SampleMax");
	sigma = config_tree.get<double>("sigma");
	TOTCut = config_tree.get<int>("TOTCut");
	chisquare_threshold = config_tree.get<double>("chisquare_threshold");
	if(gROOT->IsBatch()) is_batch = true;
	else is_batch = config_tree.get<bool>("batch");
	if(!is_batch){
		cout << "Starting Graphic Mode" << endl;
		live_graphic_display = config_tree.get<bool>("live_graphic_display");
		root_interpreter = new TRint("Rint",0,0,0,0,true);
	}
	else{
		cout << "Starting Batch Mode" << endl;
		gROOT->SetBatch(true);
		live_graphic_display = false;
		root_interpreter = 0;
	}
	can_continue = true;
	thread_number = config_tree.get<unsigned short>("thread_number");
	if(thread_number<1) thread_number = 1;
	pthread_mutex_init(&raw_data_mutex,NULL);
	pthread_mutex_init(&ped_data_mutex,NULL);
	pthread_mutex_init(&corr_data_mutex,NULL);
	pthread_mutex_init(&event_mutex,NULL);
	pthread_mutex_init(&ray_mutex,NULL);
	pthread_mutex_init(&deviation_mutex,NULL);
	raw_data_treated = 0;
	ped_data_treated = 0;
	corr_data_treated = 0;
	event_treated = 0;
	ray_treated = 0;
	deviation_treated = 0;
	sigIntHandler.sa_handler = signal_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);
}
Tomography::~Tomography(){
	save_canvases();
	if(root_interpreter) delete root_interpreter;
	is_batch = true;
	live_graphic_display = false;
	singleton_instance = 0;
	can_continue = false;
	thread_number = 0;
	pthread_mutex_destroy(&raw_data_mutex);
	pthread_mutex_destroy(&ped_data_mutex);
	pthread_mutex_destroy(&corr_data_mutex);
	pthread_mutex_destroy(&event_mutex);
	pthread_mutex_destroy(&ray_mutex);
	pthread_mutex_destroy(&deviation_mutex);
	cout << "Terminating Tomography !" << endl;
}
void Tomography::Quit(){
	if(singleton_instance != 0){
		delete singleton_instance;
		singleton_instance = 0;
	}
}
Tomography * Tomography::get_instance(){
	if(!singleton_instance){
		singleton_instance = new Tomography();
	}
	return singleton_instance;
}
void Tomography::Init(ptree config_tree_){
	if(!singleton_instance){
		singleton_instance = new Tomography(config_tree_);
	}
}
void Tomography::Init(string config_tree_file){
	if(!singleton_instance){
		ptree config_tree_;
		read_json(config_tree_file,config_tree_);
		singleton_instance = new Tomography(config_tree_);
	}
}
void Tomography::signal_handler(int s){
	cout << "\nCaught signal " << s << endl;
	cout << endl;
	get_instance()->can_continue = false;
}
int Tomography::get_Nsample() const{
	return Nsample;
}
double Tomography::get_XY_size() const{
	return XY_size;
}
int Tomography::get_SampleMin() const{
	return SampleMin;
}
int Tomography::get_SampleMax() const{
	return SampleMax;
}
int Tomography::get_TOTCut() const{
	return TOTCut;
}
double Tomography::get_chisquare_threshold() const{
	return chisquare_threshold;
}
double Tomography::get_sigma() const{
	return sigma;
}
bool Tomography::get_live_graphic_display() const{
	return live_graphic_display;
}
bool Tomography::get_is_batch() const{
	return is_batch;
}
bool Tomography::get_can_continue() const{
	return can_continue;
}
void Tomography::save_canvases(){
	time_t current_time = time(NULL);
	char buffer[100];
	strftime(buffer,100,"%y%m%d_%HH%M",localtime(&current_time));
	string base_name = config_tree.get<string>("metadata");
	if(base_name.size() == 0) return;
	base_name += "_";
	base_name += buffer;
	base_name += "_";
	for(int i=0;i<gROOT->GetListOfCanvases()->GetSize();i++){
		TCanvas * current_canvas = dynamic_cast<TCanvas*>(gROOT->GetListOfCanvases()->At(i));
		current_canvas->SaveAs((base_name + current_canvas->GetName() + ".png").c_str());
		current_canvas->SaveAs((base_name + current_canvas->GetName() + ".C").c_str());
	}
	if(gROOT->GetListOfCanvases()->GetSize() > 0) write_json(base_name + "config.json",config_tree);
}
void Tomography::Run(){
	TomoCout.stop();
	if(root_interpreter) root_interpreter->Run(true);
}

struct raw_data Tomography::get_next_raw_data(){
	pthread_mutex_lock(&raw_data_mutex);
	struct raw_data next_data = raw_data_queue.front();
	raw_data_queue.pop();
	pthread_mutex_unlock(&raw_data_mutex);
	return next_data;
}
void Tomography::push_next_raw_data(struct raw_data new_data){
	pthread_mutex_lock(&raw_data_mutex);
	raw_data_queue.push(new_data);
	pthread_mutex_unlock(&raw_data_mutex);
	raw_data_treated++;
}
bool Tomography::is_raw_data_empty() const{
	return raw_data_queue.empty();
}
struct ped_data Tomography::get_next_ped_data(){
	pthread_mutex_lock(&ped_data_mutex);
	struct ped_data next_data = ped_data_queue.front();
	ped_data_queue.pop();
	pthread_mutex_unlock(&ped_data_mutex);
	return next_data;
}
void Tomography::push_next_ped_data(struct ped_data new_data){
	pthread_mutex_lock(&ped_data_mutex);
	ped_data_queue.push(new_data);
	pthread_mutex_unlock(&ped_data_mutex);
	ped_data_treated++;
}
bool Tomography::is_ped_data_empty() const{
	return ped_data_queue.empty();
}
struct ped_data Tomography::get_next_corr_data(){
	pthread_mutex_lock(&corr_data_mutex);
	struct ped_data next_data = corr_data_queue.front();
	corr_data_queue.pop();
	pthread_mutex_unlock(&corr_data_mutex);
	return next_data;
}
void Tomography::push_next_corr_data(struct ped_data new_data){
	pthread_mutex_lock(&corr_data_mutex);
	corr_data_queue.push(new_data);
	pthread_mutex_unlock(&corr_data_mutex);
	corr_data_treated++;
}
bool Tomography::is_corr_data_empty() const{
	return corr_data_queue.empty();
}
struct event_data Tomography::get_next_event_data(){
	pthread_mutex_lock(&event_mutex);
	struct event_data next_data = event_queue.front();
	event_queue.pop();
	pthread_mutex_unlock(&event_mutex);
	return next_data;
}
void Tomography::push_next_event_data(struct event_data new_data){
	pthread_mutex_lock(&event_mutex);
	event_queue.push(new_data);
	pthread_mutex_unlock(&event_mutex);
	event_treated++;
}
bool Tomography::is_event_data_empty() const{
	return event_queue.empty();
}
struct ray_data Tomography::get_next_ray_data(){
	pthread_mutex_lock(&ray_mutex);
	struct ray_data next_data = ray_queue.front();
	ray_queue.pop();
	pthread_mutex_unlock(&ray_mutex);
	return next_data;
}
void Tomography::push_next_ray_data(struct ray_data new_data){
	pthread_mutex_lock(&ray_mutex);
	ray_queue.push(new_data);
	pthread_mutex_unlock(&ray_mutex);
	ray_treated++;
}
bool Tomography::is_ray_data_empty() const{
	return ray_queue.empty();
}
struct deviation_data Tomography::get_next_deviation_data(){
	pthread_mutex_lock(&deviation_mutex);
	struct deviation_data next_data = deviation_queue.front();
	deviation_queue.pop();
	pthread_mutex_unlock(&deviation_mutex);
	return next_data;
}
void Tomography::push_next_deviation_data(struct deviation_data new_data){
	pthread_mutex_lock(&deviation_mutex);
	deviation_queue.push(new_data);
	pthread_mutex_unlock(&deviation_mutex);
	deviation_treated++;
}
bool Tomography::is_deviation_data_empty() const{
	return deviation_queue.empty();
}
unsigned short Tomography::get_thread_number() const{
	return thread_number;
}
string Tomography::init_count() const{
	ostringstream outstring;
	outstring << left << setw(21) << "raw event" << "|" << setw(21) << "ped event" << "|" << setw(21) << "corr event" << "|" << setw(21) << "demux event" << "|" << setw(21) << "tracked event - abs" << "|" << setw(21) << "tracked event - dev";
	return outstring.str();
}
string Tomography::print_count() const{
	ostringstream outstring;
	outstring << right << setw(11) << raw_data_treated << " - " << left << setw(7) << raw_data_queue.size() << "|" << right << setw(11) << ped_data_treated << " - " << left << setw(7) << ped_data_queue.size() << "|" << right << setw(11) << corr_data_treated << " - " << left << setw(7) << corr_data_queue.size() << "|" << right << setw(11) << event_treated << " - " << left << setw(7) << event_queue.size() << "|" << right << setw(11) << ray_treated << " - " << left << setw(7) << ray_queue.size() << "|" << right << setw(11) << deviation_treated << " - " << left << setw(7) << deviation_queue.size();
	return outstring.str();
}