#define tomography_cpp
#include "tomography.h"
#include "detector.h"
#include "MT_tomography.h"

#include <sstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <utility>
#include <TROOT.h>
#include <TCanvas.h>
#include <TThread.h>
#include <Math/MinimizerOptions.h>
#include <Fit/FitConfig.h>
#include <TVirtualFitter.h>
#include <TError.h>

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
	os << Tomography::Static_Detector[det]->Name();
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

Tomography::Tomography(){
	cout << "Warning ! Dummy Tomography instanciated" << endl;
	root_interpreter = 0;
	config_tree = ptree();
	Nsample = -1;
	XY_size = 0;
	SampleMin = -1;
	SampleMax = -1;
	sigma = 0;
	first_down_layer = -1;
	TOTCut = -1;
	noise_RMS = ADC_max+1;
	chisquare_threshold = 0;
	live_graphic_display = false;
	is_batch = true;
	thread_number = 0;
	thread_number = 0;
	raw_data_treated = 0;
	ped_data_treated = 0;
	corr_data_treated = 0;
	event_treated = 0;
	ray_treated = 0;
	deviation_treated = 0;
	TThread::Initialize();
	gErrorIgnoreLevel = 1001;
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	ROOT::Fit::FitConfig::SetDefaultMinimizer("Minuit2");
	TVirtualFitter::SetDefaultFitter("Minuit2");
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
	first_down_layer = config_tree.get<int>("first_down_layer");
	TOTCut = config_tree.get<int>("TOTCut");
	noise_RMS = config_tree.get<double>("noise_RMS");
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
	raw_data_treated = 0;
	ped_data_treated = 0;
	corr_data_treated = 0;
	event_treated = 0;
	ray_treated = 0;
	deviation_treated = 0;
	TThread::Initialize();
	gErrorIgnoreLevel = 1001;
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	ROOT::Fit::FitConfig::SetDefaultMinimizer("Minuit2");
	TVirtualFitter::SetDefaultFitter("Minuit2");
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
	cout << "Terminating Tomography !" << endl;
}
void Tomography::Quit(){
	Display_Thread::Quit();
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
bool Tomography::get_is_up(int layer) const{
	return (layer < first_down_layer);
}
double Tomography::get_chisquare_threshold() const{
	return chisquare_threshold;
}
double Tomography::get_sigma() const{
	return sigma;
}
double Tomography::get_noise_RMS() const{
	return noise_RMS;
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
	Display_Thread::Quit();
	if(root_interpreter) root_interpreter->Run(true);
}

unsigned short Tomography::get_thread_number() const{
	return thread_number;
}
/*
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
*/
