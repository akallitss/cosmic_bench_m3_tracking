#define livedisplay_cpp

#include "liveDisplay.h"
#include "datareader.h"
#include "event.h"

#include <unistd.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <limits>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>

using std::ifstream;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::cout;
using std::endl;
using std::flush;
using std::bitset;
using std::numeric_limits;

using boost::property_tree::ptree;

//inotify stuff
#define MAX_EVENTS 1024
#define LEN_NAME 16
#define EVENT_SIZE  ( sizeof (struct inotify_event) )
#define BUF_LEN     ( MAX_EVENTS * ( EVENT_SIZE + LEN_NAME ))

liveDisplay::liveDisplay(): CosmicBench(){
	filenames.clear();
	max_event = -1;
	inotify_descriptor = -1;
	file_descriptor = -1;
	inotify_started = false;
	current_file = "";
	MG_N = 0;
	CM_N = 0;
	det_type_by_asic.clear();
	det_n_by_asic.clear();
}

liveDisplay::~liveDisplay(){

}

liveDisplay::liveDisplay(string config_file){
	ptree config_tree;
	read_json(config_file, config_tree);
	filenames.clear();
	max_event = config_tree.get<long>("max_event");
	inotify_descriptor = -1;
	file_descriptor = -1;
	inotify_started = false;
	current_file = "";
	CosmicBench::Init(config_tree);
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::CM;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("cm_n");
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::MG;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("mg_n");
	}
	if(CM_N!=0) cout << "warning, CosMultis are not fully supported !" << endl;
	electronic_type = Tomography::str_to_elec(config_tree.get<string>("electronic_type"));
	data_file_basename = config_tree.get<string>("data_file_basename");
	ifstream pedFile((config_tree.get<string>("Ped")).c_str());
	/*
	Pedestal_CM = new float[CM_N][DataReader::Nchannel_CM];
	Pedestal_MG = new float[MG_N][DataReader::Nchannel_MG];
	for(int i=0;i<CM_N;i++){
		for(int j=0;j<DataReader::Nchannel_CM;j++){
			pedFile >> i >> j >> Pedestal_CM[i][j];
		}
	}
	for(int i=0;i<MG_N;i++){
		for(int j=0;j<DataReader::Nchannel_MG;j++){
			pedFile >> i >> j >> Pedestal_MG[i][j];
		}
	}
	*/
	Pedestal[Tomography::CM] = vector<vector<float> >(CM_N,vector<float>(DataReader::Nchannel_CM,0));
	Pedestal[Tomography::MG] = vector<vector<float> >(MG_N,vector<float>(DataReader::Nchannel_MG,0));
	for(int i=0;i<CM_N;i++){
		for(int j=0;j<DataReader::Nchannel_CM;j++){
			pedFile >> i >> j >> Pedestal[Tomography::CM][i][j];
		}
	}
	for(int i=0;i<MG_N;i++){
		for(int j=0;j<DataReader::Nchannel_MG;j++){
			pedFile >> i >> j >> Pedestal[Tomography::MG][i][j];
		}
	}
	pedFile.close();
	max_file_size = config_tree.get<long>("max_file_size");
	add_files(config_tree.get<int>("data_file_first"),config_tree.get<int>("data_file_last"));
}

void liveDisplay::add_file(string filename){
	filenames.push_back(filename);
}
void liveDisplay::add_files(int first,int last){
	string extension = "";
	if(electronic_type == Tomography::Feminos) extension = ".aqs";
	else if(electronic_type == Tomography::Dream) extension = "_01.fdf";
	else extension = ".txt";
	for(int i=first;i<=last;i++){
		ostringstream name;
		name << data_file_basename << setfill('0') << setw(3) << i << extension;
		filenames.push_back(name.str());
	}
}
int liveDisplay::start_inotify(string filename){
	inotify_descriptor = inotify_init();
	current_file = filename;
	file_descriptor = inotify_add_watch(inotify_descriptor,current_file.c_str(), IN_ALL_EVENTS);
	if(inotify_descriptor>=0 && file_descriptor>=0) inotify_started = true;
	if(inotify_descriptor<0) return inotify_descriptor;
	return file_descriptor;
}
int liveDisplay::pause_inotify(){
	if(!inotify_started) return -1;
	inotify_rm_watch(inotify_descriptor,file_descriptor);
	close(inotify_descriptor);
	return 1;
}
int liveDisplay::resume_inotify(){
	if(!inotify_started) return -1;
	inotify_descriptor = inotify_init();
	file_descriptor = inotify_add_watch(inotify_descriptor,current_file.c_str(), IN_ALL_EVENTS);
	if(inotify_descriptor<0) return inotify_descriptor;
	return file_descriptor;
}
int liveDisplay::stop_inotify(){
	if(!inotify_started) return -1;
	inotify_rm_watch(inotify_descriptor,file_descriptor);
	close(inotify_descriptor);
	inotify_started = false;
	return 1;
}

unsigned int liveDisplay::read_inotify(){
	unsigned int global_mask = 0x00000000;
	if(!inotify_started) return global_mask;
	char buffer[BUF_LEN];
	int length = read(inotify_descriptor,buffer, BUF_LEN);
	int i = 0;
	while(i<length){
		struct inotify_event * event = (struct inotify_event*) &buffer[i];
		global_mask |= event->mask;
		i+= EVENT_SIZE + event->len;
	}
	return global_mask;
}

void liveDisplay::flux_map(double z){
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * cMap = new TCanvas("flux_map");
	TCanvas * cDisplay = new TCanvas("event_display");
	TCanvas * cStats = new TCanvas("stats");
	int bin_n = 100;
	double z_max = numeric_limits<double>::min();
	double z_min = numeric_limits<double>::max();
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		double current_z = (*det_it)->get_z();
		if(current_z>z_max) z_max = current_z;
		if(current_z<z_min) z_min = current_z;
	}
	double x_min = -Tomography::get_instance()->get_XY_size()/2.;
	double x_max = Tomography::get_instance()->get_XY_size()/2.;
	if(z>z_max){
		x_min -= Tomography::get_instance()->get_XY_size()*(z - z_max)/(z_max - z_min);
		x_max = - x_min;
	}
	else if(z<z_min){
		x_min -= Tomography::get_instance()->get_XY_size()*(z_min - z)/(z_max - z_min);
		x_max = - x_min;
	}
	double width = x_max - x_min;
	x_min -= 0.05*width;
	x_max += 0.05*width;
	map<Tomography::det_type,int> detector_div;
	detector_div[Tomography::CM] = 2;
	detector_div[Tomography::MG] = 2;
	map<Tomography::det_type,int> Nchannel;
	Nchannel[Tomography::CM] = DataReader::Nchannel_CM;
	Nchannel[Tomography::MG] = DataReader::Nchannel_MG;
	long eventReconstructed = 0;
	long eventSuitable = 0;
	long processed = 0;
	double evttime = 0;
	time_t last_time = time(NULL);
	TH2D * flux_map_h = new TH2D("flux_map","flux_map",bin_n,x_min,x_max,bin_n,x_min,x_max);
	flux_map_h->SetStats(0);
	TPaveText * stat_text = new TPaveText(0,0,1,1);
	DataReader * current_data_reader = NULL;
	long event_nb = 0;
	if(electronic_type == Tomography::Feminos){
		current_data_reader = new FeminosDataReader("live",det_type_by_asic,det_n_by_asic);
	}
	else if(electronic_type == Tomography::Dream){
		current_data_reader = new DreamDataReader("live",det_type_by_asic,det_n_by_asic);
	}
	else return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for(vector<string>::iterator filename_it = filenames.begin();filename_it!=filenames.end() && Tomography::get_instance()->get_can_continue();++filename_it){
		cout << *filename_it << endl;
		ifstream data_file(filename_it->c_str(),ifstream::binary);
		bool is_open = data_file.is_open();
		if(!is_open){
			cout << "can't open file ! Switch to next one" << endl;
			continue;
		}
		start_inotify(*filename_it);
		int current_pos = data_file.tellg();
		while(is_open && Tomography::get_instance()->get_can_continue()){
			if(max_event>0 && processed>max_event) break;
			current_pos = data_file.tellg();
			data_file.seekg(0,data_file.end);
			bool is_complete = !(data_file.tellg()<max_file_size);
			data_file.seekg(current_pos, data_file.beg);
			unsigned int read_mask = 0x00000000;
			if(!(current_pos>=max_file_size || is_complete)) read_mask = read_inotify();
			else read_mask = IN_MODIFY | IN_CLOSE;
			if((read_mask & IN_CLOSE) || is_complete) is_open = false;
			if(!(read_mask & IN_MODIFY)) continue;
			while(data_file.good() && Tomography::get_instance()->get_can_continue()){
				current_pos = data_file.tellg();
				map<Tomography::det_type,vector<vector<vector<double> > > > event_ampl = current_data_reader->read_event(&data_file,event_nb,evttime);
				if(event_ampl.size()==0){
					data_file.clear();
					data_file.seekg(current_pos, data_file.beg);
					break;
				}
				/*
				//pedestal sub
				for(unsigned int i=0;i<MG_N;i++){
					for(int j=0;j<DataReader::Nchannel_MG;j++){
						for(int k=0;k<DataReader::Nsample;k++){
							event_ampl["MG"][i][j][k] -= Pedestal_MG[i][j];
						}
					}
				}
				for(unsigned int i=0;i<CM_N;i++){
					for(int j=0;j<DataReader::Nchannel_CM;j++){
						for(int k=0;k<DataReader::Nsample;k++){
							event_ampl["CM"][i][j][k] -= Pedestal_CM[i][j];
						}
					}
				}
				//common noise sub
				for(unsigned int i=0;i<MG_N;i++){
					for(int k=0;k<DataReader::Nsample;k++){
						for(int det_div=0;det_div<detector_div_MG;det_div++){
							int strip_nb = DataReader::Nchannel_MG/detector_div_MG;
							int strip_offset = det_div*strip_nb;
							vector<float> current_sample(strip_nb,0);
							for(int j=0;j<strip_nb;j++){
								current_sample[j] = event_ampl["MG"][i][j+strip_offset][k];
							}
							sort(current_sample.begin(),current_sample.end());
							float median = current_sample[strip_nb/2];
							for(int j=0;j<strip_nb;j++){
								event_ampl["MG"][i][j+strip_offset][k] -= median;
							}
						}
					}
				}
				for(unsigned int i=0;i<CM_N;i++){
					for(int k=0;k<DataReader::Nsample;k++){
						for(int det_div=0;det_div<detector_div_CM;det_div++){
							int strip_nb = DataReader::Nchannel_CM/detector_div_CM;
							int strip_offset = det_div*strip_nb;
							vector<float> current_sample(strip_nb,0);
							for(int j=0;j<strip_nb;j++){
								current_sample[j] = event_ampl["CM"][i][j+strip_offset][k];
							}
							sort(current_sample.begin(),current_sample.end());
							float median = current_sample[strip_nb/2];
							for(int j=0;j<strip_nb;j++){
								event_ampl["CM"][i][j+strip_offset][k] -= median;
							}
						}
					}
				}
				*/
				//Ped + CMN sub
				/*
				for(int i=0;i<MG_N;i++){
					for(int k=0;k<DataReader::Nsample;k++){
						for(int det_div=0;det_div<detector_div_MG;det_div++){
							int strip_nb = DataReader::Nchannel_MG/detector_div_MG;
							int strip_offset = det_div*strip_nb;
							vector<float> current_sample(strip_nb,0);
							for(int j=0;j<strip_nb;j++){
								current_sample[j] = event_ampl["MG"][i][j+strip_offset][k] - Pedestal_MG[i][j+strip_offset];
							}
							sort(current_sample.begin(),current_sample.end());
							float median = current_sample[strip_nb/2];
							for(int j=0;j<strip_nb;j++){
								event_ampl["MG"][i][j+strip_offset][k] -= median + Pedestal_MG[i][j+strip_offset];
							}
						}
					}
				}
				for(int i=0;i<CM_N;i++){
					for(int k=0;k<DataReader::Nsample;k++){
						for(int det_div=0;det_div<detector_div_CM;det_div++){
							int strip_nb = DataReader::Nchannel_CM/detector_div_CM;
							int strip_offset = det_div*strip_nb;
							vector<float> current_sample(strip_nb,0);
							for(int j=0;j<strip_nb;j++){
								current_sample[j] = event_ampl["CM"][i][j+strip_offset][k] - Pedestal_CM[i][j+strip_offset];
							}
							sort(current_sample.begin(),current_sample.end());
							float median = current_sample[strip_nb/2];
							for(int j=0;j<strip_nb;j++){
								event_ampl["CM"][i][j+strip_offset][k] -= median + Pedestal_CM[i][j+strip_offset];
							}
						}
					}
				}
				*/
				map<Tomography::det_type,vector<vector<float> > >::iterator ped_it = Pedestal.begin();
				for(map<Tomography::det_type,vector<vector<vector<double> > > >::iterator it = event_ampl.begin();it!=event_ampl.end();++it){
					vector<vector<float> >::iterator ped_jt = (ped_it->second).begin();
					for(vector<vector<vector<double> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
						for(int k=0;k<DataReader::Nsample;k++){
							for(int det_div=0;det_div<detector_div[it->first];det_div++){
								int strip_nb = Nchannel[it->first]/detector_div[it->first] + Nchannel[it->first]%detector_div[it->first];
								int strip_offset = det_div*strip_nb;
								vector<float> current_sample(strip_nb,0);
								for(int j=0;(j<strip_nb && (j+strip_offset)<Nchannel[it->first]);j++){
									current_sample[j] = (*jt)[j+strip_offset][k] - (*ped_jt)[j+strip_offset];
								}
								sort(current_sample.begin(),current_sample.end());
								float median = current_sample[strip_nb/2];
								for(int j=0;(j<strip_nb && (j+strip_offset)<Nchannel[it->first]);j++){
									(*jt)[j+strip_offset][k] -= median + (*ped_jt)[j+strip_offset];
								}
							}
						}
						++ped_jt;
					}
					++ped_it;
				}
				vector<Event*> events;
				for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
					if((*it)->get_type() == Tomography::MG){
						MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
						events.push_back(new MG_Event(*current_det,event_ampl[Tomography::MG][current_det->get_mg_n_in_tree()],use_srf,event_nb));
						(events.back())->MultiCluster();
						(events.back())->do_cuts();
					}
					else if((*it)->get_type() == Tomography::CM){
						CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
						events.push_back(new CM_Event(*current_det,event_ampl[Tomography::CM][current_det->get_cm_n_in_tree()],use_srf,event_nb));
						(events.back())->MultiCluster();
						(events.back())->do_cuts();
					}
				}
				CosmicBenchEvent * current_full_event = new CosmicBenchEvent(this,events);
				double chisquare_threshold = 100;
				vector<Ray> currentRays = current_full_event->get_absorption_rays(chisquare_threshold);
				vector<Ray>::iterator ray_it = currentRays.begin();
				while(ray_it!= currentRays.end()){
					if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
					else ray_it = currentRays.erase(ray_it);
				}
				eventReconstructed+=currentRays.size();
				eventSuitable+=current_full_event->get_clus_N()/(CM_N+MG_N);
				processed++;
				for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
					flux_map_h->Fill(it->eval_X(z),it->eval_Y(z));
				}
				if((time(NULL)-last_time) > 5){
					current_full_event->EventDisplay(cDisplay);
					cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << processed << flush;
					cMap->cd();
					flux_map_h->Draw("colz");
					cMap->Modified();
					cMap->Update();
					stat_text->Clear();
					ostringstream text;
					text << "processed events : " << processed;
					stat_text->AddText(text.str().c_str());
					/*
					text.str("");
					text << "potential tracks : " << eventSuitable;
					stat_text->AddText(text.str().c_str());
					*/
					text.str("");
					text << "actual tracks : " << eventReconstructed;
					stat_text->AddText(text.str().c_str());
					cStats->cd();
					stat_text->Draw();
					cStats->Modified();
					cStats->Update();
					last_time = time(NULL);
				}
				delete current_full_event;
				for(vector<Event*>::iterator it = events.begin();it!=events.end();++it){
					delete *it;
				}
			}
		}
		cout << endl;
		if(max_event>0 && processed>max_event) break;
	}
	delete current_data_reader;
}
