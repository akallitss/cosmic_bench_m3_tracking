#define livedisplay_cpp

#include "liveDisplay.h"
#include "datareader.h"
#include "detector.h"
#include "event.h"

#include <string>
#include <sys/inotify.h>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>

using std::string;
using std::ifstream;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::cout;
using std::endl;
using std::flush;
using std::bitset;

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
	electronic_type = "";
	MG_N = 0;
	CM_N = 0;
	det_type_by_asic.clear();
	det_n_by_asic.clear();
}

liveDisplay::~liveDisplay(){

}

liveDisplay::liveDisplay(string config_file, int max_event_){
	ptree config_tree;
	read_json(config_file, config_tree);
	filenames.clear();
	max_event = max_event_;
	inotify_descriptor = -1;
	file_descriptor = -1;
	inotify_started = false;
	current_file = "";
	MG_N = 0;
	CM_N = 0;
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	ifstream in;
	in.open((config_tree.get<string>("RMSPed")).c_str());
	int rms_strip, det;
	//double RMS[Nstrip_MG*total_MG_N];
	vector<vector<double> > RMS;
	vector<double> empty_vector(61,0);
	for(int i=0;i<(total_MG_N+total_CM_N);i++){
		RMS.push_back(empty_vector);
	}
	int nlines=0;
	while (1) { // read the text file
		int det_n = nlines/61;
		int strip_n = nlines%61;
		if(det_n>(total_MG_N+total_CM_N-1)) break;
		double current_rms;
		in >> det >> rms_strip >> current_rms;
		RMS[det_n][strip_n] = current_rms;
		if (!in.good()) break;
		nlines++;
	}
	in.close();
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		detectors.push_back(new CM_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("is_up"),child.second.get<int>("cm_n"),child.second.get<bool>("use_thin_strip"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle")));
		detectors.back()->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		detectors.back()->set_RMS(RMS[child.second.get<int>("cm_n")]);
		(dynamic_cast<CM_Detector*>(detectors.back()))->set_ClusMaxStripAmplCut_Min_Wide(child.second.get<double>("ClusMaxStripAmplCut_Min_Wide"));
		(dynamic_cast<CM_Detector*>(detectors.back()))->set_ClusSizeCut_Max_Wide(child.second.get<double>("ClusSizeCut_Max_Wide"));
		det_type_by_asic[child.second.get<int>("asic_n")] = "CM";
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("cm_n");
		CM_N++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		detectors.push_back(new MG_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("is_up"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle"),child.second.get<double>("2D_perp_n"),child.second.get<int>("clustering_holes")));
		detectors.back()->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		detectors.back()->set_RMS(RMS[total_CM_N+child.second.get<int>("mg_n")]);
		(dynamic_cast<MG_Detector*>(detectors.back()))->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
		(dynamic_cast<MG_Detector*>(detectors.back()))->set_SRF(child.second.get<double>("SRF.offset"),child.second.get<double>("SRF.gauss"),child.second.get<double>("SRF.lorentz"),child.second.get<double>("SRF.ratio"));
		det_type_by_asic[child.second.get<int>("asic_n")] = "MG";
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("mg_n");
		MG_N++;
	}
	if((total_CM_N!=CM_N) || (total_MG_N!=MG_N)){
		cout << "problem in detectors number" << endl;
		return;
	}
	if(total_CM_N!=0) cout << "warning, CosMultis are not fully supported !" << endl;
	use_srf = config_tree.get<bool>("use_SRF");
	electronic_type = config_tree.get<string>("electronic_type");
	data_file_basename = config_tree.get<string>("data_file_basename");
	ifstream pedFile((config_tree.get<string>("Ped")).c_str());
	/*
	Pedestal_CM = new float[CM_N][DataReader::Nstrip_CM];
	Pedestal_MG = new float[MG_N][DataReader::Nstrip_MG];
	for(int i=0;i<CM_N;i++){
		for(int j=0;j<DataReader::Nstrip_CM;j++){
			pedFile >> i >> j >> Pedestal_CM[i][j];
		}
	}
	for(int i=0;i<MG_N;i++){
		for(int j=0;j<DataReader::Nstrip_MG;j++){
			pedFile >> i >> j >> Pedestal_MG[i][j];
		}
	}
	*/
	Pedestal["CM"] = vector<vector<float> >(CM_N,vector<float>(DataReader::Nstrip_CM,0));
	Pedestal["MG"] = vector<vector<float> >(MG_N,vector<float>(DataReader::Nstrip_MG,0));
	for(int i=0;i<CM_N;i++){
		for(int j=0;j<DataReader::Nstrip_CM;j++){
			pedFile >> i >> j >> Pedestal["CM"][i][j];
		}
	}
	for(int i=0;i<MG_N;i++){
		for(int j=0;j<DataReader::Nstrip_MG;j++){
			pedFile >> i >> j >> Pedestal["MG"][i][j];
		}
	}
	pedFile.close();
	max_file_size = config_tree.get<long>("max_file_size");
}

void liveDisplay::add_file(string filename){
	filenames.push_back(filename);
}
void liveDisplay::add_files(int first,int last){
	string extension = "";
	if(electronic_type == "feminos") extension = ".aqs";
	else if(electronic_type == "dream") extension = ".fdf";
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
	TCanvas * cMap = new TCanvas("flux_map");
	TCanvas * cDisplay = new TCanvas("event_display");
	TCanvas * cStats = new TCanvas("stats");
	int bin_n = 100;
	double margin = 200;
	double x_min = 0;
	double x_max = 500;
	double y_min = 0;
	double y_max = 500;
	map<string,int> detector_div;
	detector_div["CM"] = 2;
	detector_div["MG"] = 2;
	map<string,int> Nstrip;
	Nstrip["CM"] = DataReader::Nstrip_CM;
	Nstrip["MG"] = DataReader::Nstrip_MG;
	int eventReconstructed = 0;
	int eventSuitable = 0;
	int processed = 0;
	clock_t last_time = clock();
	TH2D * flux_map = new TH2D("flux_map","flux_map",bin_n,x_min-margin,x_max+margin,bin_n,y_min-margin,y_max+margin);
	flux_map->SetStats(0);
	TPaveText * stat_text = new TPaveText(0,0,1,1);
	DataReader * current_data_reader = NULL;
	int event_nb = 0;
	if(electronic_type == "feminos"){
		current_data_reader = new FeminosDataReader("live",det_type_by_asic,det_n_by_asic);
		event_nb = dynamic_cast<FeminosDataReader*>(current_data_reader)->get_first_event_nb(filenames.front());
	}
	else if(electronic_type == "dream"){
		current_data_reader = new DreamDataReader("live",det_type_by_asic,det_n_by_asic);
	}
	else return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for(vector<string>::iterator filename_it = filenames.begin();filename_it!=filenames.end();++filename_it){
		cout << *filename_it << endl;
		ifstream data_file(filename_it->c_str(),ifstream::binary);
		start_inotify(*filename_it);
		bool is_open = true;
		int current_pos = data_file.tellg();
		while(is_open){
			current_pos = data_file.tellg();
			data_file.seekg(0,data_file.end);
			bool is_complete = !(data_file.tellg()<max_file_size);
			data_file.seekg(current_pos);
			unsigned int read_mask = 0x00000000;
			if(!(current_pos>=max_file_size || is_complete)) read_mask = read_inotify();
			else read_mask = IN_MODIFY | IN_CLOSE;
			if((read_mask & IN_CLOSE) || is_complete) is_open = false;
			if(!(read_mask & IN_MODIFY)) continue;
			while(data_file.good()){
				map<string,vector<vector<vector<double> > > > event_ampl = current_data_reader->read_event(&data_file,event_nb);
				if(event_ampl.size()==0){
					data_file.seekg(current_pos);
					break;
				}
				event_nb++;
				/*
				//pedestal sub
				for(unsigned int i=0;i<MG_N;i++){
					for(int j=0;j<DataReader::Nstrip_MG;j++){
						for(int k=0;k<DataReader::Nsample;k++){
							event_ampl["MG"][i][j][k] -= Pedestal_MG[i][j];
						}
					}
				}
				for(unsigned int i=0;i<CM_N;i++){
					for(int j=0;j<DataReader::Nstrip_CM;j++){
						for(int k=0;k<DataReader::Nsample;k++){
							event_ampl["CM"][i][j][k] -= Pedestal_CM[i][j];
						}
					}
				}
				//common noise sub
				for(unsigned int i=0;i<MG_N;i++){
					for(int k=0;k<DataReader::Nsample;k++){
						for(int det_div=0;det_div<detector_div_MG;det_div++){
							int strip_nb = DataReader::Nstrip_MG/detector_div_MG;
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
							int strip_nb = DataReader::Nstrip_CM/detector_div_CM;
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
							int strip_nb = DataReader::Nstrip_MG/detector_div_MG;
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
							int strip_nb = DataReader::Nstrip_CM/detector_div_CM;
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
				map<string,vector<vector<float> > >::iterator ped_it = Pedestal.begin();
				for(map<string,vector<vector<vector<double> > > >::iterator it = event_ampl.begin();it!=event_ampl.end();++it){
					vector<vector<float> >::iterator ped_jt = (ped_it->second).begin();
					for(vector<vector<vector<double> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
						for(int k=0;k<DataReader::Nsample;k++){
							for(int det_div=0;det_div<detector_div[it->first];det_div++){
								int strip_nb = Nstrip[it->first]/detector_div[it->first] + Nstrip[it->first]%detector_div[it->first];
								int strip_offset = det_div*strip_nb;
								vector<float> current_sample(strip_nb,0);
								for(int j=0;(j<strip_nb && (j+strip_offset)<Nstrip[it->first]);j++){
									current_sample[j] = (*jt)[j+strip_offset][k] - (*ped_jt)[j+strip_offset];
								}
								sort(current_sample.begin(),current_sample.end());
								float median = current_sample[strip_nb/2];
								for(int j=0;(j<strip_nb && (j+strip_offset)<Nstrip[it->first]);j++){
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
					if((*it)->get_type() == "MG"){
						MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
						events.push_back(new MG_Event(*current_det,event_ampl["MG"][current_det->get_mg_n_in_tree()],use_srf,event_nb));
						(events.back())->MultiCluster();
						(events.back())->do_cuts();
					}
					else if((*it)->get_type() == "CM"){
						CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
						events.push_back(new CM_Event(*current_det,event_ampl["CM"][current_det->get_cm_n_in_tree()],use_srf,event_nb));
						(events.back())->MultiCluster();
						(events.back())->do_cuts();
					}
				}
				CosmicBenchEvent * current_full_event = new CosmicBenchEvent(this,events);
				vector<Ray> currentRays = current_full_event->get_absorption_rays();
				eventReconstructed+=currentRays.size();
				eventSuitable+=current_full_event->get_clus_N()/(CM_N+MG_N);
				processed++;
				for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
					if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
						flux_map->Fill(it->eval_X(z),it->eval_Y(z));
					}
				}
				if((static_cast<float>(clock()-last_time)/CLOCKS_PER_SEC) >1.){
					current_full_event->EventDisplay(cDisplay);
					cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << processed << flush;
					cMap->cd();
					flux_map->Draw("colz");
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
					last_time = clock();
				}
				delete current_full_event;
				for(vector<Event*>::iterator it = events.begin();it!=events.end();++it){
					delete *it;
				}
			}
		}
		cout << endl;
	}
	delete current_data_reader;
}