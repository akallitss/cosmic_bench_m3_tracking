#define datareader_cpp
#include "datareader.h"
#include "detector.h"
#include "ElecReader.h"
#include "Tsignal_W.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TMath.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

using std::cout;
using std::endl;
using std::flush;
using std::ostringstream;
using std::ofstream;

using TMath::Min;

DataReader::DataReader(){
	reader = NULL;
	outTree = NULL;
	DAQtype = Tomography::unknown_elec;
	det_type_by_asic.clear();
	det_n_by_asic.clear();
	max_event = -1;
	mapping = NULL;
	Nevent = -1;
	evttime = 0;
	StripAmpl.clear();
	PedName.clear();
	RMSName.clear();
	Ped.clear();
	det_N.clear();
}
void DataReader::Init(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string PedName_, string RMSName_, string outFileName, long max_event_){
	det_type_by_asic = det_type_by_asic_;
	det_n_by_asic = det_n_by_asic_;
	for(map<int,Tomography::det_type>::iterator type_it = det_type_by_asic.begin();type_it!=det_type_by_asic.end();++type_it){
		det_N[type_it->second]++;
	}
	if(outFileName.size()>0){
		outTree = new Tsignal_W(outFileName,det_N);
	}
	else outTree = NULL;
	reader = NULL;
	DAQtype = Tomography::unknown_elec;
	max_event = max_event_;
	mapping = NULL;
	Nevent = -1;
	evttime = 0;
	for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
		if(type_it->second > 0){
			StripAmpl[type_it->first] = vector<vector<vector<float> > >(type_it->second,vector<vector<float> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<float>(Tomography::get_instance()->get_Nsample(),0)));
			Ped[type_it->first] = vector<vector<float> >(type_it->second,vector<float>(Tomography::Static_Detector[type_it->first]->get_Nchannel(),0));
		}
	}
	PedName = PedName_;
	RMSName = RMSName_;
}
DataReader::DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, map<int,int> feu_id_to_n_, int first_index_, int last_index_, string PedName_, string RMSName_){
	Init(det_type_by_asic_,det_n_by_asic_,PedName_,RMSName_);
	DAQtype = Tomography::Dream;
	reader = new DreamElecReader(base_name_,feu_id_to_n_,first_index_,last_index_);
	mapping = &Dream_mapping;
}
DataReader::DataReader(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string base_name_, vector<int> fem_id, int first_index_, int last_index_, string PedName_, string RMSName_){
	Init(det_type_by_asic_,det_n_by_asic_,PedName_,RMSName_);
	DAQtype = Tomography::Feminos;
	reader = new FeminosElecReader(base_name_,fem_id,first_index_,last_index_);
	mapping = &Feminos_mapping;
}
DataReader::DataReader(ptree config_tree, bool save_to_disk){
	CosmicBench current_CB(config_tree);
	det_N.clear();
	for(int i=0;i<current_CB.get_det_N_tot();i++){
		Detector * current_det = current_CB.get_detector(i);
		det_type_by_asic[current_det->get_asic_n()] = current_det->get_type();
		det_n_by_asic[current_det->get_asic_n()] = current_det->get_n_in_tree();
		det_N[current_det->get_type()]++;
	}
	DAQtype = Tomography::str_to_elec(config_tree.get<string>("electronic_type"));
	string data_file_basename = config_tree.get<string>("data_file_basename");
	string signalName = config_tree.get<string>("signal_file");
	PedName = config_tree.get<string>("Ped");
	RMSName = config_tree.get<string>("RMSPed");
	max_event = config_tree.get<long>("max_event");
	int first_index = config_tree.get<int>("data_file_first");
	int last_index = config_tree.get<int>("data_file_last");
	if(save_to_disk){
		outTree = new Tsignal_W(signalName,det_N);
	}
	else outTree = NULL;
	if(DAQtype==Tomography::Dream){
		map<int,int> feu_id_to_n;
		BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("FEU")){
			feu_id_to_n[child.second.get<int>("id")] = child.second.get<int>("n");
		}
		reader = new DreamElecReader(data_file_basename,feu_id_to_n,first_index,last_index);
		mapping = &Dream_mapping;
	}
	else if(DAQtype==Tomography::Feminos){
		vector<int> fem_id;
		for(map<int,int>::iterator map_it=det_n_by_asic.begin();map_it!=det_n_by_asic.end();++map_it){
			bool exists = false;
			for(vector<int>::iterator vec_it=fem_id.begin();vec_it!=fem_id.end();++vec_it){
				if((*vec_it) == ((map_it->first)/Tomography::Nasic_Feminos)) exists = true;
			}
			if(!exists) fem_id.push_back((map_it->first)/Tomography::Nasic_Feminos);
		}
		reader = new FeminosElecReader(data_file_basename,fem_id,first_index,last_index);
		mapping = &Feminos_mapping;
	}
	else{
		cout << "unknown electronic type !" << endl;
		reader = NULL;
		mapping = NULL;
	}
	Nevent = -1;
	evttime = 0;
	for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
		if(type_it->second > 0){
			StripAmpl[type_it->first] = vector<vector<vector<float> > >(type_it->second,vector<vector<float> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<float>(Tomography::get_instance()->get_Nsample(),0)));
			Ped[type_it->first] = vector<vector<float> >(type_it->second,vector<float>(Tomography::Static_Detector[type_it->first]->get_Nchannel(),0));
		}
	}
}
DataReader::~DataReader(){
	if(reader != NULL) delete reader;
	if(outTree != NULL){
		outTree->Write();
		outTree->CloseFile();
		delete outTree;
	}
	DAQtype = Tomography::unknown_elec;
	det_type_by_asic.clear();
	det_n_by_asic.clear();
	StripAmpl.clear();
	Nevent = -1;
	evttime = 0;
	PedName.clear();
	RMSName.clear();
	Ped.clear();
	det_N.clear();
}
void DataReader::process(){
	if(reader == NULL){
		cout << "Data Reader not initialized !" << endl;
		return;
	}
	long event_nb = 0;
	if(outTree != NULL){
		outTree->Reset_raw();
	}
	while((!(reader->is_end())) && !((event_nb>max_event)*(max_event>0)) && Tomography::get_instance()->get_can_continue()){
		if((event_nb%100) == 0) cout << "\r" << "event processed : " << event_nb << flush;
		process_event();
		event_nb++;
	}
	cout << "\r" << "event processed : " << event_nb << endl;
	if(outTree != NULL){
		outTree->Write();
	}
}
void DataReader::process_event(){
	if(reader == NULL){
		cout << "Data Reader not initialized !" << endl;
		return;
	}
	reader->read_next_event();
	Nevent = reader->get_event_n();
	evttime = reader->get_evttime();
	for(map<int,int>::iterator map_it=det_n_by_asic.begin();map_it!=det_n_by_asic.end();++map_it){
		for(int j=0;j<Tomography::Nchannel;j++){
			int channel = mapping(det_type_by_asic[map_it->first],j);
			if(channel<0) continue;
			if(static_cast<unsigned int>(channel)>=StripAmpl[det_type_by_asic[map_it->first]][map_it->second].size()) continue;
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				StripAmpl[det_type_by_asic[map_it->first]][map_it->second][channel][k] = reader->get_data(map_it->first,j,k);
			}
		}
	}
	if(outTree != NULL){
		outTree->fillTree_raw(Nevent,evttime,StripAmpl);
	}
}
void DataReader::compute_ped(){
	for(map<Tomography::det_type,unsigned short>::const_iterator type_it = det_N.begin();type_it!=det_N.end();++type_it){
		for(unsigned int i=0;i<type_it->second;i++){
			for(int j=0;j<(Tomography::Static_Detector[type_it->first]->get_Nchannel());j++){
				Ped[type_it->first][i][j] = 0;
			}
		}
	}
	long nentries = outTree->T->GetEntriesFast();
	outTree->disable_data_branches();
	outTree->enable_raw_branches();
	for(long n=0;n<nentries && Tomography::get_instance()->get_can_continue();n++){
		StripAmpl = outTree->read_raw<float>(n);
		for(map<Tomography::det_type,vector<vector<float> > >::iterator type_it=Ped.begin();type_it!=Ped.end();++type_it){
			for(unsigned int i=0;i<(type_it->second).size();i++){
				for(unsigned int j=0;j<(type_it->second)[i].size();j++){
					vector<float> current_strip(Tomography::get_instance()->get_Nsample(),0);
					for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
						current_strip[k] = StripAmpl[type_it->first][i][j][k];
					}
					sort(current_strip.begin(),current_strip.end());
					(type_it->second)[i][j] += current_strip[Tomography::get_instance()->get_Nsample()/2];
				}
			}
		}
		if((n%100) == 0) cout << "\r" << "computing pedestal (" << n << "/" << nentries << ")" << flush;
	}
	outTree->enable_all_branches();
	cout << "\r" << "computing pedestal (" << nentries << "/" << nentries << ")" << endl;
	cout << "writing it to file : " << PedName << "..." << flush;
	ofstream pedFile(PedName.c_str());
	for(map<Tomography::det_type,vector<vector<float> > >::iterator type_it=Ped.begin();type_it!=Ped.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				(type_it->second)[i][j] /= nentries;
				pedFile << i << " " << j << " " << (type_it->second)[i][j] << "\n";
			}
		}
	}
	pedFile.close();
	cout << "!" << endl;
}
void DataReader::read_ped(){
	cout << "Reading pedestal from : " << PedName << "..." << endl;
	map<Tomography::det_type,vector<vector<double> > > tmp_Ped = CosmicBench::read_pedfile(PedName,det_N);
	for(map<Tomography::det_type,vector<vector<float> > >::iterator type_it=Ped.begin();type_it!=Ped.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				(type_it->second)[i][j] = tmp_Ped[type_it->first][i][j];
			}
		}
	}
	cout << "done !" << endl;
}
void DataReader::compute_RMSPed(){
	double Ymin=-500;
	double Ymax=500;
	int bin_n = 500;
	int sample_min = 1;
	int sample_max = Tomography::get_instance()->get_Nsample()-1;//Min(Nsample,4);
	Long64_t tot_event = 1000;
	long nentries = Min(outTree->T->GetEntries(),tot_event);
	map<Tomography::det_type,vector<vector<TH1F*> > > ampl_hist;
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator type_it = StripAmpl.begin();type_it!=StripAmpl.end();++type_it){
		//ampl_hist[type_it->first] = vector<vector<TH1F*> >((type_it->second).size(),vector<TH1F*>((type_it->second)[0].size(),NULL));
		ampl_hist[type_it->first].resize((type_it->second).size());
		for(unsigned int j=0;j<ampl_hist[type_it->first].size();j++){
			ampl_hist[type_it->first][j].resize((type_it->second)[0].size());
			for(unsigned int i=0;i<ampl_hist[type_it->first][j].size();i++){
				ostringstream name;
				name << "ampl_hist_" << type_it->first << j << "_" << i;
				ampl_hist[type_it->first][j][i] = new TH1F(name.str().c_str(),name.str().c_str(),bin_n,Ymin,Ymax);
			}
		}
	}
	outTree->disable_data_branches();
	outTree->enable_corr_branches();
	for(long n=0;n<nentries && Tomography::get_instance()->get_can_continue();n++){
		StripAmpl = outTree->read_corr<float>(n);
		for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist.begin();type_it!=ampl_hist.end();++type_it){
			for(unsigned int i=0;i<(type_it->second).size();i++){
				for(unsigned int j=0;j<(type_it->second)[i].size();j++){
					for(int k=sample_min;k<sample_max;k++){
						(type_it->second)[i][j]->Fill(StripAmpl[type_it->first][i][j][k]);
					}
				}
			}
		}
		if((n%100) == 0) cout << "\rcomputing RMS (" << n << "/" << nentries << ")" << flush;
	}
	outTree->enable_all_branches();
	cout << "\rcomputing RMS (" << nentries << "/" << nentries << ")" << endl;
	ofstream RMSPedFile(RMSName.c_str());
	for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist.begin();type_it!=ampl_hist.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				TFitResultPtr res = (type_it->second)[i][j]->Fit("gaus","SQN");
				RMSPedFile << i << " " << j << " " << res->Parameter(2) << "\n";
				delete (type_it->second)[i][j];
			}
		}
	}
	RMSPedFile.close();
}
void DataReader::do_ped_sub(){
	read_ped();
	long event_nb = 0;
	if(outTree != NULL){
		outTree->Reset_ped();
	}
	long total_event = outTree->T->GetEntriesFast();
	if(max_event>0 && max_event<total_event) total_event = max_event;
	outTree->disable_data_branches();
	outTree->enable_raw_branches();
	outTree->enable_ped_branches();
	while(event_nb<total_event && Tomography::get_instance()->get_can_continue()){
		if((event_nb%100) == 0) cout << "\rsubstracting pedestal (" << event_nb << "/" << total_event << ")" << flush;
		StripAmpl = outTree->read_raw<float>(event_nb);
		do_ped_sub_event();
		outTree->fillTree_ped(StripAmpl);
		event_nb++;
	}
	outTree->enable_all_branches();
	cout << "\rsubstracting pedestal (" << total_event << "/" << total_event << ")" << endl;
	if(outTree != NULL){
		outTree->Write();
	}
}
void DataReader::do_common_noise_sub(){
	long event_nb = 0;
	if(outTree != NULL){
		outTree->Reset_corr();
	}
	long total_event = outTree->T->GetEntriesFast();
	if(max_event>0 && max_event<total_event) total_event = max_event;
	outTree->disable_data_branches();
	outTree->enable_ped_branches();
	outTree->enable_corr_branches();
	while(event_nb<total_event && Tomography::get_instance()->get_can_continue()){
		if((event_nb%100) == 0) cout << "\rsubstracting common noise (" << event_nb << "/" << total_event << ")" << flush;
		StripAmpl = outTree->read_ped<float>(event_nb);
		do_common_noise_sub_event();
		outTree->fillTree_corr(StripAmpl);
		event_nb++;
	}
	outTree->enable_all_branches();
	cout << "\rsubstracting common noise (" << total_event << "/" << total_event << ")" << endl;
	if(outTree != NULL){
		outTree->Write();
	}
}
void DataReader::do_ped_sub_event(){
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator type_it=StripAmpl.begin();type_it!=StripAmpl.end();++type_it){
		vector<vector<float> >::iterator ped_det_it = Ped[type_it->first].begin();
		for(vector<vector<vector<float> > >::iterator det_it=(type_it->second).begin();det_it!=(type_it->second).end();++det_it){
			vector<float>::iterator ped_strip_it = ped_det_it->begin();
			for(vector<vector<float> >::iterator strip_it=det_it->begin();strip_it!=det_it->end();++strip_it){
				for(vector<float>::iterator sample_it=strip_it->begin();sample_it!=strip_it->end();++sample_it){
					*sample_it -= *ped_strip_it;
					//*sample_it -= Ped[type_it->first][det_it - (type_it->second).begin()][strip_it - det_it->begin()];
				}
				++ped_strip_it;
			}
			++ped_det_it;
		}
	}
}
void DataReader::do_common_noise_sub_event(){
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator it = StripAmpl.begin();it!=StripAmpl.end();++it){
		for(vector<vector<vector<float> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				for(int det_div=0;det_div<(Tomography::Static_Detector[it->first]->get_CMN_div());det_div++){
					int strip_nb = (Tomography::Static_Detector[it->first]->get_Nchannel())/(Tomography::Static_Detector[it->first]->get_CMN_div()) + (Tomography::Static_Detector[it->first]->get_Nchannel())%(Tomography::Static_Detector[it->first]->get_CMN_div());
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						current_sample[j] = (*jt)[j+strip_offset][k];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						(*jt)[j+strip_offset][k] -= median;
					}
				}
			}
		}
	}
}
void DataReader::do_ped_CMN_sub_event(){
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator it = StripAmpl.begin();it!=StripAmpl.end();++it){
		vector<vector<float> >::iterator ped_jt = Ped[it->first].begin();
		for(vector<vector<vector<float> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				for(int det_div=0;det_div<(Tomography::Static_Detector[it->first]->get_CMN_div());det_div++){
					int strip_nb = (Tomography::Static_Detector[it->first]->get_Nchannel())/(Tomography::Static_Detector[it->first]->get_CMN_div()) + (Tomography::Static_Detector[it->first]->get_Nchannel())%(Tomography::Static_Detector[it->first]->get_CMN_div());
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						current_sample[j] = (*jt)[j+strip_offset][k] - (*ped_jt)[j+strip_offset];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						(*jt)[j+strip_offset][k] -= median + (*ped_jt)[j+strip_offset];
					}
				}
			}
			++ped_jt;
		}
	}
}
template<typename S, typename T>
map<Tomography::det_type,vector<vector<vector<S> > > > DataReader::do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<T> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in){
	map<Tomography::det_type,vector<vector<vector<S> > > > data_out;
	for(typename map<Tomography::det_type,vector<vector<vector<T> > > >::iterator it = data_in.begin();it!=data_in.end();++it){
		data_out[it->first] = vector<vector<vector<S> > >((it->second).size(),vector<vector<S> >(Tomography::Static_Detector[it->first]->get_Nchannel(),vector<S>(Tomography::get_instance()->get_Nsample(),0)));
		vector<vector<float> >::iterator ped_jt = ped_in[it->first].begin();
		typename vector<vector<vector<S> > >::iterator out_it = data_out[it->first].begin();
		for(typename vector<vector<vector<T> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				for(int det_div=0;det_div<(Tomography::Static_Detector[it->first]->get_CMN_div());det_div++){
					int strip_nb = (Tomography::Static_Detector[it->first]->get_Nchannel())/(Tomography::Static_Detector[it->first]->get_CMN_div()) + (Tomography::Static_Detector[it->first]->get_Nchannel())%(Tomography::Static_Detector[it->first]->get_CMN_div());
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						current_sample[j] = (*jt)[j+strip_offset][k] - (*ped_jt)[j+strip_offset];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						(*out_it)[j+strip_offset][k] = (*jt)[j+strip_offset][k] - median - (*ped_jt)[j+strip_offset];
					}
				}
			}
			++ped_jt;
			++out_it;
		}
	}
	return data_out;
}
template map<Tomography::det_type,vector<vector<vector<float> > > > DataReader::do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<float> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
template map<Tomography::det_type,vector<vector<vector<double> > > > DataReader::do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<float> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
template map<Tomography::det_type,vector<vector<vector<float> > > > DataReader::do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<double> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
template map<Tomography::det_type,vector<vector<vector<double> > > > DataReader::do_ped_CMN_sub_event(map<Tomography::det_type,vector<vector<vector<double> > > > data_in, map<Tomography::det_type,vector<vector<float> > > ped_in);
long DataReader::get_event_n() const{
	return Nevent;
}
long DataReader::get_max_event() const{
	return max_event;
}
double DataReader::get_evttime() const{
	return evttime;
}
map<Tomography::det_type,vector<vector<vector<float> > > > DataReader::get_data() const{
	return StripAmpl;
}
map<Tomography::det_type,vector<vector<float> > > DataReader::get_Ped() const{
	return Ped;
}
bool DataReader::is_end() const{
	return (reader->is_end());
}
int DataReader::Dream_mapping(Tomography::det_type det,int channel){
	return Tomography::Static_Detector[det]->dream_mapping(channel);
}
int DataReader::Feminos_mapping(Tomography::det_type det,int channel){
	return Tomography::Static_Detector[det]->feminos_mapping(channel);
}