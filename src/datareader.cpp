#define datareader_cpp
#include "datareader.h"
#include "tomography.h"
#include "detector.h"
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>

#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TMath.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

using std::cout;
using std::endl;
using std::flush;
using std::map;
using std::vector;
using std::ostringstream;
using std::ofstream;

using TMath::Min;

using boost::property_tree::ptree;

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
}
void DataReader::Init(map<int,Tomography::det_type> det_type_by_asic_, map<int,int> det_n_by_asic_, string PedName_, string RMSName_, string outFileName, long max_event_){
	det_type_by_asic = det_type_by_asic_;
	det_n_by_asic = det_n_by_asic_;
	int CM_N = 0;
	int MG_N = 0;
	for(map<int,Tomography::det_type>::iterator type_it = det_type_by_asic.begin();type_it!=det_type_by_asic.end();++type_it){
		if(type_it->second == Tomography::MG) MG_N++;
		else if(type_it->second == Tomography::CM) CM_N++;
	}
	if(CM_N>0) cout << "Warning CM are not totally supported !" << endl;
	if(outFileName.size()>0){
		outTree = new Tsignal_W(outFileName,CM_N,MG_N);
	}
	else outTree = NULL;
	reader = NULL;
	DAQtype = Tomography::unknown_elec;
	max_event = max_event_;
	mapping = NULL;
	Nevent = -1;
	evttime = 0;
	if(MG_N>0){
		StripAmpl[Tomography::MG] = vector<vector<vector<float> > >(MG_N,vector<vector<float> >(MG_Detector::Nstrip,vector<float>(Tomography::Nsample,0)));
		Ped[Tomography::MG] = vector<vector<float> >(MG_N,vector<float>(MG_Detector::Nstrip,0));
	}
	if(CM_N>0){
		StripAmpl[Tomography::CM] = vector<vector<vector<float> > >(CM_N,vector<vector<float> >(CM_Detector::Nstrip,vector<float>(Tomography::Nsample,0)));
		Ped[Tomography::CM] = vector<vector<float> >(CM_N,vector<float>(CM_Detector::Nstrip,0));
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
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	int CM_N = 0;
	int MG_N = 0;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::CM;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("cm_n");
		CM_N++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		det_type_by_asic[child.second.get<int>("asic_n")] = Tomography::MG;
		det_n_by_asic[child.second.get<int>("asic_n")] = child.second.get<int>("mg_n");
		MG_N++;
	}
	if((total_CM_N!=CM_N) || (total_MG_N!=MG_N)){
		cout << "problem in detectors number" << endl;
		return;
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
		outTree = new Tsignal_W(signalName,CM_N,MG_N);
	}
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
	if(MG_N>0){
		StripAmpl[Tomography::MG] = vector<vector<vector<float> > >(MG_N,vector<vector<float> >(MG_Detector::Nstrip,vector<float>(Tomography::Nsample,0)));
		Ped[Tomography::MG] = vector<vector<float> >(MG_N,vector<float>(MG_Detector::Nstrip,0));
	}
	if(CM_N>0){
		StripAmpl[Tomography::CM] = vector<vector<vector<float> > >(CM_N,vector<vector<float> >(CM_Detector::Nstrip,vector<float>(Tomography::Nsample,0)));
		Ped[Tomography::CM] = vector<vector<float> >(CM_N,vector<float>(CM_Detector::Nstrip,0));
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
	while((!(reader->is_end())) && !((event_nb>max_event)*(max_event>0))){
		process_event();
		event_nb++;
	}
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
			for(int k=0;k<Tomography::Nsample;k++){
				StripAmpl[det_type_by_asic[map_it->first]][map_it->second][channel][k] = reader->get_data(map_it->first,j,k);
			}
		}
	}
	if(outTree != NULL){
		outTree->fillTree_raw(Nevent,evttime,StripAmpl[Tomography::MG],StripAmpl[Tomography::CM]);
	}
}
void DataReader::compute_ped(){
	int total_CM_N = StripAmpl[Tomography::CM].size();
	int total_MG_N = StripAmpl[Tomography::MG].size();
	for(int i=0;i<total_MG_N;i++){
		for(int j=0;j<MG_Detector::Nstrip;j++){
			Ped[Tomography::MG][i][j] = 0;
		}
	}
	for(int i=0;i<total_CM_N;i++){
		for(int j=0;j<CM_Detector::Nstrip;j++){
			Ped[Tomography::CM][i][j] = 0;
		}
	}
	long nentries = outTree->T->GetEntriesFast();
	for(long n=0;n<nentries;n++){
		StripAmpl = outTree->read_raw(n);
		for(int i=0;i<total_MG_N;i++){
			for(int j=0;j<MG_Detector::Nstrip;j++){
				vector<float> current_strip(Tomography::Nsample,0);
				for(int k=0;k<Tomography::Nsample;k++){
					current_strip[k] = StripAmpl[Tomography::MG][i][j][k];
				}
				sort(current_strip.begin(),current_strip.end());
				Ped[Tomography::MG][i][j] += current_strip[Tomography::Nsample/2];
			}
		}
		for(int i=0;i<total_CM_N;i++){
			for(int j=0;j<CM_Detector::Nstrip;j++){
				vector<float> current_strip(Tomography::Nsample,0);
				for(int k=0;k<Tomography::Nsample;k++){
					current_strip[k] = StripAmpl[Tomography::CM][i][j][k];
				}
				sort(current_strip.begin(),current_strip.end());
				Ped[Tomography::CM][i][j] += current_strip[Tomography::Nsample/2];
			}
		}
		if((n%100) == 0) cout << "\r" << "computing pedestal (" << n << "/" << nentries << ")" << flush;
	}
	cout << "\r" << "computing pedestal (" << nentries << "/" << nentries << ")" << endl;
	cout << "writing it to file..." << flush;
	ofstream pedFile(PedName.c_str());
	for(int i=0;i<total_CM_N;i++){
		for(int j=0;j<CM_Detector::Nstrip;j++){
			Ped[Tomography::CM][i][j] /= nentries;
			pedFile << i << " " << j << " " << Ped[Tomography::CM][i][j] << "\n";
		}
	}
	for(int i=0;i<total_MG_N;i++){
		for(int j=0;j<MG_Detector::Nstrip;j++){
			Ped[Tomography::MG][i][j] /= nentries;
			pedFile << i << " " << j << " " << Ped[Tomography::MG][i][j] << "\n";
		}
	}
	pedFile.close();
	cout << "!" << endl;
}
void DataReader::read_ped(){
	ifstream in;
	in.open(PedName.c_str());
	int ped_strip, det;
	Ped.clear();
	int total_CM_N = StripAmpl[Tomography::CM].size();
	int total_MG_N = StripAmpl[Tomography::MG].size();
	Ped[Tomography::CM] = vector<vector<float> >(total_CM_N,vector<float>(CM_Detector::Nstrip,0));
	Ped[Tomography::MG] = vector<vector<float> >(total_MG_N,vector<float>(MG_Detector::Nstrip,0));
	int n_lines = 0;
	while(in.good() && n_lines<((total_CM_N*CM_Detector::Nstrip)+(total_MG_N*MG_Detector::Nstrip))){
		float current_ped;
		in >> det >> ped_strip >> current_ped;
		Tomography::det_type current_type = (n_lines<(total_CM_N*CM_Detector::Nstrip)) ? Tomography::CM : Tomography::MG;
		if(det<0 || det>(total_MG_N+total_CM_N-1)){
			cout << "problem reading Ped file" << endl;
		}
		else if(det<total_CM_N && ped_strip>63){
			cout << "problem reading Ped file" << endl;
		}
		else if(ped_strip>60){
			cout << "problem reading Ped file" << endl;
		}
		Ped[current_type][det][ped_strip] = current_ped;
		n_lines++;
	}
	in.close();
}
void DataReader::compute_RMSPed(){
	double Ymin=-500;
	double Ymax=500;
	int bin_n = 500;
	int sample_min = 1;
	int sample_max = Tomography::Nsample-1;//Min(Nsample,4);
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
	for(long n=0;n<nentries;n++){
		StripAmpl = outTree->read_corr(n);
		for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist.begin();type_it!=ampl_hist.end();++type_it){
			for(unsigned int i=0;i<(type_it->second).size();i++){
				for(unsigned int j=0;j<(type_it->second)[i].size();j++){
					for(int k=sample_min;k<sample_max;k++){
						(type_it->second)[i][j]->Fill(StripAmpl[type_it->first][i][j][k]);
					}
				}
			}
		}
	}
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
	while((!((event_nb>max_event)*(max_event>0))) && event_nb<outTree->T->GetEntriesFast()){
		StripAmpl = outTree->read_raw(event_nb);
		do_ped_sub_event();
		outTree->fillTree_ped(StripAmpl[Tomography::MG],StripAmpl[Tomography::CM]);
		event_nb++;
	}
	if(outTree != NULL){
		outTree->Write();
	}
}
void DataReader::do_common_noise_sub(){
	long event_nb = 0;
	if(outTree != NULL){
		outTree->Reset_corr();
	}
	while((!((event_nb>max_event)*(max_event>0))) && event_nb<outTree->T->GetEntriesFast()){
		StripAmpl = outTree->read_ped(event_nb);
		do_common_noise_sub_event();
		outTree->fillTree_corr(StripAmpl[Tomography::MG],StripAmpl[Tomography::CM]);
		event_nb++;
	}
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
	map<Tomography::det_type,int> Nstrip;
	Nstrip[Tomography::CM] = CM_Detector::Nstrip;
	Nstrip[Tomography::MG] = MG_Detector::Nstrip;
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator it = StripAmpl.begin();it!=StripAmpl.end();++it){
		for(vector<vector<vector<float> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::Nsample;k++){
				for(int det_div=0;det_div<(Tomography::CMN_div.find(it->first))->second;det_div++){
					int strip_nb = Nstrip[it->first]/(Tomography::CMN_div.find(it->first))->second + Nstrip[it->first]%(Tomography::CMN_div.find(it->first))->second;
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<Nstrip[it->first]);j++){
						current_sample[j] = (*jt)[j+strip_offset][k];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<Nstrip[it->first]);j++){
						(*jt)[j+strip_offset][k] -= median;
					}
				}
			}
		}
	}
}
void DataReader::do_ped_CMN_sub_event(){
	map<Tomography::det_type,int> Nstrip;
	Nstrip[Tomography::CM] = CM_Detector::Nstrip;
	Nstrip[Tomography::MG] = MG_Detector::Nstrip;
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator it = StripAmpl.begin();it!=StripAmpl.end();++it){
		vector<vector<float> >::iterator ped_jt = Ped[it->first].begin();
		for(vector<vector<vector<float> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::Nsample;k++){
				for(int det_div=0;det_div<(Tomography::CMN_div.find(it->first))->second;det_div++){
					int strip_nb = Nstrip[it->first]/((Tomography::CMN_div.find(it->first))->second) + Nstrip[it->first]%((Tomography::CMN_div.find(it->first))->second);
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
	}
}
long DataReader::get_event_n(){
	return Nevent;
}
double DataReader::get_evttime(){
	return evttime;
}
map<Tomography::det_type,vector<vector<vector<float> > > > DataReader::get_data(){
	return StripAmpl;
}
int DataReader::Dream_mapping(Tomography::det_type det,int channel){
	if(det == Tomography::MG){
		return (channel + 1 - (2*(channel%2)));
	}
	return channel;
}
int DataReader::Feminos_mapping(Tomography::det_type det,int channel){
	if(det == Tomography::MG){
		int tmpchan = channel - 2 - (channel>13) - (channel>24) - (channel>47) - (channel>58);
		if(tmpchan>15 && tmpchan<48) return tmpchan;
		else return (tmpchan + 1 - (2*(tmpchan%2)));
	}
	return channel;
}