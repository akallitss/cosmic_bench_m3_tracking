#define signal_cpp
#include "signal.h"
#include "event.h"
#include "Tanalyse.h"
#include "ray.h"
#include "datareader.h"

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <sstream>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLine.h>

//Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

using std::ifstream;
using std::string;
using std::cout;
using std::endl;
using std::flush;
using std::map;
using std::pair;
using std::ostringstream;
using std::setw;
using std::setfill;

//boost
using boost::property_tree::ptree;

Signal::Signal(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	electronic_type = config_tree.get<string>("electronic_type");
	data_file_basename = config_tree.get<string>("data_file_basename");
	signalName = config_tree.get<string>("signal_file");
	PedName = config_tree.get<string>("Ped");
	RMSName = config_tree.get<string>("RMSPed");
	max_event = config_tree.get<int>("max_event");
	data_file_first = config_tree.get<int>("data_file_first");
	data_file_last = config_tree.get<int>("data_file_last");
	cout << signalName << endl;
	cout << RMSName << endl;
	CM_n = 0;
	MG_n = 0;
	ifstream fIn_test(signalName.c_str());
	bool exists = fIn_test.good();
	fIn_test.close();
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	ifstream in;
	in.open(RMSName.c_str());
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
		CM_n++;
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
		MG_n++;
	}
	if((total_CM_N!=CM_n) || (total_MG_N!=MG_n)){
		cout << "problem in detectors number" << endl;
		return;
	}
	CM_N = CM_n;
	MG_N = MG_n;
	if(total_CM_N!=0) cout << "warning, CosMultis are not fully supported !" << endl;
	if(exists){
		TFile *fIn = new TFile(signalName.c_str(),"READ");
		TTree * treeIn = (TTree*)(fIn->Get("T"));
		Init(treeIn,CM_n,MG_n);
	}
	else{
		Init(0,CM_n,MG_n);
		cout << "Waring, signal file is missing !" << endl;
	}
	analyseTree = config_tree.get<string>("Tree");
	use_srf = config_tree.get<bool>("use_SRF");
}
Signal::~Signal(){

}
void Signal::MultiCluster(){
	cout << "destination file : " << analyseTree << endl;
	Tanalyse * analyseFile = new Tanalyse(analyseTree,CM_n,MG_n);
	long nentries = fChain->GetEntriesFast();
	for(int i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);
		vector<MG_Event> mg_events;
		vector<CM_Event> cm_events;

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == "MG"){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				mg_events.push_back(MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent));
				mg_events.back().MultiCluster();
			}
			else if((*it)->get_type() == "CM"){
				CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
				cm_events.push_back(CM_Event(*current_det,get_cm_ampl(current_det->get_cm_n_in_tree()),use_srf,Nevent));
				cm_events.back().MultiCluster();
			}
		}
		
		analyseFile->fillTree(Nevent,0,mg_events,cm_events);
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	analyseFile->Write();
	analyseFile->CloseFile();
	//delete analyseFile;
}

void Signal::ElecToAnalyse(){
	cout << "Reading Pedestal : " << PedName << endl;
	map<string,vector<vector<float> > > Pedestal;
	Pedestal["CM"] = vector<vector<float> >(CM_N,vector<float>(DataReader::Nstrip_CM,0));
	Pedestal["MG"] = vector<vector<float> >(MG_N,vector<float>(DataReader::Nstrip_MG,0));
	ifstream pedFile(PedName.c_str());
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
	cout << "destination file : " << analyseTree << endl;
	Tanalyse * analyseFile = new Tanalyse(analyseTree,CM_n,MG_n);
	Nevent = 0;
	int event_nb = 0;
	string extension = "";
	DataReader * current_data_reader;
	if(electronic_type == "feminos"){
		extension = ".aqs";
		current_data_reader = new FeminosDataReader("tmp",det_type_by_asic,det_n_by_asic);
		ostringstream name;
		name << data_file_basename << setfill('0') << setw(3) << data_file_first << extension;
		Nevent = dynamic_cast<FeminosDataReader*>(current_data_reader)->get_first_event_nb(name.str());
	}
	else if(electronic_type == "dream"){
		extension = "_01.fdf";
		current_data_reader = new DreamDataReader("tmp",det_type_by_asic,det_n_by_asic);
	}
	else return;
	map<string,int> detector_div;
	detector_div["CM"] = 2;
	detector_div["MG"] = 2;
	map<string,int> Nstrip;
	Nstrip["CM"] = DataReader::Nstrip_CM;
	Nstrip["MG"] = DataReader::Nstrip_MG;
	for(int i=data_file_first;i<=data_file_last;i++){
		int event_nb_file = 0;
		ostringstream name;
		name << data_file_basename << setfill('0') << setw(3) << i << extension;
		ifstream data_file(name.str().c_str(),ifstream::binary);
		while(data_file.good()){
			map<string,vector<vector<vector<double> > > > event_ampl = current_data_reader->read_event(&data_file,Nevent);
			if(event_ampl.size()==0) break;
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
			vector<MG_Event> mg_events;
			vector<CM_Event> cm_events;
			for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
				if((*it)->get_type() == "MG"){
					MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
					mg_events.push_back(MG_Event(*current_det,event_ampl["MG"][current_det->get_mg_n_in_tree()],use_srf,event_nb));
					mg_events.back().MultiCluster();
				}
				else if((*it)->get_type() == "CM"){
					CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
					cm_events.push_back(CM_Event(*current_det,event_ampl["CM"][current_det->get_cm_n_in_tree()],use_srf,event_nb));
					cm_events.back().MultiCluster();
				}
			}
			analyseFile->fillTree(Nevent,0,mg_events,cm_events);
			Nevent++;
			event_nb++;
			event_nb_file++;
			if(event_nb_file%100 == 0) cout << "\r" << "processing " << name.str() << " : " << event_nb_file << " (total : " << event_nb << ")" << flush;
		}
		cout << "\r" << "processing " << name.str() << " : " << event_nb_file << " (total : " << event_nb << ")" << endl;
	}
	analyseFile->Write();
	analyseFile->CloseFile();
}

void Signal::HoughTracking(int event_nb){
	if(CM_n!=0){
		cout << "not implemented with CM" << endl;
		return;
	}
	long nentries = fChain->GetEntriesFast();
	if(event_nb<0 || event_nb>=nentries){
		cout << "invalid event number" << endl;
		return;
	}
	LoadTree(event_nb);
	GetEntry(event_nb);
	vector<MG_Cluster> all_cluster;
	vector<Event*> events;
	double max_z = -10000;
	double min_z = 10000;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == "MG"){
			MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
			MG_Event current_event = MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent);
			current_event.HoughCluster();
			vector<MG_Cluster> current_cluster = current_event.get_clusters();
			vector<MG_Cluster>::iterator clus_it = current_cluster.begin();
			while(clus_it!=current_cluster.end()){
				if(!(clus_it->is_suitable(current_det))){
					current_cluster.erase(clus_it);
					clus_it = current_cluster.begin();
				}
				else ++clus_it;
			}
			cout << "number of suitable cluster for detector n°" << current_det->get_mg_n_in_tree() << " : " << current_cluster.size() << endl;
			all_cluster.insert(all_cluster.end(),current_cluster.begin(),current_cluster.end());
			events.push_back(new MG_Event(current_event));
			(events.back())->MultiCluster();
		}
		if((*it)->get_z()>max_z) max_z = (*it)->get_z();
		if((*it)->get_z()<min_z) min_z = (*it)->get_z();
	}
	//max_z+=10;
	//min_z-=10;
	int bin_n = 250;
	double min_coord = -100;
	double max_coord = 600;
	TH2D * hough_space_X = new TH2D("hough_space_X","hough_space_X",bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
	TH2D * hough_space_Y = new TH2D("hough_space_Y","hough_space_Y",bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
	int suitable_clus_n = 0;
	//draw clusters in hough space
	bin_n = 2*bin_n;
	for(vector<MG_Cluster>::iterator it=all_cluster.begin();it!=all_cluster.end();++it){
		suitable_clus_n++;
		if(!it->get_is_up()){
			for(int i=0;i<bin_n;i++){
				double current_coord_up = min_coord +i*(max_coord-min_coord)/bin_n;
				double current_coord_down = current_coord_up + (it->get_pos_mm() - current_coord_up)*(min_z-max_z)/(it->get_z()-max_z);
				if(it->get_is_X()){
					hough_space_X->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
				}
				else{
					hough_space_Y->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
				}
			}
		}
		else{
			for(int i=0;i<bin_n;i++){
				double current_coord_down = min_coord +i*(max_coord-min_coord)/bin_n;
				double current_coord_up = current_coord_down + (it->get_pos_mm() - current_coord_down)*(max_z-min_z)/(it->get_z()-min_z);
				if(it->get_is_X()){
					hough_space_X->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
				}
				else{
					hough_space_Y->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
				}
			}
		}
	}
	cout << "total number of suitable cluster : " << suitable_clus_n << endl;
	//draw "normal" ray in hough space
	CosmicBenchEvent * thisEvent = new CosmicBenchEvent(this,events);
	vector<Ray> rays = thisEvent->get_absorption_rays();
	delete thisEvent;
	TGraph * rays_X = new TGraph();
	TGraph * rays_Y = new TGraph();
	int X_point_nb = 0;
	int Y_point_nb = 0;
	rays_X->SetMarkerStyle(24);
	rays_Y->SetMarkerStyle(24);
	rays_X->SetMarkerSize(2);
	rays_Y->SetMarkerSize(2);
	for(vector<Ray>::iterator it = rays.begin();it!=rays.end();++it){
		rays_X->SetPoint(X_point_nb,it->eval_X(min_z),it->eval_X(max_z));
		rays_Y->SetPoint(Y_point_nb,it->eval_Y(min_z),it->eval_Y(max_z));
		X_point_nb++;
		Y_point_nb++;
	}
	cout << "number of reconstructed rays : " << rays.size() << endl;
	gStyle->SetPalette(55,0);
	TCanvas * cHough = new TCanvas();
	cHough->Divide(2);
	cHough->cd(1);
	hough_space_X->Draw("colz");
	if(rays.size()>0) rays_X->Draw("sameP");
	cHough->cd(2);
	hough_space_Y->Draw("colz");
	if(rays.size()>0) rays_Y->Draw("sameP");
}
map<int,TProfile*> Signal::SignalOverNoise(){
	map<int,TProfile*> global_signal_over_noise;
	if(CM_n!=0){
		cout << "not implemented with CM" << endl;
		return global_signal_over_noise;
	}
	map<int,TCanvas*> cDisplay;
	map<int,TH1D*> global_signal;
	map<int,TH1D*> global_noise;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == "MG"){
			MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
			ostringstream name;
			name << "Multigen_" << current_det->get_mg_n_in_tree();
			//cDisplay[current_det->get_mg_n_in_tree()] = new TCanvas(name.str().c_str());
			//cDisplay[current_det->get_mg_n_in_tree()]->Divide(3);
			name << "_";
			global_signal[current_det->get_mg_n_in_tree()] = new TH1D((name.str() + "signal").c_str(),(name.str() + "signal").c_str(),500,0,4000);
			global_noise[current_det->get_mg_n_in_tree()] = new TH1D((name.str() + "noise").c_str(),(name.str() + "noise").c_str(),100,0,100);
			for(int i=0;i<61;i++){
				global_noise[current_det->get_mg_n_in_tree()]->Fill(current_det->get_RMS(i));
			}
			global_signal_over_noise[current_det->get_mg_n_in_tree()] = new TProfile((name.str() + "S/B").c_str(),(name.str() + "S/B").c_str(),61,0,61);
		}
	}
	long nentries = fChain->GetEntriesFast();
	for(int i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == "MG"){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				MG_Event current_event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent);
				current_event.MultiCluster();
				vector<MG_Cluster> current_cluster = current_event.get_clusters();
				for(vector<MG_Cluster>::iterator jt=current_cluster.begin();jt!=current_cluster.end();++jt){
					double current_signal = jt->get_maxStripAmpl();
					double current_noise = current_det->get_RMS(MG_Detector::StripToChannel[jt->get_maxStrip()]);
					global_signal[current_det->get_mg_n_in_tree()]->Fill(current_signal);
					global_signal_over_noise[current_det->get_mg_n_in_tree()]->Fill(MG_Detector::StripToChannel[jt->get_maxStrip()],current_signal/current_noise);
				}
			}
		}
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	/*
	for(map<int,TCanvas*>::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
		it->second->cd(1);
		TFitResultPtr res_signal = global_signal[it->first]->Fit("landau","SQ");
		global_signal[it->first]->Draw();
		it->second->cd(2);
		TFitResultPtr res_noise = global_noise[it->first]->Fit("gaus","SQ");
		global_noise[it->first]->Draw();
		it->second->cd(3);
		TLine * average_SoN = new TLine(0,(res_signal->Parameter(1))/(res_noise->Parameter(1)),61,(res_signal->Parameter(1))/(res_noise->Parameter(1)));
		average_SoN->SetLineStyle(2);
		average_SoN->SetLineColor(2);
		global_signal_over_noise[it->first]->Draw();
		average_SoN->Draw();
	}
	*/
	return global_signal_over_noise;
}