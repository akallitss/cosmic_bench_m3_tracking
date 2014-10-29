#define signal_cpp
#include "signal.h"
#include "event.h"
#include "Tanalyse.h"

#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <utility>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

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

//boost
using boost::property_tree::ptree;

Signal::Signal(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	cout << config_tree.get<string>("signal_file") << endl;
	cout << config_tree.get<string>("RMSPed") << endl;
	CM_n = 0;
	MG_n = 0;
	TFile *fIn = new TFile((config_tree.get<string>("signal_file")).c_str(),"READ");
	TTree * treeIn = (TTree*)(fIn->Get("T"));
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
		CM_n++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		detectors.push_back(new MG_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("is_up"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle")));
		detectors.back()->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		detectors.back()->set_RMS(RMS[total_CM_N+child.second.get<int>("mg_n")]);
		(dynamic_cast<MG_Detector*>(detectors.back()))->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
		(dynamic_cast<MG_Detector*>(detectors.back()))->set_SRF(child.second.get<double>("SRF.offset"),child.second.get<double>("SRF.gauss"),child.second.get<double>("SRF.lorentz"),child.second.get<double>("SRF.ratio"));
		MG_n++;
	}
	if((total_CM_N!=CM_n) || (total_MG_N!=MG_n)){
		cout << "problem in detectors number" << endl;
		return;
	}
	if(total_CM_N!=0) cout << "warning, CosMultis are not fully supported !" << endl;
	Init(treeIn,CM_n,MG_n);
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
				mg_events.push_back(MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf));
				mg_events.back().MultiCluster();
			}
			else if((*it)->get_type() == "CM"){
				CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
				cm_events.push_back(CM_Event(*current_det,get_cm_ampl(current_det->get_cm_n_in_tree()),use_srf));
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
	double max_z = -10000;
	double min_z = 10000;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == "MG"){
			MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
			MG_Event current_event = MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf);
			current_event.HoughCluster();
			vector<MG_Cluster> current_cluster = current_event.get_clusters();
			all_cluster.insert(all_cluster.end(),current_cluster.begin(),current_cluster.end());
		}
		if((*it)->get_z()>max_z) max_z = (*it)->get_z();
		if((*it)->get_z()<min_z) min_z = (*it)->get_z();
	}
	max_z+=10;
	min_z-=10;
	int bin_n = 500;
	double min_coord = -100;
	double max_coord = 600;
	TH2D * hough_space_X = new TH2D("hough_space_X","hough_space_X",2*bin_n,min_coord,max_coord,2*bin_n,min_coord,max_coord);
	TH2D * hough_space_Y = new TH2D("hough_space_Y","hough_space_Y",2*bin_n,min_coord,max_coord,2*bin_n,min_coord,max_coord);
	int suitable_clus_n = 0;
	for(vector<MG_Cluster>::iterator it=all_cluster.begin();it!=all_cluster.end();++it){
		if(it->get_size()<2) continue;
		if(!(it->is_suitable(dynamic_cast<MG_Detector*>(detectors[it->find_det(detectors)])))) continue;
		suitable_clus_n++;
		for(int i=0;i<bin_n;i++){
			double current_coord_up = min_coord +i*max_coord/bin_n;
			double current_coord_down = current_coord_up + (it->get_pos_mm() - current_coord_up)*(min_z-max_z)/(it->get_z()-max_z);
			if(it->get_is_X()){
				hough_space_X->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
			}
			else{
				hough_space_Y->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
			}
		}
	}
	cout << "total number of suitable cluster : " << suitable_clus_n << endl;
	gStyle->SetPalette(55,0);
	TCanvas * cHough = new TCanvas();
	cHough->Divide(2);
	cHough->cd(1);
	hough_space_X->Draw("colz");
	cHough->cd(2);
	hough_space_Y->Draw("colz");
}