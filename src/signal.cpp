#define signal_cpp
#include "signal.h"
#include "event.h"
#include "Tanalyse.h"
#include "ray.h"
#include "datareader.h"
#include "point.h"

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <sstream>
#include <limits>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLine.h>
#include <TMath.h>
#include <TF1.h>

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
using std::numeric_limits;

using TMath::CeilNint;
using TMath::Min;
using TMath::Sqrt;
using TMath::ATan;

//boost
using boost::property_tree::ptree;

Signal::Signal(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	electronic_type = Tomography::str_to_elec(config_tree.get<string>("electronic_type"));
	data_file_basename = config_tree.get<string>("data_file_basename");
	signalName = config_tree.get<string>("signal_file");
	PedName = config_tree.get<string>("Ped");
	RMSName = config_tree.get<string>("RMSPed");
	max_event = config_tree.get<long>("max_event");
	data_file_first = config_tree.get<int>("data_file_first");
	data_file_last = config_tree.get<int>("data_file_last");
	cout << signalName << endl;
	cout << RMSName << endl;
	ifstream fIn_test(signalName.c_str());
	exists = fIn_test.good();
	fIn_test.close();
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
	if(exists){
		fIn = new TFile(signalName.c_str(),"READ");
		TTree * treeIn = (TTree*)(fIn->Get("T"));
		Tsignal::Init(treeIn,CM_N,MG_N);
	}
	else{
		Tsignal::Init(0,CM_N,MG_N);
		cout << "Waring, signal file is missing !" << endl;
	}
	analyseTree = config_tree.get<string>("Tree");
	use_srf = config_tree.get<bool>("use_SRF");
}
Signal::~Signal(){
	if(exists){
		delete fIn;
	}
}
void Signal::MultiCluster(){
	cout << "destination file : " << analyseTree << endl;
	Tanalyse * analyseFile = new Tanalyse(analyseTree,CM_N,MG_N);
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);
		vector<MG_Event> mg_events;
		vector<CM_Event> cm_events;

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				mg_events.push_back(MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent));
				mg_events.back().MultiCluster();
			}
			else if((*it)->get_type() == Tomography::CM){
				CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
				cm_events.push_back(CM_Event(*current_det,get_cm_ampl(current_det->get_cm_n_in_tree()),use_srf,Nevent));
				cm_events.back().MultiCluster();
			}
		}
		
		analyseFile->fillTree(Nevent,evttime,mg_events,cm_events);
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	analyseFile->Write();
	analyseFile->CloseFile();
	//delete analyseFile;
}

void Signal::ElecToAnalyse(){
	cout << "Reading Pedestal : " << PedName << endl;
	map<Tomography::det_type,vector<vector<float> > > Pedestal;
	Pedestal[Tomography::CM] = vector<vector<float> >(CM_N,vector<float>(DataReader::Nstrip_CM,0));
	Pedestal[Tomography::MG] = vector<vector<float> >(MG_N,vector<float>(DataReader::Nstrip_MG,0));
	ifstream pedFile(PedName.c_str());
	for(int i=0;i<CM_N;i++){
		for(int j=0;j<DataReader::Nstrip_CM;j++){
			pedFile >> i >> j >> Pedestal[Tomography::CM][i][j];
		}
	}
	for(int i=0;i<MG_N;i++){
		for(int j=0;j<DataReader::Nstrip_MG;j++){
			pedFile >> i >> j >> Pedestal[Tomography::MG][i][j];
		}
	}
	pedFile.close();
	cout << "destination file : " << analyseTree << endl;
	Tanalyse * analyseFile = new Tanalyse(analyseTree,CM_N,MG_N);
	Nevent = 0;
	long event_nb = 0;
	string extension = "";
	DataReader * current_data_reader;
	if(electronic_type == Tomography::Feminos){
		extension = ".aqs";
		current_data_reader = new FeminosDataReader("tmp",det_type_by_asic,det_n_by_asic);
		ostringstream name;
		name << data_file_basename << setfill('0') << setw(3) << data_file_first << extension;
		Nevent = dynamic_cast<FeminosDataReader*>(current_data_reader)->get_first_event_nb(name.str());
	}
	else if(electronic_type == Tomography::Dream){
		extension = "_01.fdf";
		current_data_reader = new DreamDataReader("tmp",det_type_by_asic,det_n_by_asic);
	}
	else return;
	map<Tomography::det_type,int> detector_div;
	detector_div[Tomography::CM] = 2;
	detector_div[Tomography::MG] = 2;
	map<Tomography::det_type,int> Nstrip;
	Nstrip[Tomography::CM] = DataReader::Nstrip_CM;
	Nstrip[Tomography::MG] = DataReader::Nstrip_MG;
	for(int i=data_file_first;i<=data_file_last;i++){
		int event_nb_file = 0;
		ostringstream name;
		name << data_file_basename << setfill('0') << setw(3) << i << extension;
		ifstream data_file(name.str().c_str(),ifstream::binary);
		while(data_file.good()){
			map<Tomography::det_type,vector<vector<vector<double> > > > event_ampl = current_data_reader->read_event(&data_file,Nevent,false);
			if(event_ampl.size()==0) break;
			map<Tomography::det_type,vector<vector<float> > >::iterator ped_it = Pedestal.begin();
			for(map<Tomography::det_type,vector<vector<vector<double> > > >::iterator it = event_ampl.begin();it!=event_ampl.end();++it){
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
				if((*it)->get_type() == Tomography::MG){
					MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
					mg_events.push_back(MG_Event(*current_det,event_ampl[Tomography::MG][current_det->get_mg_n_in_tree()],use_srf,event_nb));
					mg_events.back().MultiCluster();
				}
				else if((*it)->get_type() == Tomography::CM){
					CM_Detector * current_det = dynamic_cast<CM_Detector*>(*it);
					cm_events.push_back(CM_Event(*current_det,event_ampl[Tomography::CM][current_det->get_cm_n_in_tree()],use_srf,event_nb));
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

void Signal::HoughTracking(long event_nb){
	if(CM_N!=0){
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
	map<bool,map<double,vector<MG_Cluster> > > all_cluster;
	vector<Event*> events;
	double max_z = -10000;
	double min_z = 10000;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == Tomography::MG){
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
			all_cluster[current_det->get_is_X()][current_det->get_z()].insert(all_cluster[current_det->get_is_X()][current_det->get_z()].end(),current_cluster.begin(),current_cluster.end());
			events.push_back(new MG_Event(current_event));
			(events.back())->MultiCluster();
		}
		if((*it)->get_z()>max_z) max_z = (*it)->get_z();
		if((*it)->get_z()<min_z) min_z = (*it)->get_z();
	}
	//max_z+=10;
	//min_z-=10;
	int bin_n = 500;
	double min_coord = -6*Tomography::XY_size/10.;
	double max_coord = 6*Tomography::XY_size/10.;
	TH2D * hough_space_X = new TH2D("hough_space_X","hough_space_X",bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
	TH2D * hough_space_Y = new TH2D("hough_space_Y","hough_space_Y",bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
	int suitable_clus_n = 0;
	//draw clusters in hough space
	bin_n = 2*bin_n;
	map<bool,map<double,int> > sizes;
	for(map<bool,map<double,vector<MG_Cluster> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
		for(map<double,vector<MG_Cluster> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			for(vector<MG_Cluster>::iterator it=(kt->second).begin();it!=(kt->second).end();++it){
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
			sizes[jt->first][kt->first] = (kt->second).size();
		}
	}
	TGraph * int_X = new TGraph();
	TGraph * int_Y = new TGraph();
	int X_int_nb = 0;
	int Y_int_nb = 0;
	int_X->SetMarkerStyle(24);
	int_Y->SetMarkerStyle(24);
	int_X->SetMarkerSize(2);
	int_Y->SetMarkerSize(2);
	int_X->SetMarkerColor(2);
	int_Y->SetMarkerColor(2);
	map<bool,Point_2D> hough_ray;
	for(map<bool,map<double,vector<MG_Cluster> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
		vector<map<double,int> > comb = CosmicBenchEvent::combinaisons(sizes[jt->first]);
		double smallest_distance = numeric_limits<double>::max();
		bool found = false;
		Point_2D best_intersection;
		for(vector<map<double,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
			vector<Point_2D> intersections;
			for(map<double,int>::iterator map_it = kt->begin();map_it!=kt->end();++map_it){
				Line_2D first_line;
				MG_Cluster first_cluster = (jt->second)[map_it->first][map_it->second];
				if(!(first_cluster.get_is_up())) first_line = Line_2D(Point_2D(min_coord + (first_cluster.get_pos_mm() - min_coord)*(min_z-max_z)/(first_cluster.get_z()-max_z),min_coord),Point_2D(max_coord + (first_cluster.get_pos_mm() - max_coord)*(min_z-max_z)/(first_cluster.get_z()-max_z),max_coord));
				else first_line = Line_2D(Point_2D(min_coord,min_coord + (first_cluster.get_pos_mm() - min_coord)*(max_z-min_z)/(first_cluster.get_z()-min_z)),Point_2D(max_coord,max_coord + (first_cluster.get_pos_mm() - max_coord)*(max_z-min_z)/(first_cluster.get_z()-min_z)));
				map<double,int>::iterator map_jt = map_it;
				map_jt++;
				while(map_jt!=kt->end()){
					Line_2D second_line;
					MG_Cluster second_cluster = (jt->second)[map_jt->first][map_jt->second];
					if(!(second_cluster.get_is_up())) second_line = Line_2D(Point_2D(min_coord + (second_cluster.get_pos_mm() - min_coord)*(min_z-max_z)/(second_cluster.get_z()-max_z),min_coord),Point_2D(max_coord + (second_cluster.get_pos_mm() - max_coord)*(min_z-max_z)/(second_cluster.get_z()-max_z),max_coord));
					else second_line = Line_2D(Point_2D(min_coord,min_coord + (second_cluster.get_pos_mm() - min_coord)*(max_z-min_z)/(second_cluster.get_z()-min_z)),Point_2D(max_coord,max_coord + (second_cluster.get_pos_mm() - max_coord)*(max_z-min_z)/(second_cluster.get_z()-min_z)));
					intersections.push_back(first_line.intersection(second_line));
					++map_jt;
				}
			}
			if(intersections.size()>1){
				double biggest_distance = 0;
				for(vector<Point_2D>::iterator vec_it = intersections.begin();vec_it!=intersections.end();++vec_it){
					vector<Point_2D>::iterator vec_jt = vec_it;
					vec_jt++;
					while(vec_jt!=intersections.end()){
						double current_distance = ((*vec_it) - (*vec_jt)).norm();
						if(current_distance> biggest_distance) biggest_distance = current_distance;
						++vec_jt;
					}
				}
				if(biggest_distance<smallest_distance){
					smallest_distance = biggest_distance;
					found = true;
					vector<Point_2D>::iterator vec_it = intersections.begin();
					best_intersection = *vec_it;
					++vec_it;
					while(vec_it!=intersections.end()){
						best_intersection += *vec_it;
						++vec_it;
					}
					best_intersection /= intersections.size();
				}
			}
		}
		if(found){
			hough_ray[jt->first] = best_intersection;
			if(jt->first){
				int_X->SetPoint(X_int_nb,best_intersection.get_X(),best_intersection.get_Y());
				X_int_nb++;
			}
			else{
				int_Y->SetPoint(Y_int_nb,best_intersection.get_X(),best_intersection.get_Y());
				Y_int_nb++;
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
	gStyle->SetNumberContours(512);
	TCanvas * cHough = new TCanvas();
	cHough->Divide(2);
	cHough->cd(1);
	hough_space_X->Draw("colz");
	if(rays.size()>0) rays_X->Draw("sameP");
	if(X_int_nb>0) int_X->Draw("sameP");
	cHough->cd(2);
	hough_space_Y->Draw("colz");
	if(rays.size()>0) rays_Y->Draw("sameP");
	if(Y_int_nb>0) int_Y->Draw("sameP");
}
map<int,TProfile*> Signal::SignalOverNoise(){
	map<int,TProfile*> global_signal_over_noise;
	if(CM_N!=0){
		cout << "not implemented with CM" << endl;
		return global_signal_over_noise;
	}
	map<int,TCanvas*> cDisplay;
	map<int,TH1D*> global_signal;
	map<int,TH1D*> global_noise;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == Tomography::MG){
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
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
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
	return global_signal_over_noise;
}

void Signal::SignalOverNoiseDisplay(){
	map<int,TProfile*> global_signal_over_noise;
	if(CM_N!=0){
		cout << "not implemented with CM" << endl;
		return;
	}
	map<int,TCanvas*> cDisplay;
	map<int,TH1D*> global_signal;
	map<int,TH1D*> global_noise;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_type() == Tomography::MG){
			MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
			ostringstream name;
			name << "Multigen_" << current_det->get_mg_n_in_tree();
			cDisplay[current_det->get_mg_n_in_tree()] = new TCanvas(name.str().c_str());
			cDisplay[current_det->get_mg_n_in_tree()]->Divide(3);
			name << "_";
			global_signal[current_det->get_mg_n_in_tree()] = new TH1D((name.str() + "signal").c_str(),(name.str() + "signal").c_str(),500,0,4000);
			global_noise[current_det->get_mg_n_in_tree()] = new TH1D((name.str() + "noise").c_str(),(name.str() + "noise").c_str(),100,0,100);
			for(int i=0;i<61;i++){
				global_noise[current_det->get_mg_n_in_tree()]->Fill(current_det->get_RMS(i));
			}
			global_signal_over_noise[current_det->get_mg_n_in_tree()] = new TProfile((name.str() + "S/B").c_str(),(name.str() + "S/B").c_str(),61,0,61);
		}
	}
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
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
		if(i%5000 == 0 && Tomography::live_graphic_display){
			for(map<int,TCanvas*>::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
				it->second->cd(1);
				global_signal[it->first]->Draw();
				it->second->cd(2);
				global_noise[it->first]->Draw();
				it->second->cd(3);
				global_signal_over_noise[it->first]->Draw();
				it->second->Modified();
				it->second->Update();
			}
		}
	}
	cout << "\r" << nentries << "/" << nentries << endl;
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
		TLine * mean_SoN = new TLine(0,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()),61,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()));
		mean_SoN->SetLineStyle(2);
		mean_SoN->SetLineColor(4);
		global_signal_over_noise[it->first]->Draw();
		average_SoN->Draw();
		mean_SoN->Draw();
		it->second->Modified();
		it->second->Update();
		cout << Tomography::MG << it->first << endl;
		cout << "    mean S/B : " << global_signal_over_noise[it->first]->GetMean(2) << endl;
		cout << "    mean S/mean B : " << (global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()) << endl;
		cout << "    sigma S/B : " << global_signal_over_noise[it->first]->GetRMS(2) << endl;
		cout << "    delta mean S/B : " << global_signal_over_noise[it->first]->GetMean(11) << endl;
	}
}

void Signal::EventDisplay(int evn_min, int evn_max){
	int column_nb = CeilNint((MG_N+CM_N)/2.);
	TCanvas * cDisplay = new TCanvas("cDisplay","cDisplay",800,600);
	cDisplay->Divide(column_nb,2);
	vector<TGraph*> signal_shape;
	for(int i=0;i<((CM_N*DataReader::Nstrip_CM)+(MG_N*DataReader::Nstrip_MG));i++){
		signal_shape.push_back(new TGraph());
	}
	long nentries = Min(fChain->GetEntriesFast(),static_cast<Long64_t>(evn_max));
	if(evn_min>nentries) return;
	for(long i=evn_min;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);
		for(int j=0;j<CM_N;j++){
			for(int k=0;k<DataReader::Nstrip_CM;k++){
				int index = k+(j*DataReader::Nstrip_CM);
				int point_nb = 0;
				for(int l=0;l<DataReader::Nsample;l++){
					signal_shape[index]->SetPoint(point_nb,l,StripAmpl_CM_corr[j][k][l]);
					point_nb++;
				}
				cDisplay->cd(j+1);
				if(k==0){
					signal_shape[index]->GetHistogram()->SetMinimum(-100);
					signal_shape[index]->GetHistogram()->SetMaximum(1200);
					signal_shape[index]->Draw("AL");
				}
				else signal_shape[index]->Draw("L");
			}
		}
		for(int j=0;j<MG_N;j++){
			for(int k=0;k<DataReader::Nstrip_MG;k++){
				int index = k + (j*DataReader::Nstrip_MG) + (CM_N*DataReader::Nstrip_CM);
				int point_nb = 0;
				for(int l=0;l<DataReader::Nsample;l++){
					signal_shape[index]->SetPoint(point_nb,l,StripAmpl_MG_corr[j][k][l]);
					point_nb++;
				}
				cDisplay->cd(CM_N+j+1);
				if(k==0){
					signal_shape[index]->GetHistogram()->SetMinimum(-100);
					signal_shape[index]->GetHistogram()->SetMaximum(1200);
					signal_shape[index]->Draw("AL");
				}
				else signal_shape[index]->Draw("L");
			}
		}
		cDisplay->Modified();
		cDisplay->Update();
		gSystem->Sleep(1000);
	}
}

void Signal::SignalDispersion(){
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * cDisplay = new TCanvas("cDisplay");
	cDisplay->Divide(2,2);
	TCanvas * cControl = new TCanvas("cControl");
	cControl->Divide(3);
	TCanvas * cAnalyse = new TCanvas("cAnalyse");
	cAnalyse->Divide(2);
	int nbin_time = Tomography::Nsample;
	TH2D * signalShape_X_pos = new TH2D("signalShape_X_pos","signalShape_X_pos",60,-6,6,nbin_time,0,nbin_time);
	TH2D * signalShape_Y_pos = new TH2D("signalShape_Y_pos","signalShape_Y_pos",60,-15,15,nbin_time,0,nbin_time);
	TH2D * signalShape_X_neg = new TH2D("signalShape_X_neg","signalShape_X_neg",60,-6,6,nbin_time,0,nbin_time);
	TH2D * signalShape_Y_neg = new TH2D("signalShape_Y_neg","signalShape_Y_neg",60,-15,15,nbin_time,0,nbin_time);
	TGraph * mean_X_pos = new TGraph();
	TGraph * mean_Y_pos = new TGraph();
	TGraph * mean_X_neg = new TGraph();
	TGraph * mean_Y_neg = new TGraph();
	mean_X_pos->SetMarkerStyle(2);
	mean_Y_pos->SetMarkerStyle(2);
	mean_X_neg->SetMarkerStyle(2);
	mean_Y_neg->SetMarkerStyle(2);

	TH2D * angle_shape_corr_X = new TH2D("angle_shape_corr_X","angle_shape_corr_X",50,-0.5,0.5,50,-1,1);
	TH2D * angle_shape_corr_Y = new TH2D("angle_shape_corr_Y","angle_shape_corr_Y",50,-0.5,0.5,50,-1,1);

	TH1D * clus_pos = new TH1D("clus_pos","clus_pos",1024,0,1023);
	TH1D * clus_size = new TH1D("clus_size","clus_size",50,-1,47);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	double min_angle = -1;
	double max_angle = 10;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	if(CM_N>0){
		cout << "does not support CM detectors" << endl;
		return;
	}
	for(long ientry=0;ientry<nentries;ientry++){
		LoadTree(ientry);
		GetEntry(ientry);
		/*
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				MG_Event current_event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent);
				current_event.MultiCluster();
				current_event.do_cuts();
				int current_n = current_event.get_n_in_tree();
				bool current_X = current_event.get_is_X();
				vector<MG_Cluster> current_clusters = current_event.get_clusters();
				for(vector<MG_Cluster>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
					clus_pos->Fill(clus_it->get_pos());
					clus_size->Fill(clus_it->get_size());
					int left = ((clus_it->get_pos() - clus_it->get_size()) > 0) ? (clus_it->get_pos() - clus_it->get_size()) : 0;
					int right = ((clus_it->get_pos() + clus_it->get_size()) < 1023) ? (clus_it->get_pos() + clus_it->get_size()) : 1023;
					for(int i=left;i<=right;i++){
						for(int j=0;j<Tomography::Nsample;j++){
							if(current_X){
								signalShape_X->Fill(i - clus_it->get_pos(),j,StripAmpl_MG_corr[current_n][MG_Detector::StripToChannel[i]][j]);
							}
							else{
								signalShape_Y->Fill(i - clus_it->get_pos(),j,StripAmpl_MG_corr[current_n][MG_Detector::StripToChannel[i]][j]);
							}
						}
					}
				}
			}
		}
		*/

		vector<Event*> events;
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				events.push_back(new MG_Event(*current_det,get_mg_ampl(current_det->get_mg_n_in_tree()),use_srf,Nevent));
				events.back()->MultiCluster();
				events.back()->do_cuts();
			}
		}
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,events);
		vector<Ray> current_rays = CBEvent->get_absorption_rays();
		for(vector<Ray>::iterator ray_it = current_rays.begin();ray_it!=current_rays.end();++ray_it){
			double slope = ATan(Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X())));
			if(slope<min_angle || slope>max_angle) continue;
			ray_slope->Fill(slope);
			vector<Cluster*> current_clusters = ray_it->get_clus();
			for(vector<Cluster*>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				clus_pos->Fill((*clus_it)->get_pos());
				clus_size->Fill((*clus_it)->get_size());
				int left = (((*clus_it)->get_pos() - (*clus_it)->get_size()) > 0) ? ((*clus_it)->get_pos() - (*clus_it)->get_size()) : 0;
				int right = (((*clus_it)->get_pos() + (*clus_it)->get_size()) < 1023) ? ((*clus_it)->get_pos() + (*clus_it)->get_size()) : 1023;
				int current_n = (*clus_it)->get_n_in_tree();
				int current_X = (*clus_it)->get_is_X();
				TGraph * gFit = new TGraph();
				TF1 * lFit = new TF1("lfit","pol1(0)",5,15);
				lFit->SetParameters(0,0);
				lFit->SetParLimits(0,-5,5);
				lFit->SetParLimits(1,-1,1);
				for(int j=0;j<Tomography::Nsample;j++){
					double mean = 0;
					double ampl = 0;
					for(int i=left;i<=right;i++){
						double current_ampl = StripAmpl_MG_corr[current_n][MG_Detector::StripToChannel[i]][j];
						double current_pos = i - (*clus_it)->get_pos();
						mean = (mean*ampl + current_pos*current_ampl)/(ampl+current_ampl);
						ampl += current_ampl;
						if(current_X){
							if(ray_it->get_slope_X() > 0) signalShape_X_pos->Fill(current_pos,j,current_ampl);
							else signalShape_X_neg->Fill(current_pos,j,current_ampl);
						}
						else{
							if(ray_it->get_slope_Y() > 0) signalShape_Y_pos->Fill(current_pos,j,current_ampl);
							else signalShape_Y_neg->Fill(current_pos,j,current_ampl);
						}
					}
					gFit->SetPoint(j,j,mean);
				}
				gFit->Fit(lFit,"QNR");
				if(current_X) angle_shape_corr_X->Fill(ray_it->get_slope_X(), lFit->GetParameter(1));
				else angle_shape_corr_Y->Fill(ray_it->get_slope_Y(), lFit->GetParameter(1));
				delete gFit; delete lFit;
				delete *clus_it;
			}
		}
		delete CBEvent;
		for(vector<Event*>::iterator ev_it = events.begin();ev_it!=events.end();++ev_it){
			delete *ev_it;
		}
		events.clear();
		if(ientry%100 == 0) cout << "\r" << ientry << "/" << nentries << flush;
		if((ientry%5000 == 0) && Tomography::live_graphic_display){
			for(int i=1;i<=nbin_time;i++){
				mean_X_pos->SetPoint(i-1,signalShape_X_pos->ProjectionX("_px",i,i)->GetMean(),i-0.5);
				mean_Y_pos->SetPoint(i-1,signalShape_Y_pos->ProjectionX("_px",i,i)->GetMean(),i-0.5);
				mean_X_neg->SetPoint(i-1,signalShape_X_neg->ProjectionX("_px",i,i)->GetMean(),i-0.5);
				mean_Y_neg->SetPoint(i-1,signalShape_Y_neg->ProjectionX("_px",i,i)->GetMean(),i-0.5);
			}
			cDisplay->cd(1);
			signalShape_X_pos->Draw("COLZ");
			mean_X_pos->Draw("PSAME");
			cDisplay->cd(2);
			signalShape_Y_pos->Draw("COLZ");
			mean_Y_pos->Draw("PSAME");
			cDisplay->cd(3);
			signalShape_X_neg->Draw("COLZ");
			mean_X_neg->Draw("PSAME");
			cDisplay->cd(4);
			signalShape_Y_neg->Draw("COLZ");
			mean_Y_neg->Draw("PSAME");
			cDisplay->Modified();
			cDisplay->Update();
			cControl->cd(1);
			clus_pos->Draw();
			cControl->cd(2);
			clus_size->Draw();
			cControl->cd(3);
			ray_slope->Draw();
			cControl->Modified();
			cControl->Update();
			cAnalyse->cd(1);
			angle_shape_corr_X->Draw("COLZ");
			cAnalyse->cd(2);
			angle_shape_corr_Y->Draw("COLZ");
			cAnalyse->Modified();
			cAnalyse->Update();
		}
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	for(int i=1;i<=nbin_time;i++){
		mean_X_pos->SetPoint(i-1,signalShape_X_pos->ProjectionX("_px",i,i)->GetMean(),i);
		mean_Y_pos->SetPoint(i-1,signalShape_Y_pos->ProjectionX("_px",i,i)->GetMean(),i);
		mean_X_neg->SetPoint(i-1,signalShape_X_neg->ProjectionX("_px",i,i)->GetMean(),i);
		mean_Y_neg->SetPoint(i-1,signalShape_Y_neg->ProjectionX("_px",i,i)->GetMean(),i);
	}
	cDisplay->cd(1);
	signalShape_X_pos->Draw("COLZ");
	mean_X_pos->Draw("PSAME");
	cDisplay->cd(2);
	signalShape_Y_pos->Draw("COLZ");
	mean_Y_pos->Draw("PSAME");
	cDisplay->cd(3);
	signalShape_X_neg->Draw("COLZ");
	mean_X_neg->Draw("PSAME");
	cDisplay->cd(4);
	signalShape_Y_neg->Draw("COLZ");
	mean_Y_neg->Draw("PSAME");
	cDisplay->Modified();
	cDisplay->Update();
	cControl->cd(1);
	clus_pos->Draw();
	cControl->cd(2);
	clus_size->Draw();
	cControl->cd(3);
	ray_slope->Draw();
	cControl->Modified();
	cControl->Update();
	cAnalyse->cd(1);
	angle_shape_corr_X->Draw("COLZ");
	cAnalyse->cd(2);
	angle_shape_corr_Y->Draw("COLZ");
	cAnalyse->Modified();
	cAnalyse->Update();
	cout << endl;
}