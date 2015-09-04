#define Signal_cpp
#include "Signal.h"
#include "event.h"
#include "Tanalyse_W.h"
#include "ray.h"
#include "datareader.h"
#include "point.h"
#include "Tray.h"
#include "detector.h"

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
#include <TGraphErrors.h>
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
using TMath::Abs;
using TMath::Pi;

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
	for(vector<Detector*>::const_iterator det_it=detectors.begin();det_it!=detectors.end();++det_it){
		det_type_by_asic[(*det_it)->get_asic_n()] = (*det_it)->get_type();
		det_n_by_asic[(*det_it)->get_asic_n()] = (*det_it)->get_n_in_tree();
	}
	if(exists){
		fIn = new TFile(signalName.c_str(),"READ");
		TTree * treeIn = (TTree*)(fIn->Get("T"));
		Tsignal_R::Init(treeIn,det_n);
	}
	else{
		Tsignal_R::Init(0,det_n);
		cout << "Waring, signal file is missing !" << endl;
	}
	analyseTree = config_tree.get<string>("Tree");
	use_srf = config_tree.get<bool>("use_SRF");
}

Signal::Signal(ptree config_tree){
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
	for(vector<Detector*>::const_iterator det_it=detectors.begin();det_it!=detectors.end();++det_it){
		det_type_by_asic[(*det_it)->get_asic_n()] = (*det_it)->get_type();
		det_n_by_asic[(*det_it)->get_asic_n()] = (*det_it)->get_n_in_tree();
	}
	if(exists){
		fIn = new TFile(signalName.c_str(),"READ");
		TTree * treeIn = (TTree*)(fIn->Get("T"));
		Tsignal_R::Init(treeIn,det_n);
	}
	else{
		Tsignal_R::Init(0,det_n);
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
	Tanalyse_W * analyseFile = new Tanalyse_W(analyseTree,det_n);
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Nevent",1);
	fChain->SetBranchStatus("evttime",1);
	fChain->SetBranchStatus("StripAmpl_*_corr",1);
	for(long i=0;i<nentries && Tomography::can_continue;i++){
		LoadTree(i);
		GetEntry(i);
		map<Tomography::det_type,vector<Event*> > events;

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			events[(*it)->get_type()].push_back((*it)->build_event(get_ampl((*it)->get_type(),(*it)->get_n_in_tree()),Nevent));
			(events[(*it)->get_type()].back())->MultiCluster();
		}
		
		analyseFile->fillTree(Nevent,evttime,events);
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
		for(map<Tomography::det_type,vector<Event*> >::iterator type_it = events.begin();type_it!=events.end();++type_it){
			for(vector<Event*>::iterator ev_it = (type_it->second).begin();ev_it!=(type_it->second).end();++ev_it){
				delete *ev_it;
			}
		}
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	analyseFile->Write();
	analyseFile->CloseFile();
	//delete analyseFile;
	fChain->SetBranchStatus("*",1);
}

void Signal::HoughTracking(long event_nb){
	long nentries = fChain->GetEntriesFast();
	if(event_nb<0 || event_nb>=nentries){
		cout << "invalid event number" << endl;
		return;
	}
	LoadTree(event_nb);
	GetEntry(event_nb);
	map<bool,map<double,vector<Cluster*> > > all_cluster;
	vector<Event*> events;
	double max_z = -10000;
	double min_z = 10000;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		Event * current_event = (*it)->build_event(get_ampl((*it)->get_type(),(*it)->get_n_in_tree()),Nevent);
		current_event->HoughCluster();
		vector<Cluster*> current_cluster = current_event->get_clusters();
		vector<Cluster*>::iterator clus_it = current_cluster.begin();
		while(clus_it!=current_cluster.end()){
			if(!((*it)->is_suitable(*clus_it))){
				delete *clus_it;
				clus_it = current_cluster.erase(clus_it);
			}
			else ++clus_it;
		}
		cout << "number of suitable cluster for " << (*it)->get_type() << "_" << (*it)->get_n_in_tree() << " : " << current_cluster.size();
		all_cluster[(*it)->get_is_X()][(*it)->get_z()].insert(all_cluster[(*it)->get_is_X()][(*it)->get_z()].end(),current_cluster.begin(),current_cluster.end());
		current_event->MultiCluster();
		cout << " (" << current_event->get_NClus() << ")" << endl;
		events.push_back(current_event);
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
	for(map<bool,map<double,vector<Cluster*> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
		for(map<double,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
			for(vector<Cluster*>::iterator it=(kt->second).begin();it!=(kt->second).end();++it){
				suitable_clus_n++;
				if(!(*it)->get_is_up()){
					for(int i=0;i<bin_n;i++){
						double current_coord_up = min_coord +i*(max_coord-min_coord)/bin_n;
						double current_coord_down = current_coord_up + ((*it)->get_pos_mm() - current_coord_up)*(min_z-max_z)/((*it)->get_z()-max_z);
						if((*it)->get_is_X()){
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
						double current_coord_up = current_coord_down + ((*it)->get_pos_mm() - current_coord_down)*(max_z-min_z)/((*it)->get_z()-min_z);
						if((*it)->get_is_X()){
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
	for(int drop = 0;drop<2;drop++){
		for(map<bool,map<double,vector<Cluster*> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
			if(hough_ray.count(jt->first)>0) continue;
			vector<map<double,int> > comb = CosmicBenchEvent::combinaisons(sizes[jt->first], (drop>0));
			double smallest_distance = numeric_limits<double>::max();
			bool found = false;
			Point_2D best_intersection;
			for(vector<map<double,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
				vector<Point_2D> intersections;
				for(map<double,int>::iterator map_it = kt->begin();map_it!=kt->end();++map_it){
					Line_2D first_line;
					Cluster * first_cluster = (jt->second)[map_it->first][map_it->second];
					if(!(first_cluster->get_is_up())) first_line = Line_2D(Point_2D(min_coord + (first_cluster->get_pos_mm() - min_coord)*(min_z-max_z)/(first_cluster->get_z()-max_z),min_coord),Point_2D(max_coord + (first_cluster->get_pos_mm() - max_coord)*(min_z-max_z)/(first_cluster->get_z()-max_z),max_coord));
					else first_line = Line_2D(Point_2D(min_coord,min_coord + (first_cluster->get_pos_mm() - min_coord)*(max_z-min_z)/(first_cluster->get_z()-min_z)),Point_2D(max_coord,max_coord + (first_cluster->get_pos_mm() - max_coord)*(max_z-min_z)/(first_cluster->get_z()-min_z)));
					map<double,int>::iterator map_jt = map_it;
					map_jt++;
					while(map_jt!=kt->end()){
						Line_2D second_line;
						Cluster * second_cluster = (jt->second)[map_jt->first][map_jt->second];
						if(!(second_cluster->get_is_up())) second_line = Line_2D(Point_2D(min_coord + (second_cluster->get_pos_mm() - min_coord)*(min_z-max_z)/(second_cluster->get_z()-max_z),min_coord),Point_2D(max_coord + (second_cluster->get_pos_mm() - max_coord)*(min_z-max_z)/(second_cluster->get_z()-max_z),max_coord));
						else second_line = Line_2D(Point_2D(min_coord,min_coord + (second_cluster->get_pos_mm() - min_coord)*(max_z-min_z)/(second_cluster->get_z()-min_z)),Point_2D(max_coord,max_coord + (second_cluster->get_pos_mm() - max_coord)*(max_z-min_z)/(second_cluster->get_z()-min_z)));
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
					if(biggest_distance<smallest_distance && biggest_distance<10){
						smallest_distance = biggest_distance;
						found = true;
						vector<Point_2D>::iterator vec_it = intersections.begin();
						best_intersection = *vec_it;
						++vec_it;
						while(vec_it!=intersections.end()){
							best_intersection += (*vec_it);
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
	}
	cout << "total number of suitable cluster : " << suitable_clus_n << endl;
	//draw "normal" ray in hough space
	CosmicBenchEvent * thisEvent = new CosmicBenchEvent(this,events);
	vector<Ray> rays = thisEvent->get_absorption_rays();
	delete thisEvent;
	for(vector<Event*>::iterator ev_it=events.begin();ev_it!=events.end();++ev_it){
		delete *ev_it;
	}
	events.clear();
	for(map<bool,map<double,vector<Cluster*> > >::iterator coord_it=all_cluster.begin();coord_it!=all_cluster.end();++coord_it){
		for(map<double,vector<Cluster*> >::iterator alt_it=(coord_it->second).begin();alt_it!=(coord_it->second).end();++alt_it){
			for(vector<Cluster*>::iterator clus_it=(alt_it->second).begin();clus_it!=(alt_it->second).end();++clus_it){
				delete *clus_it;
			}
			(alt_it->second).clear();
		}
		(coord_it->second).clear();
	}
	all_cluster.clear();
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
	cHough->Modified();
	cHough->Update();
}
map<Tomography::det_type,map<int,TProfile*> > Signal::SignalOverNoise(){
	map<Tomography::det_type,map<int,TProfile*> > global_signal_over_noise;
	map<Tomography::det_type,map<int,TH1D*> > global_signal;
	map<Tomography::det_type,map<int,TH1D*> > global_noise;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		ostringstream name;
		name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
		name << "_";
		global_signal[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH1D((name.str() + "signal").c_str(),(name.str() + "signal").c_str(),500,0,Tomography::ADC_max);
		global_noise[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH1D((name.str() + "noise").c_str(),(name.str() + "noise").c_str(),100,0,100);
		for(int i=0;i<(*it)->get_Nchannel();i++){
			global_noise[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill((*it)->get_RMS(i));
		}
		global_signal_over_noise[(*it)->get_type()][(*it)->get_n_in_tree()] = new TProfile((name.str() + "S/B").c_str(),(name.str() + "S/B").c_str(),(*it)->get_Nchannel(),0,(*it)->get_Nchannel());
	}
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries && Tomography::can_continue;i++){
		LoadTree(i);
		GetEntry(i);
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			Event * current_event = (*it)->build_event(get_ampl((*it)->get_type(),(*it)->get_n_in_tree()),Nevent);
			current_event->MultiCluster();
			vector<Cluster*> current_cluster = current_event->get_clusters();
			for(vector<Cluster*>::iterator jt=current_cluster.begin();jt!=current_cluster.end();++jt){
				double current_signal = (*jt)->get_maxStripAmpl();
				double current_noise = (*it)->get_RMS((*it)->StripToChannel((*jt)->get_maxStrip()));
				global_signal[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill(current_signal);
				global_signal_over_noise[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill((*it)->StripToChannel((*jt)->get_maxStrip()),current_signal/current_noise);
				delete *jt;
			}
			delete current_event;
		}
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	return global_signal_over_noise;
}

void Signal::SignalOverNoiseDisplay(){
	map<Tomography::det_type,map<int,TProfile*> > global_signal_over_noise;
	map<Tomography::det_type,map<int,TCanvas*> > cDisplay;
	map<Tomography::det_type,map<int,TH1D*> > global_signal;
	map<Tomography::det_type,map<int,TH1D*> > global_noise;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		ostringstream name;
		name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
		cDisplay[(*it)->get_type()][(*it)->get_n_in_tree()] = new TCanvas(name.str().c_str());
		cDisplay[(*it)->get_type()][(*it)->get_n_in_tree()]->Divide(3);
		name << "_";
		global_signal[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH1D((name.str() + "signal").c_str(),(name.str() + "signal").c_str(),500,0,Tomography::ADC_max);
		global_noise[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH1D((name.str() + "noise").c_str(),(name.str() + "noise").c_str(),100,0,100);
		for(int i=0;i<(*it)->get_Nchannel();i++){
			global_noise[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill((*it)->get_RMS(i));
		}
		global_signal_over_noise[(*it)->get_type()][(*it)->get_n_in_tree()] = new TProfile((name.str() + "S/B").c_str(),(name.str() + "S/B").c_str(),(*it)->get_Nchannel(),0,(*it)->get_Nchannel());
	}
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries && Tomography::can_continue;i++){
		LoadTree(i);
		GetEntry(i);
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			Event * current_event = (*it)->build_event(get_ampl((*it)->get_type(),(*it)->get_n_in_tree()),Nevent);
			current_event->MultiCluster();
			vector<Cluster*> current_cluster = current_event->get_clusters();
			for(vector<Cluster*>::iterator jt=current_cluster.begin();jt!=current_cluster.end();++jt){
				double current_signal = (*jt)->get_maxStripAmpl();
				double current_noise = (*it)->get_RMS((*it)->StripToChannel((*jt)->get_maxStrip()));
				global_signal[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill(current_signal);
				global_signal_over_noise[(*it)->get_type()][(*it)->get_n_in_tree()]->Fill((*it)->StripToChannel((*jt)->get_maxStrip()),current_signal/current_noise);
				delete *jt;
			}
			delete current_event;
		}
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
		if(i%5000 == 0 && Tomography::live_graphic_display){
			for(map<Tomography::det_type,map<int,TCanvas*> >::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
				for(map<int,TCanvas*>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
					jt->second->cd(1);
					global_signal[it->first][jt->first]->Draw();
					jt->second->cd(2);
					global_noise[it->first][jt->first]->Draw();
					jt->second->cd(3);
					global_signal_over_noise[it->first][jt->first]->Draw();
					jt->second->Modified();
					jt->second->Update();
				}
			}
		}
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	for(map<Tomography::det_type,map<int,TCanvas*> >::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
		for(map<int,TCanvas*>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			jt->second->cd(1);
			TFitResultPtr res_signal = global_signal[it->first][jt->first]->Fit("landau","SQ");
			global_signal[it->first][jt->first]->Draw();
			jt->second->cd(2);
			TFitResultPtr res_noise = global_noise[it->first][jt->first]->Fit("gaus","SQ");
			global_noise[it->first][jt->first]->Draw();
			jt->second->cd(3);
			TLine * average_SoN = new TLine(0,(res_signal->Parameter(1))/(res_noise->Parameter(1)),Tomography::Static_Detector[it->first]->get_Nchannel(),(res_signal->Parameter(1))/(res_noise->Parameter(1)));
			average_SoN->SetLineStyle(2);
			average_SoN->SetLineColor(2);
			TLine * mean_SoN = new TLine(0,(global_signal[it->first][jt->first]->GetMean())/(global_noise[it->first][jt->first]->GetMean()),Tomography::Static_Detector[it->first]->get_Nchannel(),(global_signal[it->first][jt->first]->GetMean())/(global_noise[it->first][jt->first]->GetMean()));
			mean_SoN->SetLineStyle(2);
			mean_SoN->SetLineColor(4);
			global_signal_over_noise[it->first][jt->first]->Draw();
			average_SoN->Draw();
			mean_SoN->Draw();
			jt->second->Modified();
			jt->second->Update();
			cout << it->first << "_" << jt->first << endl;
			cout << "    mean S/B : " << global_signal_over_noise[it->first][jt->first]->GetMean(2) << endl;
			cout << "    mean S/mean B : " << (global_signal[it->first][jt->first]->GetMean())/(global_noise[it->first][jt->first]->GetMean()) << endl;
			cout << "    sigma S/B : " << global_signal_over_noise[it->first][jt->first]->GetRMS(2) << endl;
			cout << "    delta mean S/B : " << global_signal_over_noise[it->first][jt->first]->GetMean(11) << endl;
		}
	}
}

void Signal::EventDisplay(int evn_min, int evn_max){
	int column_nb = CeilNint((get_det_N_tot())/2.);
	TCanvas * cDisplay = new TCanvas("cDisplay","cDisplay",800,600);
	cDisplay->Divide(column_nb,2);
	map<pair<Tomography::det_type,int>,vector<TGraph*> > signal_shape;
	map<pair<Tomography::det_type,int>,vector<TF1*> > rising_fit;
	for(vector<Detector*>::const_iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		for(int i=0;i<(*det_it)->get_Nchannel();i++){
			signal_shape[pair<Tomography::det_type,int>((*det_it)->get_type(),(*det_it)->get_n_in_tree())].push_back(new TGraph());
			rising_fit[pair<Tomography::det_type,int>((*det_it)->get_type(),(*det_it)->get_n_in_tree())].push_back(new TF1("rising_fit","pol1(0)",0,Tomography::Nsample));
		}
	}
	long nentries = Min(fChain->GetEntriesFast(),static_cast<Long64_t>(evn_max));
	if(evn_min>nentries) return;
	for(long i=evn_min;i<nentries && Tomography::can_continue;i++){
		LoadTree(i);
		GetEntry(i);
		int det_id = 1;
		for(vector<Detector*>::const_iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
			vector<vector<double> > current_ampl = get_ampl((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			pair<Tomography::det_type,int> current_key((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			for(int j=0;j<(*det_it)->get_Nchannel();j++){
				double max_ampl = -Tomography::ADC_max;
				double sample_max = -1;
				for(int k=0;k<Tomography::Nsample;k++){
					signal_shape[current_key][j]->SetPoint(k,k,current_ampl[j][k]);
					if(current_ampl[j][k]>max_ampl){
						max_ampl = current_ampl[j][k];
						sample_max = k;
					}
				}

				int k = sample_max;
				double mean_xx = 0;
				double mean_xy = 0;
				double mean_x = 0;
				double mean_y = 0;
				while(k>=Tomography::SampleMin && current_ampl[j][k]>(Tomography::sigma*(*det_it)->get_RMS(j))){
					mean_xx += k*k;
					mean_x += k;
					mean_xy += k*current_ampl[j][k];
					mean_y += current_ampl[j][k];
					k--;
				}
				int point_n = sample_max - k;
				mean_xx /= point_n;
				mean_xy /= point_n;
				mean_y /= point_n;
				mean_x /= point_n;
				double slope = (mean_xy - mean_x*mean_y)/(mean_xx - mean_x*mean_x);
				rising_fit[current_key][j]->SetParameters(mean_y - slope*mean_x,slope);
				cDisplay->cd(det_id);
				if(j==0){
					signal_shape[current_key][j]->GetHistogram()->SetMinimum(-100);
					signal_shape[current_key][j]->GetHistogram()->SetMaximum(1200);
					signal_shape[current_key][j]->Draw("AL");
					if(point_n>2) rising_fit[current_key][j]->Draw("same");
				}
				else{
					signal_shape[current_key][j]->Draw("L");
					if(point_n>2) rising_fit[current_key][j]->Draw("same");
				}
			}
			det_id++;
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
	cControl->Divide(2,2);
	TCanvas * cAnalyse = new TCanvas("cAnalyse");
	cAnalyse->Divide(2,2);
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
	TH2D * slope_corr_X = new TH2D("slope_corr_X","slope_corr_X",50,-1,1,50,-1,1);
	TH2D * slope_corr_Y = new TH2D("slope_corr_Y","slope_corr_Y",50,-1,1,50,-1,1);

	TH1D * clus_pos = new TH1D("clus_pos","clus_pos",1024,0,Tomography::XY_size);
	TH1D * clus_size = new TH1D("clus_size","clus_size",50,-1,47);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	TH1D * ray_phi = new TH1D("phi","phi",100,-Pi(),Pi());
	TH1D * ray_slope_X = new TH1D("slope_X","slope_X",100,0,1);
	TH1D * ray_slope_Y = new TH1D("slope_Y","slope_Y",100,0,1);
	ray_slope->SetLineColor(1);
	ray_slope_X->SetLineColor(2);
	ray_slope_Y->SetLineColor(3);
	double min_angle = -1;
	double max_angle = 10;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long ientry=0;ientry<nentries && Tomography::can_continue;ientry++){
		LoadTree(ientry);
		GetEntry(ientry);

		vector<Event*> events;
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			events.push_back((*it)->build_event(get_ampl((*it)->get_type(),(*it)->get_n_in_tree()),Nevent));
			events.back()->MultiCluster();
			events.back()->do_cuts();
		}
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,events);
		vector<Ray> current_rays = CBEvent->get_absorption_rays();
		for(vector<Ray>::iterator ray_it = current_rays.begin();ray_it!=current_rays.end();++ray_it){
			double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
			if(slope<min_angle || slope>max_angle) continue;
			ray_slope->Fill(ATan(slope));
			ray_slope_X->Fill(ATan(Abs(ray_it->get_slope_X())));
			ray_slope_Y->Fill(ATan(Abs(ray_it->get_slope_Y())));
			double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
			if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
			ray_phi->Fill(phi);
			vector<Cluster*> current_clusters = ray_it->get_clus();
			vector<double> clus_slopes_X;
			vector<double> clus_slopes_Y;
			for(vector<Cluster*>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				clus_pos->Fill((*clus_it)->get_pos()*(Tomography::Static_Detector[(*clus_it)->get_type()]->get_StripPitch()));
				clus_size->Fill((*clus_it)->get_size());
				int left = (((*clus_it)->get_pos() - (*clus_it)->get_size()) > 0) ? ((*clus_it)->get_pos() - (*clus_it)->get_size()) : 0;
				int max_strip = detectors[(*clus_it)->find_det(detectors)]->get_Nstrip();
				int right = (((*clus_it)->get_pos() + (*clus_it)->get_size()) < max_strip) ? ((*clus_it)->get_pos() + (*clus_it)->get_size()) : max_strip;
				int current_n = (*clus_it)->get_n_in_tree();
				int current_X = (*clus_it)->get_is_X();
				TGraphErrors * gFit = new TGraphErrors();
				int gFit_nb = 0;
				TF1 * lFit = new TF1("lfit","pol1(0)",5,15);
				lFit->SetParameters(0,0);
				lFit->SetParLimits(0,-10,10);
				lFit->SetParLimits(1,-2,2);
				vector<vector<double> > det_ampl = get_ampl((*clus_it)->get_type(),current_n);
				for(int j=0;j<Tomography::Nsample;j++){
					double mean = 0;
					double ampl = 0;
					double max_ampl = 0;
					for(int i=left;i<=right;i++){
						double current_ampl = det_ampl[Tomography::Static_Detector[(*clus_it)->get_type()]->StripToChannel(i)][j];
						double current_pos = i - (*clus_it)->get_pos();
						mean = (mean*ampl + current_pos*current_ampl)/(ampl+current_ampl);
						ampl += current_ampl;
						if(current_ampl > max_ampl) max_ampl = current_ampl;
						if(current_X){
							if(ray_it->get_slope_X() > 0) signalShape_X_pos->Fill(current_pos,j,current_ampl);
							else signalShape_X_neg->Fill(current_pos,j,current_ampl);
						}
						else{
							if(ray_it->get_slope_Y() > 0) signalShape_Y_pos->Fill(current_pos,j,current_ampl);
							else signalShape_Y_neg->Fill(current_pos,j,current_ampl);
						}
					}
					if(Abs(mean)<=(*clus_it)->get_size() && j>4 && j<16){
						gFit->SetPoint(gFit_nb,j,mean);
						double error_y = (ampl>0) ? 512./ampl : (right - left);
						gFit->SetPointError(gFit_nb,0.5,error_y);
						gFit_nb++;
					}
				}
				if(gFit_nb>1){
					gFit->Fit(lFit,"QNRFM");
					if(current_X){
						clus_slopes_X.push_back(lFit->GetParameter(1));
						angle_shape_corr_X->Fill(ray_it->get_slope_X(), lFit->GetParameter(1));
					}
					else{
						clus_slopes_Y.push_back(lFit->GetParameter(1));
						angle_shape_corr_Y->Fill(ray_it->get_slope_Y(), lFit->GetParameter(1));
					}
				}
				delete gFit; delete lFit;
				delete *clus_it;
			}
			for(vector<double>::iterator slope_it = clus_slopes_X.begin();slope_it!=clus_slopes_X.end();++slope_it){
				for(vector<double>::iterator slope_jt = slope_it+1;slope_jt!=clus_slopes_X.end();++slope_jt){
					slope_corr_X->Fill(*slope_it,*slope_jt);
				}
			}
			for(vector<double>::iterator slope_it = clus_slopes_Y.begin();slope_it!=clus_slopes_Y.end();++slope_it){
				for(vector<double>::iterator slope_jt = slope_it+1;slope_jt!=clus_slopes_Y.end();++slope_jt){
					slope_corr_Y->Fill(*slope_it,*slope_jt);
				}
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
			ray_slope_X->Draw("SAME");
			ray_slope_Y->Draw("SAME");
			cControl->cd(4);
			ray_phi->Draw();
			cControl->Modified();
			cControl->Update();
			cAnalyse->cd(1);
			angle_shape_corr_X->Draw("COLZ");
			cAnalyse->cd(2);
			angle_shape_corr_Y->Draw("COLZ");
			cAnalyse->cd(3);
			slope_corr_X->Draw("COLZ");
			cAnalyse->cd(4);
			slope_corr_Y->Draw("COLZ");
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
	ray_slope_X->Draw("SAME");
	ray_slope_Y->Draw("SAME");
	cControl->cd(4);
	ray_phi->Draw();
	cControl->Modified();
	cControl->Update();
	cAnalyse->cd(1);
	angle_shape_corr_X->Draw("COLZ");
	cAnalyse->cd(2);
	angle_shape_corr_Y->Draw("COLZ");
	cAnalyse->cd(3);
	slope_corr_X->Draw("COLZ");
	cAnalyse->cd(4);
	slope_corr_Y->Draw("COLZ");
	cAnalyse->Modified();
	cAnalyse->Update();
	cout << endl;
}