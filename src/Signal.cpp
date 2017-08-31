#define Signal_cpp
#include "Signal.h"
#include "event.h"
#include "Tanalyse_W.h"
#include "ray.h"
#include "datareader.h"
#include "point.h"
#include "Tray.h"
#include "cluster.h"

#include "MT_tomography.h"
#include "task/read_signal_task.h"
#include "task/ped_task.h"
#include "task/multicluster_task.h"
#include "task/write_analyse_task.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>
#include <sstream>
#include <limits>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLine.h>
#include <TMath.h>
#include <TF1.h>

//Boost
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

using std::ifstream;
using std::cout;
using std::endl;
using std::flush;
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
}

Signal::~Signal(){
	if(exists){
		delete fIn;
	}
}
void Signal::MultiCluster_raw(){
	cout << "reading pedfile : " << PedName << endl;
	map<Tomography::det_type,vector<vector<double> > > current_ped = CosmicBench::read_pedfile(PedName,det_N);
	cout << "reading RMSfile : " << RMSName << endl;
	map<Tomography::det_type,vector<vector<double> > > current_RMS = CosmicBench::read_pedfile(RMSName,det_N);
	Tanalyse_W * analyseFile = new Tanalyse_W(analyseTree,det_n);
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Nevent",1);
	fChain->SetBranchStatus("evttime",1);
	for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
		if(type_it->second > 0){
			ostringstream name;
			name << "StripAmpl_" << type_it->first;
			fChain->SetBranchStatus(name.str().c_str(),1);
		}
	}

	Output_Task<event_data> * to_write = new Write_Analyse_Task(analyseFile);
	Input_Task * to_do = new Read_Signal_Task<raw_data>(nentries,this, new Ped_Corr_Task(current_ped, new Multicluster_Task(this,to_write)));
	vector<Thread*> threads;
	threads.push_back(new Writer_Thread(to_write));
	(threads.back())->start();
	threads.push_back(new Reader_Thread(to_do));
	(threads.back())->start();
	const unsigned short n_thread = (Tomography::get_instance()->get_thread_number() > threads.size()) ? (Tomography::get_instance()->get_thread_number() - threads.size()) : 1;
	cout << "1 | " << n_thread << " | 1" << endl;
	for(unsigned short i=0;i<n_thread;i++){
		threads.push_back(new Worker_Thread());
		(threads.back())->start();
	}
	Display_Thread::get_instance()->start_count();
	//cout << Tomography::get_instance()->init_count() << "|" << setw(7) << "tasks" << endl;
	bool has_working_thread = true;
	while(has_working_thread && Tomography::get_instance()->get_can_continue()){
		has_working_thread = false;
		for(unsigned short i=0;i<threads.size();i++){
			if(threads[i]->is_working()){
				has_working_thread = true;
				break;
			}
		}
		usleep(10000);
	}
	for(unsigned short i=0;i<threads.size();i++){
		threads[i]->stop();
		delete threads[i];
	}
	Display_Thread::get_instance()->stop_count();

	analyseFile->Write();
	analyseFile->CloseFile();
	//delete analyseFile;
	fChain->SetBranchStatus("*",1);
	delete to_do;
}
void Signal::MultiCluster(){
	cout << "destination file : " << analyseTree << endl;
	Tanalyse_W * analyseFile = new Tanalyse_W(analyseTree,det_n);
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Nevent",1);
	fChain->SetBranchStatus("evttime",1);
	fChain->SetBranchStatus("StripAmpl_*_corr",1);
	for(long i=0;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);
		map<Tomography::det_type,vector<Event*> > events;

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			events[(*it)->get_type()].push_back((*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime));
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
	const int max_hole_nb = 11;
	double max_z = -10000;
	double min_z = 10000;
	int bin_n = 500;
	double min_coord = -10*Tomography::get_instance()->get_XY_size()/10.;
	double max_coord = 10*Tomography::get_instance()->get_XY_size()/10.;
	TH2D * hough_space_X[max_hole_nb];
	TH2D * hough_space_Y[max_hole_nb];
	//vector<TH2D*> hough_space_X(max_hole_nb);
	//vector<TH2D*> hough_space_Y(max_hole_nb);
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * cHough = new TCanvas();
	cHough->Divide(4,max_hole_nb/2);
	TGraph * int_X[max_hole_nb];
	TGraph * int_Y[max_hole_nb];
	//vector<TGraph*> int_X(max_hole_nb);
	//vector<TGraph*> int_Y(max_hole_nb);
	int X_int_nb[max_hole_nb];
	int Y_int_nb[max_hole_nb];
	//vector<int> X_int_nb(max_hole_nb,0);
	//vector<int> Y_int_nb(max_hole_nb,0);
	for(int i = 0; i<max_hole_nb; i++){
		cout << "--- Hole number : " << i << " ---" << endl;
		map<bool,map<int,vector<Cluster*> > > all_cluster;
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			Event * current_event = (*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime);
			current_event->HoughCluster(i);
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
			all_cluster[(*it)->get_is_X()][(*it)->get_layer()].insert(all_cluster[(*it)->get_is_X()][(*it)->get_z()].end(),current_cluster.begin(),current_cluster.end());
			current_event->MultiCluster();
			cout << " (" << current_event->get_NClus() << ")" << endl;
			delete current_event;
			if((*it)->get_z()>max_z) max_z = (*it)->get_z();
			if((*it)->get_z()<min_z) min_z = (*it)->get_z();
		}
		ostringstream hist_name;
		hist_name << "hough_space_" << i << "_";
		hough_space_X[i] = new TH2D((hist_name.str() + "X").c_str(),(hist_name.str() + "X").c_str(),bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
		hough_space_Y[i] = new TH2D((hist_name.str() + "Y").c_str(),(hist_name.str() + "Y").c_str(),bin_n,min_coord,max_coord,bin_n,min_coord,max_coord);
		//max_z+=10;
		//min_z-=10;
		int suitable_clus_n = 0;
		//draw clusters in hough space
		map<bool,map<int,int> > sizes;
		for(map<bool,map<int,vector<Cluster*> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
			for(map<int,vector<Cluster*> >::iterator kt = (jt->second).begin();kt!=(jt->second).end();++kt){
				for(vector<Cluster*>::iterator it=(kt->second).begin();it!=(kt->second).end();++it){
					suitable_clus_n++;
					if(!(Tomography::get_instance()->get_is_up((*it)->get_layer()))){
						for(int j=0;j<bin_n;j++){
							double current_coord_up = min_coord +j*(max_coord-min_coord)/bin_n;
							double current_coord_down = current_coord_up + ((*it)->get_pos_mm() - current_coord_up)*(min_z-max_z)/((*it)->get_z()-max_z);
							if((*it)->get_is_X()){
								hough_space_X[i]->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
							}
							else{
								hough_space_Y[i]->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
							}
						}
					}
					else{
						for(int j=0;j<bin_n;j++){
							double current_coord_down = min_coord +j*(max_coord-min_coord)/bin_n;
							double current_coord_up = current_coord_down + ((*it)->get_pos_mm() - current_coord_down)*(max_z-min_z)/((*it)->get_z()-min_z);
							if((*it)->get_is_X()){
								hough_space_X[i]->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
							}
							else{
								hough_space_Y[i]->Fill(current_coord_down,current_coord_up);// possible wieght : it->get_ampl()
							}
						}
					}
				}
				sizes[jt->first][kt->first] = (kt->second).size();
			}
		}
		int_X[i] = new TGraph();
		int_Y[i] = new TGraph();
		X_int_nb[i] = 0;
		Y_int_nb[i] = 0;
		int_X[i]->SetMarkerStyle(24);
		int_Y[i]->SetMarkerStyle(24);
		int_X[i]->SetMarkerSize(2);
		int_Y[i]->SetMarkerSize(2);
		int_X[i]->SetMarkerColor(2);
		int_Y[i]->SetMarkerColor(2);
		map<bool,Point_2D> hough_ray;
		for(int drop = 0;drop<2;drop++){
			for(map<bool,map<int,vector<Cluster*> > >::iterator jt = all_cluster.begin();jt!=all_cluster.end();++jt){
				if(hough_ray.count(jt->first)>0) continue;
				vector<map<int,int> > comb = CosmicBenchEvent::combinaisons(sizes[jt->first], (drop>0));
				double smallest_distance = numeric_limits<double>::max();
				bool found = false;
				Point_2D best_intersection;
				for(vector<map<int,int> >::iterator kt = comb.begin();kt!=comb.end();++kt){
					vector<Point_2D> intersections;
					for(map<int,int>::iterator map_it = kt->begin();map_it!=kt->end();++map_it){
						Line_2D first_line;
						Cluster * first_cluster = (jt->second)[map_it->first][map_it->second];
						if(!(Tomography::get_instance()->get_is_up(first_cluster->get_layer()))) first_line = Line_2D(Point_2D(min_coord + (first_cluster->get_pos_mm() - min_coord)*(min_z-max_z)/(first_cluster->get_z()-max_z),min_coord),Point_2D(max_coord + (first_cluster->get_pos_mm() - max_coord)*(min_z-max_z)/(first_cluster->get_z()-max_z),max_coord));
						else first_line = Line_2D(Point_2D(min_coord,min_coord + (first_cluster->get_pos_mm() - min_coord)*(max_z-min_z)/(first_cluster->get_z()-min_z)),Point_2D(max_coord,max_coord + (first_cluster->get_pos_mm() - max_coord)*(max_z-min_z)/(first_cluster->get_z()-min_z)));
						map<int,int>::iterator map_jt = map_it;
						map_jt++;
						while(map_jt!=kt->end()){
							Line_2D second_line;
							Cluster * second_cluster = (jt->second)[map_jt->first][map_jt->second];
							if(!(Tomography::get_instance()->get_is_up(second_cluster->get_layer()))) second_line = Line_2D(Point_2D(min_coord + (second_cluster->get_pos_mm() - min_coord)*(min_z-max_z)/(second_cluster->get_z()-max_z),min_coord),Point_2D(max_coord + (second_cluster->get_pos_mm() - max_coord)*(min_z-max_z)/(second_cluster->get_z()-max_z),max_coord));
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
						int_X[i]->SetPoint(X_int_nb[i],best_intersection.get_X(),best_intersection.get_Y());
						X_int_nb[i]++;
					}
					else{
						int_Y[i]->SetPoint(Y_int_nb[i],best_intersection.get_X(),best_intersection.get_Y());
						Y_int_nb[i]++;
					}
				}
			}
		}
		cout << "total number of suitable cluster : " << suitable_clus_n << endl;
		for(map<bool,map<int,vector<Cluster*> > >::iterator coord_it=all_cluster.begin();coord_it!=all_cluster.end();++coord_it){
			for(map<int,vector<Cluster*> >::iterator alt_it=(coord_it->second).begin();alt_it!=(coord_it->second).end();++alt_it){
				for(vector<Cluster*>::iterator clus_it=(alt_it->second).begin();clus_it!=(alt_it->second).end();++clus_it){
					delete *clus_it;
				}
				(alt_it->second).clear();
			}
			(coord_it->second).clear();
		}
		all_cluster.clear();
	}
	//draw "normal" ray in hough space
	vector<Event*> current_events;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		current_events.push_back((*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime));
		(current_events.back())->MultiCluster();
	}
	CosmicBenchEvent * thisEvent = new CosmicBenchEvent(this,current_events);
	vector<Ray> rays = thisEvent->get_absorption_rays();
	delete thisEvent;
	for(vector<Event*>::iterator ev_it=current_events.begin();ev_it!=current_events.end();++ev_it){
		delete *ev_it;
	}
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
	
	for(int i=0;i<max_hole_nb;i++){
		cHough->cd((2*i)+1);
		hough_space_X[i]->Draw("colz");
		if(rays.size()>0) rays_X->Draw("sameP");
		if(X_int_nb[i]>0) int_X[i]->Draw("sameP");
		cHough->cd((2*i)+2);
		hough_space_Y[i]->Draw("colz");
		if(rays.size()>0) rays_Y->Draw("sameP");
		if(Y_int_nb[i]>0) int_Y[i]->Draw("sameP");
	}
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
	for(long i=0;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			Event * current_event = (*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime);
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
	for(long i=0;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			Event * current_event = (*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime);
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
		if(i%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
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

void Signal::DebugHoles(long event_nb){
	TCanvas * cDisplay = new TCanvas("cDisplay","cDisplay",800,600);
	cDisplay->Divide(1,get_det_N_tot());
	long nentries = fChain->GetEntriesFast();
	if(event_nb<0 || event_nb>=nentries){
		cout << "invalid event number" << endl;
		return;
	}
	LoadTree(event_nb);
	GetEntry(event_nb);
	
	for(unsigned int i=0;i<detectors.size();i++){
		Event * current_event = detectors[i]->build_event(get_ampl<double>(detectors[i]->get_type(),detectors[i]->get_n_in_tree()),Nevent,evttime);
		TH1D * current_TOT_hist = current_event->get_TOT_hist();
		current_event->MultiCluster();
		vector<Cluster*> current_cluster = current_event->get_clusters();
		delete current_event;
		for(vector<Cluster*>::iterator clus_it=current_cluster.begin();clus_it!=current_cluster.end();++clus_it){
			cout << (*clus_it)->print() << endl;
			delete *clus_it;
		}
		cDisplay->cd(i+1);
		current_TOT_hist->Draw();
	}
	cDisplay->Modified();
	cDisplay->Update();
}

void Signal::EventDisplay(int evn_min, int evn_max, Tomography::signal_type signal_correction){
	int column_nb = CeilNint(2*Sqrt((get_det_N_tot())/3.));
	TCanvas * cDisplay = new TCanvas("cDisplay","cDisplay",800,600);
	cDisplay->Divide(column_nb,1+(get_det_N_tot()/column_nb));
	map<pair<Tomography::det_type,int>,vector<TGraph*> > signal_shape;
	map<pair<Tomography::det_type,int>,vector<TF1*> > rising_fit;
	for(vector<Detector*>::const_iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		for(int i=0;i<(*det_it)->get_Nchannel();i++){
			signal_shape[pair<Tomography::det_type,int>((*det_it)->get_type(),(*det_it)->get_n_in_tree())].push_back(new TGraph());
			rising_fit[pair<Tomography::det_type,int>((*det_it)->get_type(),(*det_it)->get_n_in_tree())].push_back(new TF1("rising_fit","pol1(0)",0,Tomography::get_instance()->get_Nsample()));
		}
	}
	long nentries = Min(fChain->GetEntriesFast(),static_cast<Long64_t>(evn_max));
	if(evn_min>nentries) return;
	for(long i=evn_min;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);
		int det_id = 1;
		ostringstream canvas_title;
		canvas_title << "cDisplay_Event_" << Nevent;
		if(signal_correction == Tomography::raw) canvas_title << "_raw";
		else if(signal_correction == Tomography::ped) canvas_title << "_ped";
		else if(signal_correction == Tomography::corr) canvas_title << "_corr";
		else canvas_title << "_unknown";
		for(vector<Detector*>::const_iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
			vector<vector<double> > current_ampl;
			if(signal_correction == Tomography::raw) current_ampl = get_ampl_raw<double>((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			else if(signal_correction == Tomography::ped) current_ampl = get_ampl_ped<double>((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			else current_ampl = get_ampl<double>((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			pair<Tomography::det_type,int> current_key((*det_it)->get_type(),(*det_it)->get_n_in_tree());
			for(int j=0;j<(*det_it)->get_Nchannel();j++){
				double max_ampl = -Tomography::ADC_max;
				double sample_max = -1;
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
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
				while(k>=Tomography::get_instance()->get_SampleMin() && current_ampl[j][k]>(Tomography::get_instance()->get_sigma()*(*det_it)->get_RMS(j))){
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
					ostringstream current_title;
					current_title << current_key.first << "_" << current_key.second << "_Event_" << Nevent;
					if(signal_correction == Tomography::raw) current_title << "_raw";
					else if(signal_correction == Tomography::ped) current_title << "_ped";
					else if(signal_correction == Tomography::corr) current_title << "_corr";
					else current_title << "_unknown";
					signal_shape[current_key][j]->SetTitle(current_title.str().c_str());
					if(signal_correction == Tomography::raw){
						signal_shape[current_key][j]->GetHistogram()->SetMinimum(0);
						signal_shape[current_key][j]->GetHistogram()->SetMaximum(4096);
					}
					else{
						signal_shape[current_key][j]->GetHistogram()->SetMinimum(-100);
						signal_shape[current_key][j]->GetHistogram()->SetMaximum(1200);
					}
					signal_shape[current_key][j]->Draw("AL");
					//if(point_n>2) rising_fit[current_key][j]->Draw("same");
				}
				else{
					signal_shape[current_key][j]->Draw("L");
					//if(point_n>2) rising_fit[current_key][j]->Draw("same");
				}
			}
			det_id++;
		}
		cDisplay->SetTitle(canvas_title.str().c_str());
		cDisplay->Modified();
		cDisplay->Update();
		gSystem->Sleep(1000);
	}
}

void Signal::SignalDispersion(){
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * cControl = new TCanvas("cControl");
	cControl->Divide(2,2);
	TCanvas * cAnalyse = new TCanvas("cAnalyse");
	cAnalyse->Divide(2,2);
	int nbin_time = Tomography::get_instance()->get_Nsample();
	map<int,TH2D*> signalShape_X;
	map<int,TH2D*> signalShape_Y;
	map<int,TGraph*> mean_X;
	map<int,TGraph*> mean_Y;
	int divisions = 11;
	double division_limit = 30;
	TCanvas * cDisplay = new TCanvas("cDisplay");
	cDisplay->Divide(divisions,2);
	for(int i=0;i<divisions;i++){
		//double division_mean = 90.*(((2.*i+1.)/divisions)-1.);
		ostringstream oss_title;
		oss_title << "signalShape_" << (division_limit*(2.*i/divisions - 1)) << "_" << (division_limit*(2.*(i+1)/divisions - 1)) << "_";
		ostringstream oss_name;
		oss_name << "signalShape_" << i << "_";
		signalShape_X[i] = new TH2D((oss_name.str() + "X").c_str(),(oss_title.str() + "X").c_str(),60,-6,6,nbin_time,0,nbin_time);
		signalShape_Y[i] = new TH2D((oss_name.str() + "Y").c_str(),(oss_title.str() + "Y").c_str(),60,-6,6,nbin_time,0,nbin_time);
		mean_X[i] = new TGraph();
		mean_Y[i] = new TGraph();
		mean_X[i]->SetMarkerStyle(2);
		mean_Y[i]->SetMarkerStyle(2);
	}
	/*
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
	*/

	TH2D * angle_shape_corr_X = new TH2D("angle_shape_corr_X","angle_shape_corr_X",50,-0.5,0.5,50,-1,1);
	TH2D * angle_shape_corr_Y = new TH2D("angle_shape_corr_Y","angle_shape_corr_Y",50,-0.5,0.5,50,-1,1);
	TH2D * slope_corr_X = new TH2D("slope_corr_X","slope_corr_X",50,-1,1,50,-1,1);
	TH2D * slope_corr_Y = new TH2D("slope_corr_Y","slope_corr_Y",50,-1,1,50,-1,1);

	TH1D * clus_pos = new TH1D("clus_pos","clus_pos",1024,0,Tomography::get_instance()->get_XY_size());
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
	for(long ientry=0;ientry<nentries && Tomography::get_instance()->get_can_continue();ientry++){
		LoadTree(ientry);
		GetEntry(ientry);

		vector<Event*> events;
		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			events.push_back((*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime));
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
				vector<vector<double> > det_ampl = get_ampl<double>((*clus_it)->get_type(),current_n);
				for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
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
							for(int k=0;k<divisions;k++){
								double min_slope = (division_limit*Pi()/180)*((k*2./divisions)-1);
								double max_slope = (division_limit*Pi()/180)*((2.*(k+1)/divisions)-1);
								double current_slope = ATan(ray_it->get_slope_X());
								if(current_slope>min_slope && current_slope<max_slope){
									signalShape_X[k]->Fill(current_pos,j,current_ampl);
									break;
								}
							}
							//if(ray_it->get_slope_X() > 0) signalShape_X_pos->Fill(current_pos,j,current_ampl);
							//else signalShape_X_neg->Fill(current_pos,j,current_ampl);
						}
						else{
							for(int k=0;k<divisions;k++){
								double min_slope = (division_limit*Pi()/180)*((k*2./divisions)-1);
								double max_slope = (division_limit*Pi()/180)*((2.*(k+1)/divisions)-1);
								double current_slope = ATan(ray_it->get_slope_Y());
								if(current_slope>min_slope && current_slope<max_slope){
									signalShape_Y[k]->Fill(current_pos,j,current_ampl);
									break;
								}
							}
							//if(ray_it->get_slope_Y() > 0) signalShape_Y_pos->Fill(current_pos,j,current_ampl);
							//else signalShape_Y_neg->Fill(current_pos,j,current_ampl);
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
		if((ientry%5000 == 0) && Tomography::get_instance()->get_live_graphic_display()){
			for(int j=0;j<divisions;j++){
				for(int i=1;i<=nbin_time;i++){
					mean_X[j]->SetPoint(i-1,signalShape_X[j]->ProjectionX("_px",i,i)->GetMean(),i-0.5);
					mean_Y[j]->SetPoint(i-1,signalShape_Y[j]->ProjectionX("_px",i,i)->GetMean(),i-0.5);
				}
				cDisplay->cd(j+1);
				signalShape_X[j]->Draw("COLZ");
				mean_X[j]->Draw("PSAME");
				cDisplay->cd(divisions+j+1);
				signalShape_Y[j]->Draw("COLZ");
				mean_Y[j]->Draw("PSAME");
			}
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
	for(int j=0;j<divisions;j++){
		for(int i=1;i<=nbin_time;i++){
			mean_X[j]->SetPoint(i-1,signalShape_X[j]->ProjectionX("_px",i,i)->GetMean(),i-0.5);
			mean_Y[j]->SetPoint(i-1,signalShape_Y[j]->ProjectionX("_px",i,i)->GetMean(),i-0.5);
		}
		cDisplay->cd(j+1);
		signalShape_X[j]->Draw("COLZ");
		mean_X[j]->Draw("PSAME");
		cDisplay->cd(divisions+j+1);
		signalShape_Y[j]->Draw("COLZ");
		mean_Y[j]->Draw("PSAME");
	}
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

void Signal::ConvClusterTest(){

	TH1D * pos_diff = new TH1D("pos_diff","pos_diff",500,-250,250);
	TH1D * max_ampl_diff = new TH1D("max_ampl_diff","max_ampl_diff",500,0,1000);
	TH1D * ampl_track_diff = new TH1D("ampl_track_diff","ampl_track_diff",500,0,1000);
	TH1D * pos_track_diff = new TH1D("pos_track_diff","pos_track_diff",500,-250,250);

	TCanvas * cTest = new TCanvas();
	cTest->Divide(2,2);

	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Nevent",1);
	fChain->SetBranchStatus("evttime",1);
	fChain->SetBranchStatus("StripAmpl_*_corr",1);
	for(long i=0;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);
		vector<Event*> events_mult;
		vector<Event*> events_conv;
		map<Tomography::det_type,map<int,vector<Cluster*> > > clusters_mult;
		map<Tomography::det_type,map<int,vector<Cluster*> > > clusters_conv;
		map<Tomography::det_type,map<int,double> > max_ampl_conv;

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			events_mult.push_back((*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime));
			(events_mult.back())->MultiCluster();
			clusters_mult[(*it)->get_type()][(*it)->get_n_in_tree()] = (events_mult.back())->get_clusters();
			events_conv.push_back((*it)->build_event(get_ampl<double>((*it)->get_type(),(*it)->get_n_in_tree()),Nevent,evttime));
			(events_conv.back())->ConvCluster();
			clusters_conv[(*it)->get_type()][(*it)->get_n_in_tree()] = (events_conv.back())->get_clusters();
			double first_ampl = 0;
			double second_ampl = 0;
			for(vector<Cluster*>::iterator clus_it = clusters_conv[(*it)->get_type()][(*it)->get_n_in_tree()].begin();clus_it!=clusters_conv[(*it)->get_type()][(*it)->get_n_in_tree()].end();++clus_it){
				if((*clus_it)->get_ampl() > first_ampl) first_ampl = (*clus_it)->get_ampl();
				else if((*clus_it)->get_ampl() > second_ampl) second_ampl = (*clus_it)->get_ampl();
			}
			max_ampl_diff->Fill(first_ampl - second_ampl);
			max_ampl_conv[(*it)->get_type()][(*it)->get_n_in_tree()] = first_ampl;
			vector<Cluster*> clusters_conv_tmp = clusters_conv[(*it)->get_type()][(*it)->get_n_in_tree()];
			for(vector<Cluster*>::iterator mult_it = clusters_mult[(*it)->get_type()][(*it)->get_n_in_tree()].begin();mult_it!=clusters_mult[(*it)->get_type()][(*it)->get_n_in_tree()].end();++mult_it){
				vector<Cluster*>::iterator best_conv_cluster = clusters_conv_tmp.end();
				double shorter_dist = 600;
				for(vector<Cluster*>::iterator clus_it = clusters_conv_tmp.begin();clus_it!=clusters_conv_tmp.end();++clus_it){
					double current_dist = (*clus_it)->get_pos_mm() - (*mult_it)->get_pos_mm();
					if(Abs(current_dist) < Abs(shorter_dist)){
						shorter_dist = current_dist;
						best_conv_cluster = clus_it;
					}
				}
				if(best_conv_cluster!=clusters_conv_tmp.end()){
					pos_diff->Fill(shorter_dist);
					clusters_conv_tmp.erase(best_conv_cluster);
				}
				delete *mult_it;
			}
		}
		CosmicBenchEvent * CBE_mult = new CosmicBenchEvent(this,events_mult);
		vector<Ray> current_rays = CBE_mult->get_absorption_rays();
		for(vector<Ray>::iterator ray_it = current_rays.begin();ray_it!=current_rays.end();++ray_it){
			vector<Cluster*> current_clusters = ray_it->get_clus();
			for(vector<Cluster*>::iterator clus_it=current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				int det_id = (*clus_it)->find_det(detectors);
				double min_pos = 500;
				double min_ampl = 0;
				for(vector<Cluster*>::iterator conv_it = clusters_conv[detectors[det_id]->get_type()][detectors[det_id]->get_n_in_tree()].begin();conv_it != clusters_conv[detectors[det_id]->get_type()][detectors[det_id]->get_n_in_tree()].end();++conv_it){
					double current_pos = Abs((*conv_it)->get_pos_mm() - (*clus_it)->get_pos_mm());
					if(current_pos < min_pos){
						min_pos = current_pos;
						min_ampl = (*conv_it)->get_ampl();
					}
				}
				pos_track_diff->Fill(min_pos);
				ampl_track_diff->Fill(max_ampl_conv[detectors[det_id]->get_type()][detectors[det_id]->get_n_in_tree()] - min_ampl);
			}
		}
		
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
		if(i%5000 == 0){
			cTest->cd(1);
			pos_diff->Draw();
			cTest->cd(2);
			max_ampl_diff->Draw();
			cTest->cd(3);
			pos_track_diff->Draw();
			cTest->cd(4);
			ampl_track_diff->Draw();
			cTest->Modified();
			cTest->Update();
		}
		for(vector<Event*>::iterator ev_it = events_mult.begin();ev_it!=events_mult.end();++ev_it){
			delete *ev_it;
		}
		for(vector<Event*>::iterator ev_it = events_conv.begin();ev_it!=events_conv.end();++ev_it){
			delete *ev_it;
		}
		delete CBE_mult;
		for(map<Tomography::det_type,map<int,vector<Cluster*> > >::iterator type_it=clusters_conv.begin();type_it!=clusters_conv.end();++type_it){
			for(map<int,vector<Cluster*> >::iterator det_it=(type_it->second).begin();det_it!=(type_it->second).end();++det_it){
				for(vector<Cluster*>::iterator clus_it=(det_it->second).begin();clus_it!=(det_it->second).end();++clus_it){
					delete *clus_it;
				}
			}
		}
	}
	cout << "\r" << nentries << "/" << nentries << endl;
	fChain->SetBranchStatus("*",1);
	cTest->cd(1);
	pos_diff->Draw();
	cTest->cd(2);
	max_ampl_diff->Draw();
	cTest->cd(3);
	pos_track_diff->Draw();
	cTest->cd(4);
	ampl_track_diff->Draw();
	cTest->Modified();
	cTest->Update();
}
void Signal::NoiseLevels(){
	double Ymin=-500;
	double Ymax=500;
	int bin_n = 500;
	int sample_min = 1;
	int sample_max = Tomography::get_instance()->get_Nsample()-1;//Min(Nsample,4);
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	map<Tomography::det_type,vector<vector<TH1F*> > > ampl_hist_raw;
	map<Tomography::det_type,vector<vector<TH1F*> > > ampl_hist_ped;
	map<Tomography::det_type,vector<vector<TH1F*> > > ampl_hist_corr;
	int tot_chan = 0;
	vector<TLine*> sep_raw;
	ostringstream det_order;
	det_order << "|";
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		int j = (*det_it)->get_n_in_tree();
		Tomography::det_type current_type = (*det_it)->get_type();
		tot_chan += ((*det_it)->get_Nchannel());
		ampl_hist_raw[current_type].push_back(vector<TH1F*>(((*det_it)->get_Nchannel()),NULL));
		ampl_hist_ped[current_type].push_back(vector<TH1F*>(((*det_it)->get_Nchannel()),NULL));
		ampl_hist_corr[current_type].push_back(vector<TH1F*>(((*det_it)->get_Nchannel()),NULL));
		for(int i=0;i<((*det_it)->get_Nchannel());i++){
			ostringstream name;
			name << "ampl_hist_" << current_type << j << "_" << i << "_";
			ampl_hist_raw[current_type][j][i] = new TH1F((name.str() + "raw").c_str(),(name.str() + "raw").c_str(),bin_n,Ymin,Ymax);
			ampl_hist_ped[current_type][j][i] = new TH1F((name.str() + "ped").c_str(),(name.str() + "ped").c_str(),bin_n,Ymin,Ymax);
			ampl_hist_corr[current_type][j][i] = new TH1F((name.str() + "corr").c_str(),(name.str() + "corr").c_str(),bin_n,Ymin,Ymax);
		}
		sep_raw.push_back(new TLine(tot_chan,0,tot_chan,4096));
		(sep_raw.back())->SetLineStyle(2);
		(sep_raw.back())->SetLineColor(2);
		det_order << " " << current_type << "_" << j << " |";
	}
	TProfile * ampl_prof_raw = new TProfile("prof_raw","prof_raw",tot_chan,0,tot_chan);
	//ampl_prof_raw->SetTitle(det_order.str().c_str());
	TProfile * ampl_prof_ped = new TProfile("prof_ped","prof_ped",tot_chan,0,tot_chan);
	TProfile * ampl_prof_corr = new TProfile("prof_corr","prof_corr",tot_chan,0,tot_chan);
	TGraph * RMS_raw = new TGraph();
	//RMS_raw->SetTitle(det_order.str().c_str());
	TGraph * RMS_ped = new TGraph();
	TGraph * RMS_corr = new TGraph();
	for(long n=0;n<nentries && Tomography::get_instance()->get_can_continue();n++){
		LoadTree(n);
		GetEntry(n);
		int n_chan = 0;
		for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist_raw.begin();type_it!=ampl_hist_raw.end();++type_it){
			for(unsigned int i=0;i<(type_it->second).size();i++){
				vector<vector<float> > current_ampl_raw = get_ampl_raw<float>(type_it->first,i);
				vector<vector<float> > current_ampl_ped = get_ampl_ped<float>(type_it->first,i);
				vector<vector<float> > current_ampl_corr = get_ampl<float>(type_it->first,i);
				for(unsigned int j=0;j<(type_it->second)[i].size();j++){
					for(int k=sample_min;k<sample_max;k++){
						(type_it->second)[i][j]->Fill(current_ampl_raw[j][k]);
						ampl_prof_raw->Fill(n_chan,current_ampl_raw[j][k]);
						ampl_hist_ped[type_it->first][i][j]->Fill(current_ampl_ped[j][k]);
						ampl_prof_ped->Fill(n_chan,current_ampl_ped[j][k]);
						ampl_hist_corr[type_it->first][i][j]->Fill(current_ampl_corr[j][k]);
						ampl_prof_corr->Fill(n_chan,current_ampl_corr[j][k]);
					}
					n_chan++;
				}
			}
		}
		if((n%100) == 0) cout << "\rcomputing RMS (" << n << "/" << nentries << ")" << flush;
	}
	cout << "\rcomputing RMS (" << nentries << "/" << nentries << ")" << endl;
	int n_chan = 0;
	for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist_raw.begin();type_it!=ampl_hist_raw.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				TFitResultPtr res = (type_it->second)[i][j]->Fit("gaus","SQN");
				if(res->IsEmpty()){
					RMS_raw->SetPoint(n_chan,n_chan,-1);
					cout << "problem fitting raw channel " << type_it->first << "_" << i << "_" << j << endl;
				}
				else{
					RMS_raw->SetPoint(n_chan,n_chan,res->Parameter(2));
				}
				n_chan++;
			}
		}
	}
	n_chan = 0;
	for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist_ped.begin();type_it!=ampl_hist_ped.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				TFitResultPtr res = (type_it->second)[i][j]->Fit("gaus","SQN");
				if(res->IsEmpty()){
					RMS_ped->SetPoint(n_chan,n_chan,-1);
					cout << "problem fitting ped channel " << type_it->first << "_" << i << "_" << j << endl;
				}
				else{
					RMS_ped->SetPoint(n_chan,n_chan,res->Parameter(2));
				}
				n_chan++;
			}
		}
	}
	n_chan = 0;
	for(map<Tomography::det_type,vector<vector<TH1F*> > >::iterator type_it = ampl_hist_corr.begin();type_it!=ampl_hist_corr.end();++type_it){
		for(unsigned int i=0;i<(type_it->second).size();i++){
			for(unsigned int j=0;j<(type_it->second)[i].size();j++){
				TFitResultPtr res = (type_it->second)[i][j]->Fit("gaus","SQN");
				if(res->IsEmpty()){
					RMS_corr->SetPoint(n_chan,n_chan,-1);
					cout << "problem fitting corr channel " << type_it->first << "_" << i << "_" << j << endl;
				}
				else{
					RMS_corr->SetPoint(n_chan,n_chan,res->Parameter(2));
				}
				n_chan++;
			}
		}
	}
	TCanvas * c_prof = new TCanvas();
	c_prof->SetTitle(det_order.str().c_str());
	TCanvas * c_RMS = new TCanvas();
	c_RMS->SetTitle(det_order.str().c_str());
	c_prof->Divide(1,3);
	c_RMS->Divide(1,3);
	
	c_prof->cd(1);
	ampl_prof_raw->Draw();
	double current_min = ampl_prof_raw->GetMinimum();
	double current_max = ampl_prof_raw->GetMaximum();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}
	c_RMS->cd(1);
	RMS_raw->Draw();
	current_min = RMS_raw->GetYaxis()->GetXmin();
	current_max = RMS_raw->GetYaxis()->GetXmax();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}
	c_prof->cd(2);
	ampl_prof_ped->Draw();
	current_min = ampl_prof_ped->GetMinimum();
	current_max = ampl_prof_ped->GetMaximum();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}
	c_RMS->cd(2);
	RMS_ped->Draw();
	current_min = RMS_ped->GetYaxis()->GetXmin();
	current_max = RMS_ped->GetYaxis()->GetXmax();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}
	c_prof->cd(3);
	ampl_prof_corr->Draw();
	current_min = ampl_prof_corr->GetMinimum();
	current_max = ampl_prof_corr->GetMaximum();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}
	c_RMS->cd(3);
	RMS_corr->Draw();
	current_min = RMS_corr->GetYaxis()->GetXmin();
	current_max = RMS_corr->GetYaxis()->GetXmax();
	for(vector<TLine*>::iterator line_it=sep_raw.begin();line_it!=sep_raw.end();++line_it){
		(*line_it)->SetY1(current_min);
		(*line_it)->SetY2(current_max);
		(*line_it)->DrawClone();
	}

}
