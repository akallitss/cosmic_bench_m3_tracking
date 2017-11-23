#define analyse_cpp
//my class
#include "analyse.h"
#include "ray.h"
#include "event.h"
#include "Tsignal_R.h"
#include "tomography.h"
#include "acceptanceFunction.h"
#include "Tray.h"
#include "cluster.h"
//ROOT
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TLine.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLine.h>
//Boost
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
//std
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>

//std
using std::cout;
using std::endl;
using std::flush;
using std::setw;
using std::vector;
using std::ostringstream;
using std::numeric_limits;
using std::max_element;
using std::left;
using std::right;
using std::pair;
//tmath
using TMath::Sqrt;
using TMath::ATan;
using TMath::MaxElement;
using TMath::Abs;
using TMath::FloorNint;
using TMath::CeilNint;
using TMath::Max;
using TMath::Min;
using TMath::Pi;
using TMath::Sin;
using TMath::Cos;
using TMath::Power;

Analyse::Analyse(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CosmicBench::Init(config_tree);
	Tanalyse_R::Init(tree,det_n);
	signal_file_name = config_tree.get<string>("signal_file");
}
Analyse::Analyse(ptree config_tree){
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CosmicBench::Init(config_tree);
	Tanalyse_R::Init(tree,det_n);
	signal_file_name = config_tree.get<string>("signal_file");
}
Analyse::~Analyse(){
	//delete f;
}
void Analyse::Residus(){
	map<Tomography::det_type,vector<TCanvas*> > c_MM;
	map<Tomography::det_type,vector<TH1D*> > hist_residus;
	int nbins = 200;
	long eventReconstructed = 0;
	long eventSuitable = 0;
	for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
		if(type_it->second > 0){
			c_MM[type_it->first] = vector<TCanvas*>(type_it->second,NULL);
			hist_residus[type_it->first] = vector<TH1D*>(type_it->second,NULL);
		}
		for(int i=0;i<type_it->second;i++){
			ostringstream name;
			name << type_it->first << "_" << i;
			c_MM[type_it->first][i] = new TCanvas(name.str().c_str(),name.str().c_str());
			name << "_residu";
			hist_residus[type_it->first][i] = new TH1D(name.str().c_str(),name.str().c_str(),nbins,-Tomography::get_instance()->get_XY_size(),Tomography::get_instance()->get_XY_size());
		}
	}
	TCanvas * c0 = new TCanvas("chiSquare","chiSquare");
	TH1D * chiSquareH = new TH1D("chiSquares","chiSquares",nbins,0,3*Tomography::get_instance()->get_XY_size());

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		//eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
				chiSquareH->Fill(it->get_chiSquare_X()+it->get_chiSquare_Y());
				eventReconstructed++;
			}
			else{
				continue;
			}
			for(vector<Detector*>::iterator jt=detectors.begin();jt!=detectors.end();++jt){
				double residu = it->get_residu(*jt);
				if(residu != numeric_limits<double>::min()){
					hist_residus[(*jt)->get_type()][(*jt)->get_n_in_tree()]->Fill(residu);
				}
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
				for(int i=0;i<type_it->second;i++){
					c_MM[type_it->first][i]->cd();
					hist_residus[type_it->first][i]->Draw();
					c_MM[type_it->first][i]->Modified();
					c_MM[type_it->first][i]->Update();
				}
			}
			c0->cd();
			chiSquareH->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
		for(int i=0;i<type_it->second;i++){
			c_MM[type_it->first][i]->cd();
			hist_residus[type_it->first][i]->Draw();
			c_MM[type_it->first][i]->Modified();
			c_MM[type_it->first][i]->Update();
		}
	}
	c0->cd();
	chiSquareH->Draw();
	c0->Modified();
	c0->Update();
}
void Analyse::Residus_time(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	LoadTree(0);
	fChain->GetEntry(0);
	double evttime_min = evttime;
	LoadTree(nentries-1);
	fChain->GetEntry(nentries-1);
	double evttime_max = evttime;


	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	map<string,TCanvas*> c_MM;
	map<string,TH2D*> muon_seen;
	map<string,TH2D*> muon_total;
	map<string,TH2D*> efficacity_2D;
	map<string,TH2D*> ampl_h;
	map<string,TGraph*> correlation;
	map<string,int> point_nb;
	map<string,double> efficacity;
	map<string,TProfile*> efficiency_time;
	map<string,TProfile*> amplitude_time;
	int nbins = 200;
	double marge = 1./10.;
	int nbins_2D = 100*(1+2*marge);
	long eventReconstructed = 0;
	double eventSuitable = 0;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
			c_MM[name.str()]->Divide(3,2);
			muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			ampl_h[name.str()] = new TH2D((name.str()+"_ampl_mean").c_str(),(name.str()+"_ampl_mean").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			efficacity_2D[name.str()]->SetStats(false);
			correlation[name.str()] = new TGraph();
			point_nb[name.str()] = 0;
			efficacity[name.str()] = 0;
			efficiency_time[name.str()] = new TProfile((name.str()+"_efficiency_time").c_str(),(name.str()+"_efficiency_time").c_str(),nentries/50,evttime_min,evttime_max);
			amplitude_time[name.str()] = new TProfile((name.str()+"_amplitude_time").c_str(),(name.str()+"_amplitude_time").c_str(),nentries/50,evttime_min,evttime_max);
			if((*it)->get_is_X()) nref_x_n++;
		}
		//if((*it)->get_is_X()) det_x_n++;
	}
	/*
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	*/
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
		}
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(3,2);
	TProfile * chisquares = new TProfile("chiSquares","chiSquares",nentries/50,evttime_min,evttime_max);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",get_det_N_tot() + 2,0,get_det_N_tot() + 2);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	TH1D * ray_phi = new TH1D("phi","phi",100,-Pi(),Pi());
	TH1D * ray_slope_X = new TH1D("slope_X","slope_X",100,0,1);
	TH1D * ray_slope_Y = new TH1D("slope_Y","slope_Y",100,0,1);
	TProfile * freq_time = new TProfile("freq_time","freq_time",nentries/50,evttime_min,evttime_max);
	ray_slope->SetLineColor(1);
	ray_slope_X->SetLineColor(2);
	ray_slope_Y->SetLineColor(3);
	if (fChain == 0) return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	double evttime_last = evttime_min;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				chisquares->Fill(evttime,ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y());
				ray_clus_n->Fill(ray_it->get_clus_n());
				double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
				ray_slope->Fill(ATan(slope));
				ray_slope_X->Fill(ATan(Abs(ray_it->get_slope_X())));
				ray_slope_Y->Fill(ATan(Abs(ray_it->get_slope_Y())));
				double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
				if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
				ray_phi->Fill(phi);
				++ray_it;
			}
			else ray_it = currentRays.erase(ray_it);
		}
		freq_time->Fill(evttime,evttime - evttime_last);
		evttime_last = evttime;
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(get_det_N_tot());
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				ostringstream name;
				name <<(*it)->get_type() << "_" << (*it)->get_n_in_tree();
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				double biggest_ampl = 0;
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					double current_ampl = (*kt)->get_maxStripAmpl();
					if(current_ampl>biggest_ampl) biggest_ampl = current_ampl;
				}
				amplitude_time[name.str()]->Fill(evttime,biggest_ampl);
				for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
					//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
					//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
					//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
					//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
					if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
					if(jt->get_clus_n()<static_cast<unsigned int>(get_layers_n()-non_ref_n)) continue;
					double residu = numeric_limits<double>::max();
					vector<Cluster*>::iterator matching_cluster = current_clusters.end();
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						(*kt)->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(*kt);
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
					if(matching_cluster == current_clusters.end()) continue;
					(*matching_cluster)->set_perp_pos_mm(*jt);
					if((*matching_cluster)->get_is_X()){
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),(*matching_cluster)->get_pos_mm());
					}
					else{
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_pos_mm());
					}
					point_nb[name.str()]++;
					if(residu<chisquare_threshold){
						muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						ampl_h[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_ampl());
					}
					double pos_X = jt->eval_X((*it)->get_z());
					double pos_Y = jt->eval_Y((*it)->get_z());
					if(pos_X<=2*Tomography::get_instance()->get_XY_size()/5. && pos_X>=-2*Tomography::get_instance()->get_XY_size()/5. && pos_Y<=2*Tomography::get_instance()->get_XY_size()/5. && pos_Y>=-2*Tomography::get_instance()->get_XY_size()/5.){
						efficiency_time[name.str()]->Fill(evttime,residu<chisquare_threshold);
					}
					delete *matching_cluster;
					current_clusters.erase(matching_cluster);
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
				for(int i=1;i<=nbins_2D;i++){
					for(int j=1;j<=nbins_2D;j++){
						int binN = muon_total[it->first]->GetBin(i,j);
						double binContent = 0;
						if(muon_total[it->first]->GetBinContent(binN) > 0) binContent = (muon_seen[it->first]->GetBinContent(binN))/(muon_total[it->first]->GetBinContent(binN));
						efficacity_2D[it->first]->SetBinContent(binN,binContent);
					}
				}
				it->second->cd(1);
				efficacity_2D[it->first]->Draw("COLZ");
				it->second->cd(2);
				muon_seen[it->first]->Draw("COLZ");
				it->second->cd(3);
				amplitude_time[it->first]->Draw();
				it->second->cd(4);
				efficiency_time[it->first]->Draw();
				it->second->cd(5);
				if(point_nb[it->first]>0) correlation[it->first]->Draw("AP");
				it->second->cd(6);
				ampl_h[it->first]->Draw("COLZ");
				it->second->Modified();
				it->second->Update();
			}
			c0->cd(1);
			ray_clus_n->Draw();
			c0->cd(2);
			chisquares->Draw();
			c0->cd(3);
			ray_slope->Draw();
			ray_slope_X->Draw("SAME");
			ray_slope_Y->Draw("SAME");
			c0->cd(4);
			ray_phi->Draw();
			c0->cd(5);
			freq_time->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << nentries << endl;
	for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
		double total_seen = 0;
		double total_passed = 0;
		for(int i=1;i<=nbins_2D;i++){
			for(int j=1;j<=nbins_2D;j++){
				int binN = muon_total[it->first]->GetBin(i,j);
				double binContent = 0;
				if(muon_total[it->first]->GetBinContent(binN) > 0) binContent = (muon_seen[it->first]->GetBinContent(binN))/(muon_total[it->first]->GetBinContent(binN));
				efficacity_2D[it->first]->SetBinContent(binN,binContent);
				if(muon_seen[it->first]->GetBinContent(binN) > 0) ampl_h[it->first]->SetBinContent(binN,(ampl_h[it->first]->GetBinContent(binN))/(muon_seen[it->first]->GetBinContent(binN)));
				double pos_X = muon_total[it->first]->GetXaxis()->GetBinCenter(i);
				double pos_Y = muon_total[it->first]->GetYaxis()->GetBinCenter(j);
				if(pos_X<=2*Tomography::get_instance()->get_XY_size()/5. && pos_X>=-2*Tomography::get_instance()->get_XY_size()/5. && pos_Y<=2*Tomography::get_instance()->get_XY_size()/5. && pos_Y>=-2*Tomography::get_instance()->get_XY_size()/5.){
					total_seen += muon_seen[it->first]->GetBinContent(binN);
					total_passed += muon_total[it->first]->GetBinContent(binN);
				}
			}
		}
		efficacity[it->first] = total_seen/total_passed;
		it->second->cd(1);
		efficacity_2D[it->first]->Draw("COLZ");
		it->second->cd(2);
		muon_seen[it->first]->Draw("COLZ");
		it->second->cd(3);
		amplitude_time[it->first]->Draw();
		it->second->cd(4);
		efficiency_time[it->first]->Draw();
		it->second->cd(5);
		if(point_nb[it->first]>0) correlation[it->first]->Draw("AP");
		it->second->cd(6);
		ampl_h[it->first]->Draw("COLZ");
		it->second->Modified();
		it->second->Update();
		cout << it->first << " efficacity : " << 100.*efficacity[it->first] << "%" << endl;
	}
	c0->cd(1);
	ray_clus_n->Draw();
	c0->cd(2);
	chisquares->Draw();
	c0->cd(3);
	ray_slope->Draw();
	ray_slope_X->Draw("SAME");
	ray_slope_Y->Draw("SAME");
	c0->cd(4);
	ray_phi->Draw();
	c0->cd(5);
	freq_time->Draw();
	c0->Modified();
	c0->Update();
}
void Analyse::Amplitude_time(){
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	LoadTree(0);
	fChain->GetEntry(0);
	double evttime_min = evttime;
	LoadTree(nentries-1);
	fChain->GetEntry(nentries-1);
	double evttime_max = evttime;


	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * c_MM = new TCanvas();
	int column_nb = CeilNint(2*Sqrt((3*get_det_N_tot())/3.));
	c_MM->Divide(column_nb,1+(3*get_det_N_tot()/column_nb));
	map<Tomography::det_type,map<int,TProfile*> > amplitude_time;
	//map<string,TProfile*> amplitude_time;
	map<Tomography::det_type,map<int,TH1D*> > ampl_h;
	//map<string,TH1D*> ampl_h;
	map<Tomography::det_type,map<int,TH2D*> > ampl_time_h;
	//map<string,TH2D*> ampl_time_h;
	int n_bins_adc = 800;
	int n_bins_time = nentries/100;
	double clock_rate = Tomography::get_instance()->get_clock_rate();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		ostringstream name;
		name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
		amplitude_time[(*it)->get_type()][(*it)->get_n_in_tree()] = new TProfile((name.str()+"_amplitude_time").c_str(),(name.str()+"_amplitude_time").c_str(),n_bins_time,evttime_min/clock_rate,evttime_max/clock_rate);
		ampl_h[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH1D((name.str()+"_ampl_h").c_str(),(name.str()+"_ampl_h").c_str(),n_bins_adc,0,4096);
		ampl_time_h[(*it)->get_type()][(*it)->get_n_in_tree()] = new TH2D((name.str()+"_ampl_time_h").c_str(),(name.str()+"_ampl_time_h").c_str(),n_bins_adc/2,0,4096,n_bins_time,evttime_min/clock_rate,evttime_max/clock_rate);
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(3);
	int n_bin_log = 200;
	double * bin_edges = new double[n_bin_log+1];
	//double log_max = 20*(evttime_max-evttime_min)/nentries;
	double log_max = clock_rate/(1e5);
	//double log_min = 1e4;
	double log_min = 0.1;//100*clock_rate/((evttime_max-evttime_min));
	for(int i=0;i<=n_bin_log;i++){
		bin_edges[i] = log_min*Power(log_max/log_min,i/static_cast<double>(n_bin_log));
	}

	TH1D * freq_time = new TH1D("freq_time","freq_time",n_bins_time,evttime_min/clock_rate,evttime_max/clock_rate);
	freq_time->Sumw2();
	TH1D * freq_h = new TH1D("freq_h","freq_h",n_bin_log,bin_edges);
	TH2D * freq_time_h = new TH2D("freq_time_h","freq_time_h",n_bin_log,bin_edges,n_bins_time,evttime_min/clock_rate,evttime_max/clock_rate);
	c0->GetPad(2)->SetLogx();
	c0->GetPad(2)->SetLogy();
	c0->GetPad(3)->SetLogx();
	c0->GetPad(3)->SetLogz();
	delete bin_edges;
	if (fChain == 0) return;
	cout << setw(20) << "total processed" << endl;
	double evttime_last = evttime_min;
	for (Long64_t jentry=1; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		freq_time->Fill(evttime/clock_rate);//clock_rate/(evttime - evttime_last));
		freq_h->Fill(clock_rate/(evttime - evttime_last));
		freq_time_h->Fill(clock_rate/(evttime - evttime_last),evttime/clock_rate);
		evttime_last = evttime;
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			//ostringstream name;
			//name <<(*it)->get_type() << "_" << (*it)->get_n_in_tree();
			Tomography::det_type current_type = (*it)->get_type();
			int current_n = (*it)->get_n_in_tree();
			vector<Cluster*> current_clusters = (*it)->get_clusters();
			for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
				double current_ampl = (*kt)->get_maxStripAmpl();
				ampl_h[current_type][current_n]->Fill(current_ampl);
				amplitude_time[current_type][current_n]->Fill(evttime/clock_rate,current_ampl);
				ampl_time_h[current_type][current_n]->Fill(current_ampl,evttime/clock_rate);
				delete *kt;
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r" << setw(20) << jentry << flush;
		if(jentry%20000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			int pad_id= 1;
			for(map<Tomography::det_type,map<int,TProfile*> >::iterator it = amplitude_time.begin();it!=amplitude_time.end();++it){
				for(map<int,TProfile*>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
					c_MM->cd(pad_id);
					jt->second->Draw();
					pad_id++;
					c_MM->cd(pad_id);
					ampl_h[it->first][jt->first]->Draw();
					pad_id++;
					c_MM->cd(pad_id);
					ampl_time_h[it->first][jt->first]->Draw("COLZ");
					pad_id++;
				}
			}
			c_MM->Modified();
			c_MM->Update();
			c0->cd(1);
			freq_time->Draw();
			c0->cd(2);
			freq_h->Draw();
			c0->cd(3);
			freq_time_h->Draw("COLZ");
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r" << setw(20) << nentries << endl;
	freq_time->Scale(n_bins_time*clock_rate/(evttime_max-evttime_min));
	int pad_id= 1;
	for(map<Tomography::det_type,map<int,TProfile*> >::iterator it = amplitude_time.begin();it!=amplitude_time.end();++it){
		for(map<int,TProfile*>::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			c_MM->cd(pad_id);
			jt->second->Draw();
			pad_id++;
			c_MM->cd(pad_id);
			ampl_h[it->first][jt->first]->Draw();
			pad_id++;
			c_MM->cd(pad_id);
			ampl_time_h[it->first][jt->first]->Draw("COLZ");
			pad_id++;
		}
	}
	c_MM->Modified();
	c_MM->Update();
	c0->cd(1);
	freq_time->Draw("E");
	c0->cd(2);
	freq_h->Draw();
	c0->cd(3);
	freq_time_h->Draw("COLZ");
	c0->Modified();
	c0->Update();
}
void Analyse::Residus_ref(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	map<string,TCanvas*> c_MM;
	map<string,TH1D*> MM_residus;
	map<string,TF1*> offset_fit;
	map<string,TH2D*> muon_seen;
	map<string,TH2D*> muon_total;
	map<string,TH2D*> efficacity_2D;
	map<string,TH2D*> ampl_h;
	map<string,TGraph*> correlation;
	map<string,TProfile*> angle_alignment;
	map<string,TF1*> angle_z_fit;
	map<string,TProfile*> resVSpos;
	map<string,TProfile*> resVSampl;
	map<string,TProfile*> resVStime;
	map<string,TProfile*> resVSangle;
	map<string,TF1*> angle_xy_fit;
	map<string,TProfile*> resVSanglePerp;
	map<string,TProfile*> resVStot;
	map<string,TProfile*> resVSsize;
	map<string,TProfile*> absResVStot;
	map<string,TProfile*> absResVSsize;
	map<string,TProfile*> absResVSampl;
	map<string,TProfile*> absResVStime;
	map<string,TProfile*> absResVSabsAngle;
	map<string,int> point_nb;
	map<string,double> efficacity;
	int nbins = 200;
	double marge = 1./10.;
	int nbins_2D = 100*(1+2*marge);
	long eventReconstructed = 0;
	double eventSuitable = 0;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
			c_MM[name.str()]->Divide(4,4);
			MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
			muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			ampl_h[name.str()] = new TH2D((name.str()+"_ampl_mean").c_str(),(name.str()+"_ampl_mean").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			efficacity_2D[name.str()]->SetStats(false);
			correlation[name.str()] = new TGraph();
			angle_alignment[name.str()] = new TProfile((name.str()+"_resVSperpPos").c_str(),(name.str()+"_resVSperpPos").c_str(),500,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,-5,5);
			resVSpos[name.str()] = new TProfile((name.str()+"_resVSpos").c_str(),(name.str()+"_resVSpos").c_str(),500,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,-5,5);
			resVSampl[name.str()] = new TProfile((name.str()+"_resVSampl").c_str(),(name.str()+"_resVSampl").c_str(),500,-100,10000,-5,5);
			resVStime[name.str()] = new TProfile((name.str()+"_resVStime").c_str(),(name.str()+"_resVStime").c_str(),38,-2,34,-5,5);
			resVSangle[name.str()] = new TProfile((name.str()+"_resVSangle").c_str(),(name.str()+"_resVSangle").c_str(),50,-0.6,0.6,-5,5);
			resVSanglePerp[name.str()] = new TProfile((name.str()+"_resVSanglePerp").c_str(),(name.str()+"_resVSanglePerp").c_str(),50,-0.6,0.6,-5,5);
			resVStot[name.str()] = new TProfile((name.str()+"_resVStot").c_str(),(name.str()+"_resVStot").c_str(),26,0,25,-5,5);
			resVSsize[name.str()] = new TProfile((name.str()+"_resVSsize").c_str(),(name.str()+"_resVSsize").c_str(),50,0,50,-5,5);
			absResVSampl[name.str()] = new TProfile((name.str()+"_absResVSampl").c_str(),(name.str()+"_absResVSampl").c_str(),500,-100,10000,0,5);
			absResVStime[name.str()] = new TProfile((name.str()+"_absResVStime").c_str(),(name.str()+"_absResVStime").c_str(),38,-2,34,0,5);
			absResVSabsAngle[name.str()] = new TProfile((name.str()+"_absResVSabsAngle").c_str(),(name.str()+"_absResVSabsAngle").c_str(),50,0,0.6,0,5);
			absResVStot[name.str()] = new TProfile((name.str()+"_absResVStot").c_str(),(name.str()+"_absResVStot").c_str(),26,0,25,0,5);
			absResVSsize[name.str()] = new TProfile((name.str()+"_absResVSsize").c_str(),(name.str()+"_absResVSsize").c_str(),50,0,50,0,5);
			point_nb[name.str()] = 0;
			efficacity[name.str()] = 0;
			offset_fit[name.str()] = new TF1("offset_fit","[3]*exp(-(x-[0])*(x-[0])/(2*[1]*[1])) + [4]*exp(-(x-[0])*(x-[0])/(2*[2]*[2]))",-5,5);
			offset_fit[name.str()]->SetParameters(0,0.5,2,nentries/10);
			offset_fit[name.str()]->SetParLimits(0,-10,10);
			offset_fit[name.str()]->SetParLimits(1,0,1);
			offset_fit[name.str()]->SetParLimits(2,0,10);
			offset_fit[name.str()]->SetParLimits(3,1,nentries);
			offset_fit[name.str()]->SetParLimits(4,1,nentries);
			angle_z_fit[name.str()] = new TF1("angle_z_fit","pol1(0)",-150,150);
			angle_z_fit[name.str()]->SetParameters(0,0);
			angle_z_fit[name.str()]->SetParLimits(0,-5,5);
			angle_z_fit[name.str()]->SetParLimits(1,-1,1);
			angle_xy_fit[name.str()] = new TF1("angle_xy_fit","pol2(0)",-0.3,0.3);
			angle_xy_fit[name.str()]->SetParameters(0,0,0);
			angle_xy_fit[name.str()]->SetParLimits(0,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(1,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(2,-20,20);
			//nref_is_X[name.str()] = (*it)->get_is_X();
			if((*it)->get_is_X()) nref_x_n++;
		}
		//if((*it)->get_is_X()) det_x_n++;
	}
	/*
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	*/
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
		}
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2,2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",get_det_N_tot() + 2,0,get_det_N_tot() + 2);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	TH1D * ray_phi = new TH1D("phi","phi",100,-Pi(),Pi());
	TH1D * ray_slope_X = new TH1D("slope_X","slope_X",100,0,1);
	TH1D * ray_slope_Y = new TH1D("slope_Y","slope_Y",100,0,1);
	ray_slope->SetLineColor(1);
	ray_slope_X->SetLineColor(2);
	ray_slope_Y->SetLineColor(3);
	if (fChain == 0) return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				chisquares->Fill(ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y());
				ray_clus_n->Fill(ray_it->get_clus_n());
				double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
				ray_slope->Fill(ATan(slope));
				ray_slope_X->Fill(ATan(Abs(ray_it->get_slope_X())));
				ray_slope_Y->Fill(ATan(Abs(ray_it->get_slope_Y())));
				double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
				if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
				ray_phi->Fill(phi);
				++ray_it;
			}
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(get_det_N_tot());
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				ostringstream name;
				name <<(*it)->get_type() << "_" << (*it)->get_n_in_tree();
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
					//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
					//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
					//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
					//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
					if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
					//if(jt->get_clus_n()<static_cast<unsigned int>(get_layers_n()-non_ref_n)) continue;
					if(Ray_2D(*jt,(*it)->get_is_X() ? 'X' : 'Y').has_layer((*it)->get_det()->get_layer())) continue;
					if(!is_compatible(*it,&(*jt))) continue;
					double residu = numeric_limits<double>::max();
					vector<Cluster*>::iterator matching_cluster = current_clusters.end();
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						(*kt)->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(*kt);
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
					if(matching_cluster == current_clusters.end()) continue;
					(*matching_cluster)->set_perp_pos_mm(*jt);
					if((*matching_cluster)->get_is_X()){
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),(*matching_cluster)->get_pos_mm());
						angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						resVSpos[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
						absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_X()),Abs(residu));
						resVSanglePerp[name.str()]->Fill(jt->get_slope_Y(),residu);
					}
					else{
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_pos_mm());
						angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						resVSpos[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
						absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_Y()),Abs(residu));
						resVSanglePerp[name.str()]->Fill(jt->get_slope_X(),residu);
					}
					point_nb[name.str()]++;
					MM_residus[name.str()]->Fill(residu);
					resVStime[name.str()]->Fill((*matching_cluster)->get_t(),residu);
					resVSampl[name.str()]->Fill((*matching_cluster)->get_ampl(),residu);
					resVStot[name.str()]->Fill((*matching_cluster)->get_TOT(),residu);
					resVSsize[name.str()]->Fill((*matching_cluster)->get_size(),residu);
					absResVStime[name.str()]->Fill((*matching_cluster)->get_t(),Abs(residu));
					absResVSampl[name.str()]->Fill((*matching_cluster)->get_ampl(),Abs(residu));
					absResVStot[name.str()]->Fill((*matching_cluster)->get_TOT(),Abs(residu));
					absResVSsize[name.str()]->Fill((*matching_cluster)->get_size(),Abs(residu));
					if(residu<chisquare_threshold){
						muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						ampl_h[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_ampl());
					}
					delete *matching_cluster;
					current_clusters.erase(matching_cluster);
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
				it->second->cd(1);
				MM_residus[it->first]->Draw();
				for(int i=1;i<=nbins_2D;i++){
					for(int j=1;j<=nbins_2D;j++){
						int binN = muon_total[it->first]->GetBin(i,j);
						double binContent = 0;
						if(muon_total[it->first]->GetBinContent(binN) > 0) binContent = (muon_seen[it->first]->GetBinContent(binN))/(muon_total[it->first]->GetBinContent(binN));
						efficacity_2D[it->first]->SetBinContent(binN,binContent);
					}
				}
				it->second->cd(2);
				efficacity_2D[it->first]->Draw("COLZ");
				it->second->cd(3);
				if(point_nb[it->first]>0) correlation[it->first]->Draw("AP");
				it->second->cd(4);
				angle_alignment[it->first]->Draw();
				it->second->cd(5);
				resVSpos[it->first]->Draw();
				it->second->cd(6);
				resVSangle[it->first]->Draw();
				it->second->cd(7);
				resVSanglePerp[it->first]->Draw();
				it->second->cd(8);
				resVStime[it->first]->Draw();
				it->second->cd(9);
				resVSampl[it->first]->Draw();
				it->second->cd(10);
				resVStot[it->first]->Draw();
				it->second->cd(11);
				resVSsize[it->first]->Draw();
				it->second->cd(12);
				absResVSsize[it->first]->Draw();
				it->second->cd(13);
				absResVStime[it->first]->Draw();
				it->second->cd(14);
				absResVSampl[it->first]->Draw();
				it->second->cd(15);
				absResVStot[it->first]->Draw();
				it->second->cd(16);
				absResVSabsAngle[it->first]->Draw();
				//ampl_h[it->first]->Draw("COLZ");
				it->second->Modified();
				it->second->Update();
			}
			c0->cd(1);
			ray_clus_n->Draw();
			c0->cd(2);
			chisquares->Draw();
			c0->cd(3);
			ray_slope->Draw();
			ray_slope_X->Draw("SAME");
			ray_slope_Y->Draw("SAME");
			c0->cd(4);
			ray_phi->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << nentries << endl;
	for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
		cout << endl;
		cout << endl;
		it->second->cd(1);
		MM_residus[it->first]->Draw();
		cout << "resolution : " << endl;
		MM_residus[it->first]->Fit(offset_fit[it->first],"R");
		double total_seen = 0;
		double total_passed = 0;
		for(int i=1;i<=nbins_2D;i++){
			for(int j=1;j<=nbins_2D;j++){
				int binN = muon_total[it->first]->GetBin(i,j);
				double binContent = 0;
				if(muon_total[it->first]->GetBinContent(binN) > 0) binContent = (muon_seen[it->first]->GetBinContent(binN))/(muon_total[it->first]->GetBinContent(binN));
				efficacity_2D[it->first]->SetBinContent(binN,binContent);
				if(muon_seen[it->first]->GetBinContent(binN) > 0) ampl_h[it->first]->SetBinContent(binN,(ampl_h[it->first]->GetBinContent(binN))/(muon_seen[it->first]->GetBinContent(binN)));
				double pos_X = muon_total[it->first]->GetXaxis()->GetBinCenter(i);
				double pos_Y = muon_total[it->first]->GetYaxis()->GetBinCenter(j);
				if(pos_X<=2*Tomography::get_instance()->get_XY_size()/5. && pos_X>=-2*Tomography::get_instance()->get_XY_size()/5. && pos_Y<=2*Tomography::get_instance()->get_XY_size()/5. && pos_Y>=-2*Tomography::get_instance()->get_XY_size()/5.){
					total_seen += muon_seen[it->first]->GetBinContent(binN);
					total_passed += muon_total[it->first]->GetBinContent(binN);
				}
			}
		}
		efficacity[it->first] = total_seen/total_passed;
		it->second->cd(2);
		efficacity_2D[it->first]->Draw("COLZ");
		it->second->cd(3);
		if(point_nb[it->first]>0) correlation[it->first]->Draw("AP");
		it->second->cd(4);
		cout << "angle Z : " << endl;
		angle_alignment[it->first]->Fit(angle_z_fit[it->first],"R");
		angle_alignment[it->first]->Draw();
		it->second->cd(5);
		resVSpos[it->first]->Draw();
		it->second->cd(6);
		cout << "angle XY : " << endl;
		resVSangle[it->first]->Fit(angle_xy_fit[it->first],"R");
		cout << "-b/2a : " << angle_xy_fit[it->first]->GetParameter(1)/(2.*angle_xy_fit[it->first]->GetParameter(2)) << endl;
		resVSangle[it->first]->Draw();
		it->second->cd(7);
		resVSanglePerp[it->first]->Draw();
		it->second->cd(8);
		resVStime[it->first]->Draw();
		it->second->cd(9);
		resVSampl[it->first]->Draw();
		it->second->cd(10);
		resVStot[it->first]->Draw();
		it->second->cd(11);
		resVSsize[it->first]->Draw();
		it->second->cd(12);
		absResVSsize[it->first]->Draw();
		it->second->cd(13);
		absResVStime[it->first]->Draw();
		it->second->cd(14);
		absResVSampl[it->first]->Draw();
		it->second->cd(15);
		absResVStot[it->first]->Draw();
		it->second->cd(16);
		absResVSabsAngle[it->first]->Draw();
		//ampl_h[it->first]->Draw("COLZ");
		it->second->Modified();
		it->second->Update();
		cout << it->first << " efficacity : " << 100.*efficacity[it->first] << "%" << endl;
	}
	c0->cd(1);
	ray_clus_n->Draw();
	c0->cd(2);
	chisquares->Draw();
	c0->cd(3);
	ray_slope->Draw();
	ray_slope_X->Draw("SAME");
	ray_slope_Y->Draw("SAME");
	c0->cd(4);
	ray_phi->Draw();
	c0->Modified();
	c0->Update();
}
#include "MT_tomography.h"
#include "task/read_analyse_task.h"
#include "task/tracking_task.h"
void Analyse::Residus_ref_MT(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	Display_Thread * MT_display = Display_Thread::get_instance();
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	map<string,TCanvas*> c_MM;
	map<string,TH1D*> MM_residus;
	map<string,TF1*> offset_fit;
	map<string,TH2D*> muon_seen;
	map<string,TH2D*> muon_total;
	map<string,TH2D*> efficacity_2D;
	map<string,TH2D*> ampl_h;
	map<string,TGraph*> correlation;
	map<string,TProfile*> angle_alignment;
	map<string,TF1*> angle_z_fit;
	map<string,TProfile*> resVSpos;
	map<string,TProfile*> resVSampl;
	map<string,TProfile*> resVStime;
	map<string,TProfile*> resVSangle;
	map<string,TF1*> angle_xy_fit;
	map<string,TProfile*> resVSanglePerp;
	map<string,TProfile*> resVStot;
	map<string,TProfile*> resVSsize;
	map<string,TProfile*> absResVStot;
	map<string,TProfile*> absResVSsize;
	map<string,TProfile*> absResVSampl;
	map<string,TProfile*> absResVStime;
	map<string,TProfile*> absResVSabsAngle;
	map<string,int> point_nb;
	map<string,double> efficacity;
	int nbins = 200;
	double marge = 1./10.;
	int nbins_2D = 100*(1+2*marge);
	long eventReconstructed = 0;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	unsigned long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
			c_MM[name.str()]->Divide(4,4);
			MT_display->register_canvas(c_MM[name.str()], 4*4);
			MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
			MT_display->register_plot(MM_residus[name.str()],name.str(),"",1);
			muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			MT_display->register_plot(efficacity_2D[name.str()],name.str(),"COLZ",2);
			ampl_h[name.str()] = new TH2D((name.str()+"_ampl_mean").c_str(),(name.str()+"_ampl_mean").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
			MT_display->register_plot(ampl_h[name.str()],name.str(),"COLZ",16);
			efficacity_2D[name.str()]->SetStats(false);
			correlation[name.str()] = new TGraph();
			MT_display->register_plot(correlation[name.str()],name.str(),"AP",3);
			angle_alignment[name.str()] = new TProfile((name.str()+"_resVSperpPos").c_str(),(name.str()+"_resVSperpPos").c_str(),500,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,-5,5);
			MT_display->register_plot(angle_alignment[name.str()],name.str(),"",4);
			resVSpos[name.str()] = new TProfile((name.str()+"_resVSpos").c_str(),(name.str()+"_resVSpos").c_str(),500,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,-5,5);
			MT_display->register_plot(resVSpos[name.str()],name.str(),"",5);
			resVSampl[name.str()] = new TProfile((name.str()+"_resVSampl").c_str(),(name.str()+"_resVSampl").c_str(),500,-100,10000,-5,5);
			MT_display->register_plot(resVSampl[name.str()],name.str(),"",9);
			resVStime[name.str()] = new TProfile((name.str()+"_resVStime").c_str(),(name.str()+"_resVStime").c_str(),38,-2,34,-5,5);
			MT_display->register_plot(resVStime[name.str()],name.str(),"",8);
			resVSangle[name.str()] = new TProfile((name.str()+"_resVSangle").c_str(),(name.str()+"_resVSangle").c_str(),50,-0.6,0.6,-5,5);
			MT_display->register_plot(resVSangle[name.str()],name.str(),"",6);
			resVSanglePerp[name.str()] = new TProfile((name.str()+"_resVSanglePerp").c_str(),(name.str()+"_resVSanglePerp").c_str(),50,-0.6,0.6,-5,5);
			MT_display->register_plot(resVSanglePerp[name.str()],name.str(),"",7);
			resVStot[name.str()] = new TProfile((name.str()+"_resVStot").c_str(),(name.str()+"_resVStot").c_str(),26,0,25,-5,5);
			MT_display->register_plot(resVStot[name.str()],name.str(),"",10);
			resVSsize[name.str()] = new TProfile((name.str()+"_resVSsize").c_str(),(name.str()+"_resVSsize").c_str(),50,0,50,-5,5);
			MT_display->register_plot(resVSsize[name.str()],name.str(),"",11);
			absResVSampl[name.str()] = new TProfile((name.str()+"_absResVSampl").c_str(),(name.str()+"_absResVSampl").c_str(),500,-100,10000,0,5);
			MT_display->register_plot(absResVSampl[name.str()],name.str(),"",14);
			absResVStime[name.str()] = new TProfile((name.str()+"_absResVStime").c_str(),(name.str()+"_absResVStime").c_str(),38,-2,34,0,5);
			MT_display->register_plot(absResVStime[name.str()],name.str(),"",13);
			//absResVSabsAngle[name.str()] = new TProfile((name.str()+"_absResVSabsAngle").c_str(),(name.str()+"_absResVSabsAngle").c_str(),50,0,0.6,0,5);
			absResVStot[name.str()] = new TProfile((name.str()+"_absResVStot").c_str(),(name.str()+"_absResVStot").c_str(),26,0,25,0,5);
			MT_display->register_plot(absResVStot[name.str()],name.str(),"",15);
			absResVSsize[name.str()] = new TProfile((name.str()+"_absResVSsize").c_str(),(name.str()+"_absResVSsize").c_str(),50,0,50,0,5);
			MT_display->register_plot(absResVSsize[name.str()],name.str(),"",12);
			point_nb[name.str()] = 0;
			efficacity[name.str()] = 0;
			offset_fit[name.str()] = new TF1("offset_fit","[3]*exp(-(x-[0])*(x-[0])/(2*[1]*[1])) + [4]*exp(-(x-[0])*(x-[0])/(2*[2]*[2]))",-5,5);
			offset_fit[name.str()]->SetParameters(0,0.5,2,nentries/10);
			offset_fit[name.str()]->SetParLimits(0,-10,10);
			offset_fit[name.str()]->SetParLimits(1,0,1);
			offset_fit[name.str()]->SetParLimits(2,0,10);
			offset_fit[name.str()]->SetParLimits(3,1,nentries);
			offset_fit[name.str()]->SetParLimits(4,1,nentries);
			angle_z_fit[name.str()] = new TF1("angle_z_fit","pol1(0)",-150,150);
			angle_z_fit[name.str()]->SetParameters(0,0);
			angle_z_fit[name.str()]->SetParLimits(0,-5,5);
			angle_z_fit[name.str()]->SetParLimits(1,-1,1);
			angle_xy_fit[name.str()] = new TF1("angle_xy_fit","pol2(0)",-0.3,0.3);
			angle_xy_fit[name.str()]->SetParameters(0,0,0);
			angle_xy_fit[name.str()]->SetParLimits(0,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(1,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(2,-20,20);
			//nref_is_X[name.str()] = (*it)->get_is_X();
			if((*it)->get_is_X()) nref_x_n++;
		}
		//if((*it)->get_is_X()) det_x_n++;
	}
	/*
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	*/
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
		}
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2,2);
	MT_display->register_canvas(c0, 2*2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	MT_display->register_plot(chisquares,"stats","",2);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",get_det_N_tot() + 2,0,get_det_N_tot() + 2);
	MT_display->register_plot(ray_clus_n,"stats","",1);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	MT_display->register_plot(ray_slope,"stats","",3);
	TH1D * ray_phi = new TH1D("phi","phi",100,-Pi(),Pi());
	MT_display->register_plot(ray_phi,"stats","",4);
	TH1D * ray_slope_X = new TH1D("slope_X","slope_X",100,0,1);
	MT_display->register_plot(ray_slope_X,"stats","SAME",3);
	TH1D * ray_slope_Y = new TH1D("slope_Y","slope_Y",100,0,1);
	MT_display->register_plot(ray_slope_Y,"stats","SAME",3);
	ray_slope->SetLineColor(1);
	ray_slope_X->SetLineColor(2);
	ray_slope_Y->SetLineColor(3);


	Buffer_Task<ray_data> * ray_list = new Buffer_Task<ray_data>();
	Input_Task * to_do = new Read_Analyse_Task(nentries,this,this, new Tracking_Abs_Task(this, ray_list));
	vector<Thread*> threads;
	threads.push_back(new Reader_Thread(to_do));
	(threads.back())->start();
	const unsigned short n_thread = (Tomography::get_instance()->get_thread_number() > threads.size()) ? (Tomography::get_instance()->get_thread_number() - threads.size()) : 1;
	*MT_display << "1 | " << n_thread << "\n";
	for(unsigned short i=0;i<n_thread;i++){
		threads.push_back(new Worker_Thread());
		(threads.back())->start();
	}
	MT_display->start_count();
	//MT_display << Tomography::get_instance()->init_count() << "|" << setw(7) << "tracks\n";
	bool has_working_thread = true;
	unsigned long jentry = 0;
	while((has_working_thread || ray_list->can_fetch_data()) && Tomography::get_instance()->get_can_continue()){
		has_working_thread = false;
		for(unsigned short i=0;i<threads.size();i++){
			if(threads[i]->is_working()){
				has_working_thread = true;
				break;
			}
		}
		if(!(ray_list->can_fetch_data())){
			usleep(1000);
			continue;
		}
		ray_data * current_rays = ray_list->fetch_data();
		vector<Ray>::iterator ray_it = (current_rays->rays).begin();
		while(ray_it != (current_rays->rays).end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				chisquares->Fill(ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y());
				ray_clus_n->Fill(ray_it->get_clus_n());
				double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
				ray_slope->Fill(ATan(slope));
				ray_slope_X->Fill(ATan(Abs(ray_it->get_slope_X())));
				ray_slope_Y->Fill(ATan(Abs(ray_it->get_slope_Y())));
				double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
				if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
				ray_phi->Fill(phi);
				++ray_it;
			}
			else ray_it = (current_rays->rays).erase(ray_it);
		}
		eventReconstructed+=(current_rays->rays).size();
		for(vector<Event*>::iterator it = (current_rays->CBevent->events).begin();it!=(current_rays->CBevent->events).end();++it){
			if(!((*it)->get_is_ref())){
				ostringstream name;
				name <<(*it)->get_type() << "_" << (*it)->get_n_in_tree();
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				for(vector<Ray>::iterator jt=(current_rays->rays).begin();jt!=(current_rays->rays).end();++jt){
					//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
					//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
					//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
					//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
					if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
					if(jt->get_clus_n()<static_cast<unsigned int>(get_layers_n()-non_ref_n)) continue;
					double residu = numeric_limits<double>::max();
					vector<Cluster*>::iterator matching_cluster = current_clusters.end();
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						(*kt)->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(*kt);
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
					if(matching_cluster == current_clusters.end()) continue;
					(*matching_cluster)->set_perp_pos_mm(*jt);
					if((*matching_cluster)->get_is_X()){
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),(*matching_cluster)->get_pos_mm());
						angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						resVSpos[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
						//absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_X()),Abs(residu));
						resVSanglePerp[name.str()]->Fill(jt->get_slope_Y(),residu);
					}
					else{
						correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_pos_mm());
						angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						resVSpos[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
						//absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_Y()),Abs(residu));
						resVSanglePerp[name.str()]->Fill(jt->get_slope_X(),residu);
					}
					point_nb[name.str()]++;
					MM_residus[name.str()]->Fill(residu);
					resVStime[name.str()]->Fill((*matching_cluster)->get_t(),residu);
					resVSampl[name.str()]->Fill((*matching_cluster)->get_ampl(),residu);
					resVStot[name.str()]->Fill((*matching_cluster)->get_TOT(),residu);
					resVSsize[name.str()]->Fill((*matching_cluster)->get_size(),residu);
					absResVStime[name.str()]->Fill((*matching_cluster)->get_t(),Abs(residu));
					absResVSampl[name.str()]->Fill((*matching_cluster)->get_ampl(),Abs(residu));
					absResVStot[name.str()]->Fill((*matching_cluster)->get_TOT(),Abs(residu));
					absResVSsize[name.str()]->Fill((*matching_cluster)->get_size(),Abs(residu));
					if(residu<chisquare_threshold){
						muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						ampl_h[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()),(*matching_cluster)->get_ampl());
					}
					int binN = efficacity_2D[name.str()]->FindBin(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
					efficacity_2D[name.str()]->SetBinContent(binN,(muon_seen[name.str()]->GetBinContent(binN))/(muon_total[name.str()]->GetBinContent(binN)));
					delete *matching_cluster;
					current_clusters.erase(matching_cluster);
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
			}
		}
		delete current_rays;

		//if(jentry%100) MT_display << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << eventReconstructed;
		jentry++;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()) MT_display->display_canvas();
	}
	for(unsigned short i=0;i<threads.size();i++){
		threads[i]->stop();
		delete threads[i];
	}

	//MT_display << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << eventReconstructed << "\n";
	MT_display->stop();

	for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
		cout << "\n\n";
		cout << "resolution : " << "\n";
		MM_residus[it->first]->Fit(offset_fit[it->first],"R");
		double total_seen = 0;
		double total_passed = 0;
		for(int i=1;i<=nbins_2D;i++){
			for(int j=1;j<=nbins_2D;j++){
				int binN = muon_total[it->first]->GetBin(i,j);
				if(muon_seen[it->first]->GetBinContent(binN) > 0) ampl_h[it->first]->SetBinContent(binN,(ampl_h[it->first]->GetBinContent(binN))/(muon_seen[it->first]->GetBinContent(binN)));
				double pos_X = muon_total[it->first]->GetXaxis()->GetBinCenter(i);
				double pos_Y = muon_total[it->first]->GetYaxis()->GetBinCenter(j);
				if(pos_X<=2*Tomography::get_instance()->get_XY_size()/5. && pos_X>=-2*Tomography::get_instance()->get_XY_size()/5. && pos_Y<=2*Tomography::get_instance()->get_XY_size()/5. && pos_Y>=-2*Tomography::get_instance()->get_XY_size()/5.){
					total_seen += muon_seen[it->first]->GetBinContent(binN);
					total_passed += muon_total[it->first]->GetBinContent(binN);
				}
			}
		}
		efficacity[it->first] = total_seen/total_passed;
		cout << "angle Z : " << "\n";
		angle_alignment[it->first]->Fit(angle_z_fit[it->first],"R");
		cout << "angle XY : " << "\n";
		resVSangle[it->first]->Fit(angle_xy_fit[it->first],"R");
		cout << "-b/2a : " << angle_xy_fit[it->first]->GetParameter(1)/(2.*angle_xy_fit[it->first]->GetParameter(2)) << "\n";
		cout << it->first << " efficacity : " << 100.*efficacity[it->first] << "%" << endl;
	}
	MT_display->display_canvas();
}
double Analyse::Residus_ref_cost(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	map<string,TH1D*> MM_residus;
	map<string,TF1*> offset_fit;
	map<string,TProfile*> angle_alignment;
	map<string,TF1*> angle_z_fit;
	map<string,TProfile*> resVSangle;
	map<string,TF1*> angle_xy_fit;
	int nbins = 200;
	double marge = 1./10.;
	long eventReconstructed = 0;
	double eventSuitable = 0;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
			angle_alignment[name.str()] = new TProfile((name.str()+"_angle").c_str(),(name.str()+"_angle").c_str(),500,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,-5,5);
			resVSangle[name.str()] = new TProfile((name.str()+"_resVSangle").c_str(),(name.str()+"_resVSangle").c_str(),50,-0.6,0.6,-5,5);
			offset_fit[name.str()] = new TF1("offset_fit","[3]*exp(-(x-[0])*(x-[0])/(2*[1]*[1])) + [4]*exp(-(x-[0])*(x-[0])/(2*[2]*[2]))",-5,5);
			offset_fit[name.str()]->SetParameters(0,0.5,2,nentries/10);
			offset_fit[name.str()]->SetParLimits(0,-10,10);
			offset_fit[name.str()]->SetParLimits(1,0,1);
			offset_fit[name.str()]->SetParLimits(2,0,10);
			offset_fit[name.str()]->SetParLimits(3,1,nentries);
			offset_fit[name.str()]->SetParLimits(4,1,nentries);
			angle_z_fit[name.str()] = new TF1("angle_z_fit","pol1(0)",-150,150);
			angle_z_fit[name.str()]->SetParameters(0,0);
			angle_z_fit[name.str()]->SetParLimits(0,-5,5);
			angle_z_fit[name.str()]->SetParLimits(1,-1,1);
			angle_xy_fit[name.str()] = new TF1("angle_xy_fit","pol2(0)",-0.3,0.3);
			angle_xy_fit[name.str()]->SetParameters(0,0,0);
			angle_xy_fit[name.str()]->SetParLimits(0,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(1,-5,5);
			angle_xy_fit[name.str()]->SetParLimits(2,-20,20);
			//nref_is_X[name.str()] = (*it)->get_is_X();
			if((*it)->get_is_X()) nref_x_n++;
		}
		//if((*it)->get_is_X()) det_x_n++;
	}
	/*
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	*/
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return 0;
		}
	}
	if(non_ref_n>2){
		cout << "too many non ref det" << endl;
	}
	if(nref_x_n>1){
		cout << "too many non ref det" << endl;
	}

	if (fChain == 0) return 0;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				++ray_it;
			}
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(get_det_N_tot());
		
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				ostringstream name;
				name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
					//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
					//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
					//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
					//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
					if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
					if(jt->get_clus_n()<static_cast<unsigned int>(get_layers_n()-non_ref_n)) continue;
					double residu = numeric_limits<double>::max();
					vector<Cluster*>::iterator matching_cluster = current_clusters.end();
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						(*kt)->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(*kt);
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					if(matching_cluster == current_clusters.end()) continue;
					(*matching_cluster)->set_perp_pos_mm(*jt);
					if((*matching_cluster)->get_is_X()){
						angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
					}
					else{
						angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
					}
					delete *matching_cluster;
					current_clusters.erase(matching_cluster);
					MM_residus[name.str()]->Fill(residu);
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << nentries << endl;
	double cost = 0;
	for(map<string,TH1D*>::iterator it = MM_residus.begin();it!=MM_residus.end();++it){
		//MM_residus[it->first]->Fit(offset_fit[it->first],"RNQ");
		angle_alignment[it->first]->Fit(angle_z_fit[it->first],"RNQ");
		resVSangle[it->first]->Fit(angle_xy_fit[it->first],"RNQ");
		//cost += (offset_fit[it->first]->GetParameter(0))*(offset_fit[it->first]->GetParameter(0));
		cost += (MM_residus[it->first]->GetMean())*(MM_residus[it->first]->GetMean());
		cost += (angle_z_fit[it->first]->GetParameter(1))*(angle_z_fit[it->first]->GetParameter(1));
		cost += (angle_z_fit[it->first]->GetParameter(2))*(angle_z_fit[it->first]->GetParameter(2));
		cost += (angle_xy_fit[it->first]->GetParameter(1))*(angle_xy_fit[it->first]->GetParameter(1));
		delete MM_residus[it->first];
		delete angle_alignment[it->first];
		delete resVSangle[it->first];
		delete offset_fit[it->first];
		delete angle_z_fit[it->first];
		delete angle_xy_fit[it->first];
	}
	return cost;
}
void Analyse::Residus_ref_2D(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas* c_MM;
	TH2D* muon_seen;
	TH2D* muon_total;
	TH2D* efficacity_2D;
	double efficacity = 0;
	int nbins = 200;
	double marge = 1./10.;
	int nbins_2D = 100*(1+2*marge);
	long eventReconstructed = 0;
	double eventSuitable = 0;
	unsigned int nref_x_n = 0;
	unsigned int det_x_n = 0;
	unsigned int det_num = get_det_N_tot();
	vector<double> det_z;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	ostringstream name;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			if(name.str().size()>0) name << "_";
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			if((*it)->get_is_X()) nref_x_n++;
			det_z.push_back((*it)->get_z());
		}
		if((*it)->get_is_X()) det_x_n++;
	}
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "you can only calculate 2D efficacity for 2D detector" << endl;
			return;
		}
	}
	if(det_z.size()>2){
		cout << "this is meant to calculate the efficacity of only one 2D detector" << endl;
		return;
	}
	double nref_z = det_z[0];
	if(det_z[1]!=nref_z){
		cout << "problem in z coordinate of the 2D detector" << endl;
		return;
	}

	c_MM = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
	c_MM->Divide(3);
	muon_seen = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
	muon_total = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
	efficacity_2D = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2,nbins_2D,-(1+marge)*Tomography::get_instance()->get_XY_size()/2,(1+marge)*Tomography::get_instance()->get_XY_size()/2);
	efficacity_2D->SetStats(false);

	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",get_det_N_tot() + 2,0,get_det_N_tot() + 2);

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		/*
		for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
			if((jt->get_chiSquare_X()+jt->get_chiSquare_Y()) < chisquare_threshold){
				chisquares->Fill(jt->get_chiSquare_X()+jt->get_chiSquare_Y());
				ray_clus_n->Fill(jt->get_clus_n());
			}
		}
		*/
		//eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(get_det_N_tot());
		
		vector<Event*> nref_event;
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				nref_event.push_back((*it)->Clone());
			}
		}
		if(nref_event.size()!=2){
			cout << "problem in event size" << endl;
			return;
		}
		delete currentCBEvent;
		for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
			if(jt->get_clus_n()<(det_num-4)){
				continue;
			}
			if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > 2.){
				continue;
			}
			eventReconstructed++;
			chisquares->Fill(jt->get_chiSquare_X()+jt->get_chiSquare_Y());
			ray_clus_n->Fill(jt->get_clus_n());
			vector<unsigned int> seen_clus_in_array;
			for(vector<Event*>::iterator it = nref_event.begin();it!=nref_event.end();++it){
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				double residu = numeric_limits<double>::max();
				vector<Cluster*>::iterator matching_cluster = current_clusters.end();
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					(*kt)->set_perp_pos_mm(*jt);
					double current_residu = jt->get_residu_ref(*kt);
					if(current_residu<residu){
						residu = current_residu;
						matching_cluster = kt;
					}
				}
				if(matching_cluster == current_clusters.end()) continue;
				(*matching_cluster)->set_perp_pos_mm(*jt);
				if(residu<chisquare_threshold){
					seen_clus_in_array.push_back(matching_cluster - current_clusters.begin());
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
				if(seen_clus_in_array.size()==2){
					muon_seen->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
					for(int i_event=0;i_event<2;i_event++){
						delete (nref_event[i_event]->clusters)[seen_clus_in_array[i_event]];
						(nref_event[i_event]->clusters).erase((nref_event[i_event]->clusters).begin()+seen_clus_in_array[i_event]);
					}
				}
			}
			muon_total->Fill(jt->eval_X(nref_z),jt->eval_Y(nref_z));
		}
		for(vector<Event*>::iterator it = (nref_event).begin();it!=(nref_event).end();++it){
			delete *it;
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(int i=1;i<=nbins_2D;i++){
				for(int j=1;j<=nbins_2D;j++){
					int binN = muon_total->GetBin(i,j);
					double binContent = 0;
					if(muon_total->GetBinContent(binN) > 0) binContent = (muon_seen->GetBinContent(binN))/(muon_total->GetBinContent(binN));
					efficacity_2D->SetBinContent(binN,binContent);
				}
			}
			c_MM->cd(1);
			efficacity_2D->Draw("COLZ");
			c_MM->cd(2);
			muon_seen->Draw("COLZ");
			c_MM->cd(3);
			muon_total->Draw("COLZ");
			c_MM->Modified();
			c_MM->Update();
			c0->cd(1);
			chisquares->Draw();
			c0->cd(2);
			ray_clus_n->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << nentries << endl;
	double total_seen = 0;
	double total_passed = 0;
	for(int i=1;i<=nbins_2D;i++){
		for(int j=1;j<=nbins_2D;j++){
			int binN = muon_total->GetBin(i,j);
			double binContent = 0;
			if(muon_total->GetBinContent(binN) > 0) binContent = (muon_seen->GetBinContent(binN))/(muon_total->GetBinContent(binN));
			efficacity_2D->SetBinContent(binN,binContent);
			double pos_X = muon_total->GetXaxis()->GetBinCenter(i);
			double pos_Y = muon_total->GetYaxis()->GetBinCenter(j);
			if(pos_X<=2*Tomography::get_instance()->get_XY_size()/5. && pos_X>=-2*Tomography::get_instance()->get_XY_size()/5. && pos_Y<=2*Tomography::get_instance()->get_XY_size()/5. && pos_Y>=-2*Tomography::get_instance()->get_XY_size()/5.){
				total_seen += muon_seen->GetBinContent(binN);
				total_passed += muon_total->GetBinContent(binN);
			}
		}
	}
	efficacity = total_seen/total_passed;
	c_MM->cd(1);
	efficacity_2D->Draw("COLZ");
	c_MM->cd(2);
	muon_seen->Draw("COLZ");
	c_MM->cd(3);
	muon_total->Draw("COLZ");
	c_MM->Modified();
	c_MM->Update();
	cout << name.str() << " efficacity : " << 100.*efficacity << "%" << endl;
	c0->cd(1);
	chisquares->Draw();
	c0->cd(2);
	ray_clus_n->Draw();
	c0->Modified();
	c0->Update();
}
void Analyse::Efficacity(){
	double chisquare_threshold = 10;
	/*
	int non_ref_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	*/
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	map<string,TCanvas*> c_MM;
	map<string,TProfile*> MM_eff;
	long eventReconstructed = 0;
	double eventSuitable = 0;
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	LoadTree(nentries-1);
	fChain->GetEntry(nentries-1);
	double maxtime = evttime;
	int nbins = nentries/100;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			if((*it)->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>((*it)->get_type(),(*it)->get_n_in_tree())] = pair<Tomography::det_type,int>((*it)->get_perp_type(),(*it)->get_perp_n());
			c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
			MM_eff[name.str()] = new TProfile((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins,0,maxtime,0,1);
			if((*it)->get_is_X()) nref_x_n++;
		}
		//if((*it)->get_is_X()) det_x_n++;
	}
	/*
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	*/
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
		}
	}
	if (fChain == 0) return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				++ray_it;
			}
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(get_det_N_tot());
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				ostringstream name;
				name <<(*it)->get_type() << "_" << (*it)->get_n_in_tree();
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
					//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
					//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
					//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
					//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
					if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
					if(jt->get_clus_n()<static_cast<unsigned int>(get_layers_n()-non_ref_n)) continue;
					double residu = numeric_limits<double>::max();
					vector<Cluster*>::iterator matching_cluster = current_clusters.end();
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						(*kt)->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(*kt);
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					if(matching_cluster == current_clusters.end()){
						MM_eff[name.str()]->Fill(evttime,0);
						continue;
					}
					delete *matching_cluster;
					current_clusters.erase(matching_cluster);
					if(residu<chisquare_threshold){
						MM_eff[name.str()]->Fill(evttime,1);
					}
					else{
						MM_eff[name.str()]->Fill(evttime,0);
					}
				}
				for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
					delete *kt;
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
				it->second->cd();
				MM_eff[it->first]->Draw();
				it->second->Modified();
				it->second->Update();
			}
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << nentries << endl;
	for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
		cout << endl;
		cout << endl;
		it->second->cd();
		MM_eff[it->first]->Draw();
		it->second->Modified();
		it->second->Update();
	}
}
void Analyse::ExportAbsorptionRays(string outFileName){
	long eventReconstructed = 0;
	long eventSuitable = 0;
	double chisquare_threshold = 100;

	double Z_Up = numeric_limits<double>::min();
	double Z_Down = numeric_limits<double>::max();
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		double current_z = (*det_it)->get_z();
		if(current_z>Z_Up) Z_Up = current_z;
		if(current_z<Z_Down) Z_Down = current_z;
	}

	Tray * rayTree = new Tray(outFileName);

	//if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		rayTree->fillTree(evn,evttime,currentRays,Z_Up,Z_Down);
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	rayTree->Write();
	rayTree->CloseFile();
}
TH2D * Analyse::AbsorptionFluxMap(double z, TCanvas * c1, double y_angle, int nbins){
	long eventReconstructed = 0;
	long eventSuitable = 0;
	double chisquare_threshold = 100;
	if(Abs(y_angle)>Pi()/2.){
		return new TH2D("error","error",1,0,1,1,0,1);
	}

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	//double z_Pb = 1553;
	if(c1 == 0){
		c1 = new TCanvas("fluxMap","fluxMap");
	}
	double z_max = numeric_limits<double>::min();
	double z_min = numeric_limits<double>::max();
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		double current_z = (*det_it)->get_z();
		if(current_z>z_max) z_max = current_z;
		if(current_z<z_min) z_min = current_z;
	}

	double x_min = -Tomography::get_instance()->get_XY_size()/2.;
	double x_max = Tomography::get_instance()->get_XY_size()/2.;
	double y_min = -Tomography::get_instance()->get_XY_size()/2.;
	double y_max = Tomography::get_instance()->get_XY_size()/2.;
	Point orig(0,0,z);
	Point norm(0,Sin(y_angle),Cos(y_angle));
	Plane proj(norm,orig);
	cout << proj.get_a() << "*x + " << proj.get_b() << "*y + " << proj.get_c() << "*z + " << proj.get_d() << " = 0" << endl;
	if(z>z_max || z<z_min){
		Line first_line(Point(-Tomography::get_instance()->get_XY_size()/2.,-Tomography::get_instance()->get_XY_size()/2.,z_min),Point(Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,z_max));
		Line second_line(Point(Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,z_min),Point(-Tomography::get_instance()->get_XY_size()/2.,-Tomography::get_instance()->get_XY_size()/2.,z_max));
		Point corner_a = proj.intersection(first_line);
		Point corner_b = proj.intersection(second_line);
		x_max = Max(Abs((corner_a-orig).get_X()),Abs((corner_b-orig).get_X()));
		y_max = Max((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
		x_min = -x_max;
		y_min = Min((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
	}
	y_max = y_max/Cos(y_angle);
	y_min = y_min/Cos(y_angle);
	double width_x = x_max - x_min;
	double width_y = y_max - y_min;
	x_max += 0.05*width_x;
	x_min -= 0.05*width_x;
	y_max += 0.05*width_y;
	y_min -= 0.05*width_y;
	cout << x_min << " " << x_max << " " << y_min << " " << y_max << endl;

	//if (fChain == 0) return fluxMapZ;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	int nbins_2D = 0;
	if(nbins>0) nbins_2D = nbins;
	else nbins_2D = Sqrt(0.02*nentries);
	TH2D * fluxMapZ = new TH2D("fluxMapZ","fluxMapZ",nbins_2D,x_min,x_max,nbins_2D,y_min,y_max);
	fluxMapZ->SetStats(0);

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			Point current_point = it->eval_plane(proj);
			fluxMapZ->Fill(current_point.get_X(),current_point.get_Y()/Cos(y_angle));
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			c1->cd();
			fluxMapZ->Draw("COLZ");
			c1->Modified();
			c1->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	c1->cd();
	fluxMapZ->Draw("COLZ");
	c1->Modified();
	c1->Update();
	return fluxMapZ;
}
void Analyse::WatToFluxMap(double z,TEllipse el, TCanvas * c1, TCanvas * c2, double y_angle){
	long eventReconstructed = 0;
	long eventSuitable = 0;
	double chisquare_threshold = 100;
	unsigned int interval_length = 50000;
	if(Abs(y_angle)>Pi()/2.){
		return;
	}

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	//double z_Pb = 1553;
	if(c1 == 0){
		c1 = new TCanvas("fluxMap","fluxMap");
	}
	c1->Divide(2);
	double z_max = numeric_limits<double>::min();
	double z_min = numeric_limits<double>::max();
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		double current_z = (*det_it)->get_z();
		if(current_z>z_max) z_max = current_z;
		if(current_z<z_min) z_min = current_z;
	}

	double x_min = -Tomography::get_instance()->get_XY_size()/2.;
	double x_max = Tomography::get_instance()->get_XY_size()/2.;
	double y_min = -Tomography::get_instance()->get_XY_size()/2.;
	double y_max = Tomography::get_instance()->get_XY_size()/2.;
	Point orig(0,0,z);
	Point norm(0,Sin(y_angle),Cos(y_angle));
	Plane proj(norm,orig);
	cout << proj.get_a() << "*x + " << proj.get_b() << "*y + " << proj.get_c() << "*z + " << proj.get_d() << " = 0" << endl;
	if(z>z_max || z<z_min){
		Line first_line(Point(-Tomography::get_instance()->get_XY_size()/2.,-Tomography::get_instance()->get_XY_size()/2.,z_min),Point(Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,z_max));
		Line second_line(Point(Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,z_min),Point(-Tomography::get_instance()->get_XY_size()/2.,-Tomography::get_instance()->get_XY_size()/2.,z_max));
		Point corner_a = proj.intersection(first_line);
		Point corner_b = proj.intersection(second_line);
		x_max = Max(Abs((corner_a-orig).get_X()),Abs((corner_b-orig).get_X()));
		y_max = Max((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
		x_min = -x_max;
		y_min = Min((corner_a-orig).get_Y(),(corner_b-orig).get_Y());
	}
	y_max = y_max/Cos(y_angle);
	y_min = y_min/Cos(y_angle);
	double width_x = x_max - x_min;
	double width_y = y_max - y_min;
	x_max += 0.05*width_x;
	x_min -= 0.05*width_x;
	y_max += 0.05*width_y;
	y_min -= 0.05*width_y;
	cout << x_min << " " << x_max << " " << y_min << " " << y_max << endl;

	//if (fChain == 0) return fluxMapZ;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	LoadTree(0);
	fChain->GetEntry(0);
	double evttime_start = evttime;
	LoadTree(nentries-1);
	fChain->GetEntry(nentries-1);
	double evttime_end = evttime;


	TH2D * fluxMapZ = new TH2D("fluxMapZ","fluxMapZ",Sqrt(0.02*nentries),x_min,x_max,Sqrt(0.02*nentries),y_min,y_max);
	fluxMapZ->SetStats(0);
	TEllipse * this_el = new TEllipse(el);
	this_el->SetFillStyle(0);
	double ellipse_angle = this_el->GetTheta()*Pi()/180.;
	double band_half_width = Sqrt((this_el->GetR1()*(this_el->GetR1())*Sin(ellipse_angle)*Sin(ellipse_angle))+(this_el->GetR2()*(this_el->GetR2())*Cos(ellipse_angle)*Cos(ellipse_angle)));
	TLine * up_line = new TLine(x_min,band_half_width+this_el->GetY1(),x_max,band_half_width+this_el->GetY1());
	TLine * down_line = new TLine(x_min,-band_half_width+this_el->GetY1(),x_max,-band_half_width+this_el->GetY1());
	unsigned int interval_n = 0;
	unsigned int event_n_interval = 0;
	unsigned int reconstructed_track = 0;
	unsigned int track_in_ellipse = 0;
	unsigned int track_in_band = 0;
	double vertical_band_half_width = Sqrt((this_el->GetR1()*(this_el->GetR1())*Cos(ellipse_angle)*Cos(ellipse_angle))+(this_el->GetR2()*(this_el->GetR2())*Sin(ellipse_angle)*Sin(ellipse_angle)));
	TLine * left_line = new TLine(-vertical_band_half_width+this_el->GetX1(),y_min,-vertical_band_half_width+this_el->GetX1(),y_max);
	TLine * right_line = new TLine(vertical_band_half_width+this_el->GetX1(),y_min,vertical_band_half_width+this_el->GetX1(),y_max);
	TH2D * tank_profile = new TH2D("tank_profile","tank_profile",nentries/interval_length,0,nentries/interval_length,Sqrt(0.02*nentries),y_min,y_max);
	unsigned int track_in_vertical_band = 0;
	if(c2 == 0){
		c2 = new TCanvas("WaterMon","WaterMon");
	}
	c2->Divide(2,2);
	TH1D * tank_tracks = new TH1D("tank_tracks","tank_tracks",nentries/interval_length,0,nentries);
	TH1D * tank_tracks_norm = new TH1D("tank_tracks_norm","tank_tracks_norm",nentries/interval_length,0,nentries/interval_length);
	TH1D * tank_tracks_norm_H = new TH1D("tank_tracks_norm_H","tank_tracks_norm_H",nentries/interval_length,0,nentries/interval_length);
	TH1D * tank_tracks_norm_V = new TH1D("tank_tracks_norm_V","tank_tracks_norm_V",nentries/interval_length,0,nentries/interval_length);
	TCanvas * c3 = new TCanvas("WaterMon_time","WaterMon_time");
	c3->Divide(1,4);
	double bin_time = 2*3600*Tomography::get_instance()->get_clock_rate();
	double time_bin_n = (evttime_end-evttime_start)/bin_time;
	TH1D * tank_tracks_time = new TH1D("tank_tracks_time","tank_tracks_time",time_bin_n,evttime_start,evttime_end);
	TH1D * out_tracks_time = new TH1D("out_tracks_time","out_tracks_time",time_bin_n,evttime_start,evttime_end);
	TH1D * ratio_tracks_time = new TH1D("ratio_tracks_time","ratio_tracks_time",time_bin_n,evttime_start,evttime_end);
	tank_tracks_time->Sumw2();
	out_tracks_time->Sumw2();
	ratio_tracks_time->Sumw2();
	double * interval_time = new double[(nentries/interval_length)+1];
	interval_time[0] = evttime_start;
	cout << setw(20) << "interval n" <<  "|" << setw(20) << "total track" <<  "|" << setw(20) << "track in ellipse" <<  "|" << setw(20) << "track in H band" <<  "|" << setw(20) << "track in V band" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			Point current_point = it->eval_plane(proj);
			Point_2D point_in_plane(current_point.get_X(),current_point.get_Y()/Cos(y_angle));
			fluxMapZ->Fill(point_in_plane.get_X(),point_in_plane.get_Y());
			reconstructed_track++;
			if(point_in_plane.is_inside(*this_el)){
				track_in_ellipse++;
				tank_tracks->Fill(jentry);
				tank_tracks_time->Fill(evttime);
			}
			else{
				out_tracks_time->Fill(evttime);
			}
			if(Abs(point_in_plane.get_Y() - this_el->GetY1()) < band_half_width) track_in_band++;
			if(Abs(point_in_plane.get_X() - this_el->GetX1()) < vertical_band_half_width){
				tank_profile->Fill(interval_n,point_in_plane.get_Y());
				track_in_vertical_band++;
			}
		}
		event_n_interval++;
		if(event_n_interval%interval_length == 0){
			cout << setw(20) << interval_n <<  "|" << setw(20) << reconstructed_track <<  "|" << setw(20) << track_in_ellipse <<  "|" << setw(20) << track_in_band <<  "|" << setw(20) << track_in_vertical_band << endl;

			for(int j=1;j<=tank_profile->GetNbinsY();j++){
				int binN = tank_profile->GetBin(interval_n,j);
				if(track_in_vertical_band>0) tank_profile->SetBinContent(binN,tank_profile->GetBinContent(binN)/track_in_vertical_band);
				else tank_profile->SetBinContent(binN,0);
			}
			if((reconstructed_track-track_in_ellipse)>0){
				tank_tracks_norm->SetBinContent(interval_n,track_in_ellipse*1./(reconstructed_track-track_in_ellipse));
				tank_tracks_norm->SetBinError(interval_n,Sqrt(track_in_ellipse*(reconstructed_track)*1./((reconstructed_track-track_in_ellipse)*(reconstructed_track-track_in_ellipse)*(reconstructed_track-track_in_ellipse))));
			}
			else{
				tank_tracks_norm->SetBinContent(interval_n,0);
				tank_tracks_norm->SetBinError(interval_n,0);
			}
			if((track_in_band-track_in_ellipse)>0){
				tank_tracks_norm_H->SetBinContent(interval_n,track_in_ellipse*1./(track_in_band-track_in_ellipse));
				tank_tracks_norm_H->SetBinError(interval_n,Sqrt(track_in_ellipse*(track_in_band)*1./((track_in_band-track_in_ellipse)*(track_in_band-track_in_ellipse)*(track_in_band-track_in_ellipse))));
			}
			else{
				tank_tracks_norm_H->SetBinContent(interval_n,0);
				tank_tracks_norm_H->SetBinError(interval_n,0);
			}
			if((track_in_vertical_band-track_in_ellipse)>0){
				tank_tracks_norm_V->SetBinContent(interval_n,track_in_ellipse*1./(track_in_vertical_band-track_in_ellipse));
				tank_tracks_norm_V->SetBinError(interval_n,Sqrt(track_in_ellipse*(track_in_vertical_band)*1./((track_in_vertical_band-track_in_ellipse)*(track_in_vertical_band-track_in_ellipse)*(track_in_vertical_band-track_in_ellipse))));
			}
			else{
				tank_tracks_norm_V->SetBinContent(interval_n,0);
				tank_tracks_norm_V->SetBinError(interval_n,0);
			}
			interval_n++;
			interval_time[interval_n] = evttime;
			event_n_interval = 0;
			reconstructed_track = 0;
			track_in_ellipse = 0;
			track_in_band = 0;
			track_in_vertical_band = 0;
		}
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			c1->cd(1);
			fluxMapZ->Draw("COLZ");
			this_el->Draw();
			up_line->Draw();
			down_line->Draw();
			left_line->Draw();
			right_line->Draw();
			c1->cd(2);
			tank_profile->Draw("COLZ");
			c1->Modified();
			c1->Update();
			c2->cd(1);
			tank_tracks->Draw("E");
			c2->cd(2);
			tank_tracks_norm->Draw("E");
			c2->cd(3);
			tank_tracks_norm_H->Draw("E");
			c2->cd(4);
			tank_tracks_norm_V->Draw("E");
			c2->Modified();
			c2->Update();
			ratio_tracks_time->Divide(tank_tracks_time,out_tracks_time);
			c3->cd(1);
			tank_tracks_time->Draw("E");
			c3->cd(2);
			out_tracks_time->Draw("E");
			c3->cd(3);
			ratio_tracks_time->Draw("E");
			c3->Modified();
			c3->Update();
		}
	}
	cout << setw(20) << interval_n <<  "|" << setw(20) << reconstructed_track <<  "|" << setw(20) << track_in_ellipse <<  "|" << setw(20) << track_in_band <<  "|" << setw(20) << track_in_vertical_band << endl;
	if(event_n_interval != 0){
		for(int j=1;j<=tank_profile->GetNbinsY();j++){
			int binN = tank_profile->GetBin(interval_n,j);
			if(track_in_vertical_band>0) tank_profile->SetBinContent(binN,tank_profile->GetBinContent(binN)/track_in_vertical_band);
			else tank_profile->SetBinContent(binN,0);
		}
		if((reconstructed_track-track_in_ellipse)>0){
			tank_tracks_norm->SetBinContent(interval_n,track_in_ellipse*1./(reconstructed_track-track_in_ellipse));
			tank_tracks_norm->SetBinError(interval_n,Sqrt(track_in_ellipse*(reconstructed_track)*1./((reconstructed_track-track_in_ellipse)*(reconstructed_track-track_in_ellipse)*(reconstructed_track-track_in_ellipse))));
		}
		else{
			tank_tracks_norm->SetBinContent(interval_n,0);
			tank_tracks_norm->SetBinError(interval_n,0);
		}
		if((track_in_band-track_in_ellipse)>0){
			tank_tracks_norm_H->SetBinContent(interval_n,track_in_ellipse*1./(track_in_band-track_in_ellipse));
			tank_tracks_norm_H->SetBinError(interval_n,Sqrt(track_in_ellipse*(track_in_band)*1./((track_in_band-track_in_ellipse)*(track_in_band-track_in_ellipse)*(track_in_band-track_in_ellipse))));
		}
		else{
			tank_tracks_norm_H->SetBinContent(interval_n,0);
			tank_tracks_norm_H->SetBinError(interval_n,0);
		}
		if((track_in_vertical_band-track_in_ellipse)>0){
			tank_tracks_norm_V->SetBinContent(interval_n,track_in_ellipse*1./(track_in_vertical_band-track_in_ellipse));
			tank_tracks_norm_V->SetBinError(interval_n,Sqrt(track_in_ellipse*(track_in_vertical_band)*1./((track_in_vertical_band-track_in_ellipse)*(track_in_vertical_band-track_in_ellipse)*(track_in_vertical_band-track_in_ellipse))));
		}
		else{
			tank_tracks_norm_V->SetBinContent(interval_n,0);
			tank_tracks_norm_V->SetBinError(interval_n,0);
		}
	}
	interval_time[nentries/interval_length] = evttime;
	c1->cd(1);
	fluxMapZ->Draw("COLZ");
	this_el->Draw();
	up_line->Draw();
	down_line->Draw();
	left_line->Draw();
	right_line->Draw();
	c1->cd(2);
	tank_profile->Draw("COLZ");
	c1->Modified();
	c1->Update();
	c2->cd(1);
	tank_tracks->Draw("E");
	c2->cd(2);
	tank_tracks_norm->Draw("E");
	c2->cd(3);
	tank_tracks_norm_H->Draw("E");
	c2->cd(4);
	tank_tracks_norm_V->Draw("E");
	c2->Modified();
	c2->Update();
	ratio_tracks_time->Divide(tank_tracks_time,out_tracks_time);
	c3->cd(1);
	tank_tracks_time->Draw("E");
	c3->cd(2);
	out_tracks_time->Draw("E");
	c3->cd(3);
	ratio_tracks_time->Draw("E");
	c3->cd(4);
	TH1D * ratio_tracks_rebinned = new TH1D("ratio_tracks_rebinned","ratio_track_rebinned",nentries/interval_length,interval_time);
	for(int i=0;i<(nentries/interval_length+1);i++){
		ratio_tracks_rebinned->SetBinContent(i,tank_tracks_norm->GetBinContent(i));
		ratio_tracks_rebinned->SetBinError(i,tank_tracks_norm->GetBinError(i));
	}
	ratio_tracks_rebinned->Draw("E");
	c3->Modified();
	c3->Update();
}
void Analyse::AbsorptionFluxMapNormTheo(double z, double bench_angle, TCanvas * c1, TCanvas * c2, TCanvas * c3, TCanvas * c4){
	long eventReconstructed = 0;
	long eventSuitable = 0;
	double chisquare_threshold = 100;

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	//gStyle->SetPalette(1);
	//double z_Pb = 1553;
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

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	int nbins = Sqrt(0.02*nentries);
	if(c1==0) c1 = new TCanvas("fluxMap","fluxMap");
	TH2D * fluxMapZ = new TH2D("fluxMapSignal","fluxMapSignal",nbins,x_min,x_max,nbins,x_min,x_max);
	fluxMapZ->SetStats(0);
	if(c2 == 0) c2 = new TCanvas("fluxMapNorm","fluxMapNorm");
	TH2D * fluxMapSigma = new TH2D("fluxMapSigma","fluxMapSigma",nbins,x_min,x_max,nbins,x_min,x_max);
	fluxMapSigma->SetStats(0);
	if(c3 == 0) c3 = new TCanvas("fluxMap_Sigma","fluxMap_Sigma");

	map<pair<pair<int,int>,pair<int,int> >,unsigned long> ray_class_n;

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			pair<pair<int,int>,pair<int,int> > current_extremal_det = it->get_extremal_det(this);
			map<pair<pair<int,int>,pair<int,int> >,unsigned long>::iterator current_ray_class = ray_class_n.find(current_extremal_det);
			if(current_ray_class!= ray_class_n.end()) (current_ray_class->second)++;
			else ray_class_n[current_extremal_det] = 1;
			fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%10000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			c1->cd();
			fluxMapZ->Draw("COLZ");
			c1->Modified();
			c1->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	
	c1->cd();
	fluxMapZ->Draw("COLZ");
	c1->Modified();
	c1->Update();
	
	TCanvas * c5 = new TCanvas("background_comp","background_comp");
	c5->Divide(CeilNint(Sqrt(ray_class_n.size())),CeilNint(Sqrt(ray_class_n.size())));

	TH2D * full_bkg = new TH2D("full_bkg","full_bkg",nbins,x_min,x_max,nbins,x_min,x_max);
	int comp_n = 1;
	cout << ray_class_n.size() << endl;
	for(map<pair<pair<int,int>,pair<int,int> >,unsigned long>::iterator class_it = ray_class_n.begin();class_it!=ray_class_n.end();++class_it){
		double size_X_Down = detectors[(class_it->first).first.first]->get_size();
		double size_X_Up = detectors[(class_it->first).first.second]->get_size();
		double size_Y_Down = detectors[(class_it->first).second.first]->get_size();
		double size_Y_Up = detectors[(class_it->first).second.second]->get_size();

		double offset_X_Down = detectors[(class_it->first).first.first]->get_offset();
		double offset_X_Up = detectors[(class_it->first).first.second]->get_offset();
		double offset_Y_Down = detectors[(class_it->first).second.first]->get_offset();
		double offset_Y_Up = detectors[(class_it->first).second.second]->get_offset();

		double z_X_Down = detectors[(class_it->first).first.first]->get_z();
		double z_X_Up = detectors[(class_it->first).first.second]->get_z();
		double z_Y_Down = detectors[(class_it->first).second.first]->get_z();
		double z_Y_Up = detectors[(class_it->first).second.second]->get_z();
		cout << "    " << z_X_Down << " | " << z_X_Up << " | " << z_Y_Down << " | " << z_Y_Up << " | " << (class_it->second) << endl;

		acceptanceFunction current_acceptance(offset_X_Up - size_X_Up/2,offset_X_Up + size_X_Up/2,offset_Y_Up - size_Y_Up/2,offset_Y_Up + size_Y_Up/2,offset_X_Down - size_X_Down/2,offset_X_Down + size_X_Down/2,offset_Y_Down - size_Y_Down/2,offset_Y_Down + size_Y_Down/2,z_X_Up,z_X_Down,z_Y_Up,z_Y_Down,bench_angle);
		TH2D * current_bkg = new TH2D(current_acceptance.plot_XY(nbins,x_min,x_max,nbins,x_min,x_max,z));
		full_bkg->Add(current_bkg,(class_it->second)/(current_bkg->Integral()));
		ostringstream current_name;
		current_name << "bkg_comp_" << comp_n << "_" << class_it->second;
		current_bkg->SetTitle(current_name.str().c_str());
		current_bkg->SetName(current_name.str().c_str());
		current_bkg->SetStats(0);
		c5->cd(comp_n);
		current_bkg->Draw("COLZ");
		comp_n++;
	}
	if(c4 == 0) c4 = new TCanvas("fluxMap_background","fluxMap_background");
	c4->cd();
	full_bkg->Draw("COLZ");
	c4->Modified();
	c4->Update();

	TH2D * copy = new TH2D(*fluxMapZ);
	copy->SetNameTitle("fluxMapDiff","fluxMapDiff");
	copy->SetStats(0);
	copy->Scale(1./copy->Integral());
	copy->Add(full_bkg,-1./full_bkg->Integral());
	copy->Scale(-1.);
	c2->cd();
	copy->Draw("COLZ");
	c2->Modified();
	c2->Update();
	c1->cd();
	fluxMapZ->Draw("COLZ");
	c1->Modified();
	c1->Update();
	double ratio = full_bkg->Integral()/fluxMapZ->Integral();
	for(int i=1;i<=nbins;i++){
		for(int j=1;j<=nbins;j++){
			int binN = fluxMapSigma->GetBin(i,j);
			double binContent = (fluxMapZ->GetBinContent(binN)*ratio - full_bkg->GetBinContent(binN))/Sqrt(full_bkg->GetBinContent(binN) + fluxMapZ->GetBinContent(binN)*ratio);
			fluxMapSigma->SetBinContent(binN,-binContent);
		}
	}
	c3->cd();
	fluxMapSigma->Draw("COLZ");
	c3->Modified();
	c3->Update();
}
void Analyse::AbsorptionFluxMapNormTheoAngle(double bench_angle, int mult, TCanvas * c1, TCanvas * c2, TCanvas * c3, TCanvas * c4){
	Display_Thread * MT_display = Display_Thread::get_instance();
	double chisquare_threshold = 100;

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	//gStyle->SetPalette(1);
	//double z_Pb = 1553;
	double z_max = numeric_limits<double>::min();
	double z_min = numeric_limits<double>::max();
	vector<double> z_det;
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		if((*det_it)->get_is_X() && ((*det_it)->get_perp_n()>0)) continue;
		double current_z = (*det_it)->get_z();
		if(current_z>z_max) z_max = current_z;
		if(current_z<z_min) z_min = current_z;
		z_det.push_back(current_z);
	}
	double phi_max = ATan(Tomography::get_instance()->get_XY_size()/(z_max - z_min));
	phi_max += 0.1*phi_max;

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	int nbins = Sqrt(0.02*nentries);
	if(c1==0) c1 = new TCanvas("fluxMap","fluxMap");
	TH2D * fluxMapZ = new TH2D("fluxMapSignal","fluxMapSignal",nbins,-phi_max,phi_max,nbins,bench_angle-phi_max,bench_angle+phi_max);
	fluxMapZ->SetStats(0);
	if(c2 == 0) c2 = new TCanvas("fluxMapNorm","fluxMapNorm");
	TH2D * fluxMapSigma = new TH2D("fluxMapSigma","fluxMapSigma",nbins,-phi_max,phi_max,nbins,bench_angle-phi_max,bench_angle+phi_max);
	fluxMapSigma->SetStats(0);
	if(c3 == 0) c3 = new TCanvas("fluxMap_Sigma","fluxMap_Sigma");
	FreeSkyFunction acceptanceEstimation(-Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,-Tomography::get_instance()->get_XY_size()/2.,Tomography::get_instance()->get_XY_size()/2.,z_det, bench_angle);
	TH2D * background = new TH2D(acceptanceEstimation.plot_PhiTheta(nbins,-phi_max,phi_max,nbins,bench_angle-phi_max,bench_angle+phi_max,mult));
	if(c4 == 0) c4 = new TCanvas("fluxMap_background","fluxMap_background");
	c4->cd();
	background->Draw("COLZ");
	c4->Modified();
	c4->Update();

	Buffer_Task<ray_data> * ray_list = new Buffer_Task<ray_data>();
	Input_Task * to_do = new Read_Analyse_Task(nentries,this,this, new Tracking_Abs_Task(this, ray_list));
	vector<Thread*> threads;
	threads.push_back(new Reader_Thread(to_do));
	(threads.back())->start();
	const unsigned short n_thread = (Tomography::get_instance()->get_thread_number() > threads.size()) ? (Tomography::get_instance()->get_thread_number() - threads.size()) : 1;
	*MT_display << "1 | " << n_thread << "\n";
	for(unsigned short i=0;i<n_thread;i++){
		threads.push_back(new Worker_Thread());
		(threads.back())->start();
	}
	usleep(1000);
	MT_display->start_count();
	//MT_display << Tomography::get_instance()->init_count() << "|" << setw(7) << "tracks\n";
	bool has_working_thread = true;
	unsigned long jentry = 0;
	while((has_working_thread || ray_list->can_fetch_data()) && Tomography::get_instance()->get_can_continue()){
		has_working_thread = false;
		for(unsigned short i=0;i<threads.size();i++){
			if(threads[i]->is_working()){
				has_working_thread = true;
				break;
			}
		}
		if(!(ray_list->can_fetch_data())){
			usleep(1000);
			continue;
		}
		ray_data * current_rays = ray_list->fetch_data();
		vector<Ray>::iterator ray_it = (current_rays->rays).begin();
		while(ray_it != (current_rays->rays).end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				fluxMapZ->Fill(ATan(ray_it->get_slope_X()),bench_angle+ATan(ray_it->get_slope_Y()));
			}
			++ray_it;
		}
		delete current_rays;
		jentry++;
		if(jentry%10000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			TH2D * copy = new TH2D(*fluxMapZ);
			copy->SetNameTitle("fluxMapDiff","fluxMapDiff");
			copy->SetStats(0);
			copy->Scale(1./copy->Integral());
			copy->Add(background,-1./background->Integral());
			copy->Scale(-1.);
			c2->cd();
			copy->Draw("COLZ");
			c2->Modified();
			c2->Update();
			delete copy;
			c1->cd();
			fluxMapZ->Draw("COLZ");
			c1->Modified();
			c1->Update();
			double ratio = background->Integral()/fluxMapZ->Integral();
			for(int i=1;i<=nbins;i++){
				for(int j=1;j<=nbins;j++){
					int binN = fluxMapSigma->GetBin(i,j);
					double binContent = (fluxMapZ->GetBinContent(binN)*ratio - background->GetBinContent(binN))/Sqrt(background->GetBinContent(binN) + fluxMapZ->GetBinContent(binN)*ratio);
					fluxMapSigma->SetBinContent(binN,-binContent);
				}
			}
			c3->cd();
			fluxMapSigma->Draw("COLZ");
			c3->Modified();
			c3->Update();
		}
		//*MT_display << jentry << " | " << Tomography::get_instance()->get_can_continue() << " | " << has_working_thread << " | " << ray_list->can_fetch_data() << "\n";
	}
	for(unsigned short i=0;i<threads.size();i++){
		threads[i]->stop();
		delete threads[i];
	}

	//MT_display << "\r" << Tomography::get_instance()->print_count() << "|" << setw(7) << eventReconstructed << "\n";
	MT_display->stop();

	TH2D * copy = new TH2D(*fluxMapZ);
	copy->SetNameTitle("fluxMapDiff","fluxMapDiff");
	copy->SetStats(0);
	copy->Scale(1./copy->Integral());
	copy->Add(background,-1./background->Integral());
	copy->Scale(-1.);
	c2->cd();
	copy->Draw("COLZ");
	c2->Modified();
	c2->Update();
	c1->cd();
	fluxMapZ->Draw("COLZ");
	c1->Modified();
	c1->Update();
	double ratio = background->Integral()/fluxMapZ->Integral();
	for(int i=1;i<=nbins;i++){
		for(int j=1;j<=nbins;j++){
			int binN = fluxMapSigma->GetBin(i,j);
			double binContent = (fluxMapZ->GetBinContent(binN)*ratio - background->GetBinContent(binN))/Sqrt(background->GetBinContent(binN) + fluxMapZ->GetBinContent(binN)*ratio);
			fluxMapSigma->SetBinContent(binN,-binContent);
		}
	}
	c3->cd();
	fluxMapSigma->Draw("COLZ");
	c3->Modified();
	c3->Update();
}
void Analyse::AbsorptionFluxMapNorm(double z,TH2D * background, int nbins, TCanvas * c1, TCanvas * c2, TCanvas * c3){
	long eventReconstructed = 0;
	long eventSuitable = 0;

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	//gStyle->SetPalette(1);
	//double z_Pb = 1553;
	if(c1==0) c1 = new TCanvas("fluxMap","fluxMap");
	TH2D * fluxMapZ = new TH2D("fluxMapSignal","fluxMapSignal",background->GetNbinsX(),background->GetXaxis()->GetXmin(),background->GetXaxis()->GetXmax(),background->GetNbinsY(),background->GetYaxis()->GetXmin(),background->GetYaxis()->GetXmax());
	fluxMapZ->SetStats(0);
	if(c2 == 0) c2 = new TCanvas("fluxMapNorm","fluxMapNorm");
	TH2D * fluxMapSigma = new TH2D("fluxMapSigma","fluxMapSigma",background->GetNbinsX(),background->GetXaxis()->GetXmin(),background->GetXaxis()->GetXmax(),background->GetNbinsY(),background->GetYaxis()->GetXmin(),background->GetYaxis()->GetXmax());
	fluxMapSigma->SetStats(0);
	if(c3 == 0) c3 = new TCanvas("fluxMap_Sigma","fluxMap_Sigma");

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		//eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
				fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
				eventReconstructed++;
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%10000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			TH2D * copy = new TH2D(*fluxMapZ);
			copy->SetNameTitle("fluxMapDiff","fluxMapDiff");
			copy->SetStats(0);
			copy->Scale(1./copy->Integral());
			copy->Add(background,-1./background->Integral());
			copy->Scale(-1.);
			c2->cd();
			copy->Draw("COLZ");
			c2->Modified();
			c2->Update();
			delete copy;
			c1->cd();
			fluxMapZ->Draw("COLZ");
			c1->Modified();
			c1->Update();
			double ratio = background->Integral()/fluxMapZ->Integral();
			for(int i=1;i<=nbins;i++){
				for(int j=1;j<=nbins;j++){
					int binN = fluxMapSigma->GetBin(i,j);
					double binContent = (fluxMapZ->GetBinContent(binN)*ratio - background->GetBinContent(binN))/Sqrt(background->GetBinContent(binN) + fluxMapZ->GetBinContent(binN)*ratio);
					fluxMapSigma->SetBinContent(binN,-binContent);
				}
			}
			c3->cd();
			fluxMapSigma->Draw("COLZ");
			c3->Modified();
			c3->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	TH2D * copy = new TH2D(*fluxMapZ);
	copy->SetNameTitle("fluxMapDiff","fluxMapDiff");
	copy->SetStats(0);
	copy->Scale(1./copy->Integral());
	copy->Add(background,-1./background->Integral());
	copy->Scale(-1.);
	c2->cd();
	copy->Draw("COLZ");
	c2->Modified();
	c2->Update();
	c1->cd();
	fluxMapZ->Draw("COLZ");
	c1->Modified();
	c1->Update();
	double ratio = background->Integral()/fluxMapZ->Integral();
	for(int i=1;i<=nbins;i++){
		for(int j=1;j<=nbins;j++){
			int binN = fluxMapSigma->GetBin(i,j);
			double binContent = (fluxMapZ->GetBinContent(binN)*ratio - background->GetBinContent(binN))/Sqrt(background->GetBinContent(binN) + fluxMapZ->GetBinContent(binN)*ratio);
			fluxMapSigma->SetBinContent(binN,-binContent);
		}
	}
	c3->cd();
	fluxMapSigma->Draw("COLZ");
	c3->Modified();
	c3->Update();
}

void Analyse::StoreRayPairs(string outFileName){

	double z_Up = get_z_Up();
	double z_Down = get_z_Down();

	TCanvas * c1 = new TCanvas();
	c1->Divide(2,2);
	TH1D * thetaXUp = new TH1D("thetaXUp","thetaXUp",50,-2,2);
	thetaXUp->SetLineColor(2);
	TH1D * thetaYUp = new TH1D("thetaYUp","thetaYUp",50,-2,2);
	thetaYUp->SetLineColor(4);
	TH1D * thetaXDown = new TH1D("thetaXDown","thetaXDown",50,-2,2);
	thetaXDown->SetLineColor(2);
	TH1D * thetaYDown = new TH1D("thetaYDown","thetaYDown",50,-2,2);
	thetaYDown->SetLineColor(4);

	TH2D * XY_correlation = new TH2D("correlation","correlation",70,-6*Tomography::get_instance()->get_XY_size()/10.,6*Tomography::get_instance()->get_XY_size()/10.,70,-6*Tomography::get_instance()->get_XY_size()/10.,6*Tomography::get_instance()->get_XY_size()/10.);

	TH1D * doca = new TH1D("doca","doca",100,0,Tomography::get_instance()->get_XY_size());

	TCanvas * c2 = new TCanvas();
	c2->Divide(3);
	TH1D * pocaX = new TH1D("pocaX","pocaX",100,-Tomography::get_instance()->get_XY_size(),Tomography::get_instance()->get_XY_size());
	TH1D * pocaY = new TH1D("pocaY","pocaY",100,-Tomography::get_instance()->get_XY_size(),Tomography::get_instance()->get_XY_size());
	
	TCanvas * c3 = new TCanvas();
	c3->Divide(3,2);
	TH1D * chi2_h = new TH1D("chisquares","chisquares",100,0,1);
	c3->GetPad(1)->SetLogy();
	TH1D * sizes = new TH1D("sizes","sizes",10,0,10);
	TH1D * suitable_h = new TH1D("suitable","suitable",10,0,10);
	suitable_h->SetLineColor(2);
	c3->GetPad(2)->SetLogy();
	TH1D * nclus_xh_h = new TH1D("nclus_xh","nclus_xh",20,0,20);
	TH1D * nclus_yh_h = new TH1D("nclus_yh","nclus_yh",20,0,20);
	TH1D * nclus_xb_h = new TH1D("nclus_xb","nclus_xb",20,0,20);
	TH1D * nclus_yb_h = new TH1D("nclus_yb","nclus_yb",20,0,20);
	
	double z_Up_tot = numeric_limits<double>::min();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_z()>z_Up_tot){
			z_Up_tot = (*it)->get_z();
		}
	}
	double z_Down_tot = numeric_limits<double>::max();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_z()<z_Down){
			z_Down_tot = (*it)->get_z();
		}
	}
	TH1D * pocaZ = new TH1D("pocaZ","pocaZ",100,z_Down_tot - (z_Up_tot-z_Down_tot),z_Up_tot + (z_Up_tot-z_Down_tot));
	TLine * downLim = new TLine(z_Down_tot,0,z_Down_tot,1);
	TLine * upLim = new TLine(z_Up_tot,0,z_Up_tot,1);
	downLim->SetLineStyle(2);
	upLim->SetLineStyle(2);
	downLim->SetLineColor(2);
	upLim->SetLineColor(2);
	long eventReconstructed = 0;
	long eventSuitable = 0;

	TFile * outFile = new TFile(outFileName.c_str(),"RECREATE");
	TTree * outTree = new TTree("T","T");

	double X_Up = 0;
	double Y_Up = 0;
	double theta_X_Up = 0;
	double theta_Y_Up = 0;
	double X_Down = 0;
	double Y_Down = 0;
	double theta_X_Down = 0;
	double theta_Y_Down = 0;

	outTree->Branch("X_Up", &X_Up, "X_Up/D");
	outTree->Branch("Y_Up", &Y_Up, "Y_Up/D");
	outTree->Branch("theta_X_Up", &theta_X_Up, "theta_X_Up/D");
	outTree->Branch("theta_Y_Up", &theta_Y_Up, "theta_Y_Up/D");
	outTree->Branch("X_Down", &X_Down, "X_Down/D");
	outTree->Branch("Y_Down", &Y_Down, "Y_Down/D");
	outTree->Branch("theta_X_Down", &theta_X_Down, "theta_X_Down/D");
	outTree->Branch("theta_Y_Down", &theta_Y_Down, "theta_Y_Down/D");
	outTree->Branch("Z_Up", &z_Up, "Z_Up/D");
	outTree->Branch("Z_Down", &z_Down, "Z_Down/D");

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		//for perfect simu
		/*
		bool multi = false;
		for(int i =0;i<MG_N;i++){
			if(MG_NClus[i]>1) multi = true;
		}
		if(multi) continue;
		*/

		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		currentCBEvent->createPairs();
		eventSuitable+=currentCBEvent->get_clus_N()/(get_det_N_tot());
		bool is_suitable = true;	
		for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
			if(currentCBEvent->get_clus_N_by_det(*it)<1) is_suitable = false;
			if(Tomography::get_instance()->get_is_up((*it)->get_layer())){
				if((*it)->get_is_X()){
					nclus_xh_h->Fill(currentCBEvent->get_clus_N_by_det(*it));
				}
				else{	
					nclus_yh_h->Fill(currentCBEvent->get_clus_N_by_det(*it));
				}
			}
			else{
				if((*it)->get_is_X()){
					nclus_xb_h->Fill(currentCBEvent->get_clus_N_by_det(*it));
				}
				else{	
					nclus_yb_h->Fill(currentCBEvent->get_clus_N_by_det(*it));
				}
			}
		}
		if(is_suitable) suitable_h->Fill(1);
		else suitable_h->Fill(0);
		eventReconstructed+=currentCBEvent->get_rayPairs_N();
		sizes->Fill(currentCBEvent->get_rayPairs_N());
		for(unsigned int i=0;i<currentCBEvent->get_rayPairs_N();i++){
			RayPair currentRayPair = currentCBEvent->get_rayPair(i);
			double current_chi = 0;
			current_chi += currentRayPair.downRay.get_chiSquare_X()*currentRayPair.downRay.get_chiSquare_X();
			current_chi += currentRayPair.downRay.get_chiSquare_Y()*currentRayPair.downRay.get_chiSquare_Y();
			current_chi += currentRayPair.upRay.get_chiSquare_X()*currentRayPair.upRay.get_chiSquare_X();
			current_chi += currentRayPair.upRay.get_chiSquare_Y()*currentRayPair.upRay.get_chiSquare_Y();
			current_chi = Sqrt(current_chi);
			chi2_h->Fill(current_chi);
			doca->Fill(currentRayPair.get_doca());
			thetaXUp->Fill(currentRayPair.upRay.get_slope_X());
			thetaYUp->Fill(currentRayPair.upRay.get_slope_Y());
			thetaXDown->Fill(currentRayPair.downRay.get_slope_X());
			thetaYDown->Fill(currentRayPair.downRay.get_slope_Y());
			X_Up = currentRayPair.get_x_up(z_Up);
			Y_Up = currentRayPair.get_y_up(z_Up);
			theta_X_Up = currentRayPair.get_theta_x_up();
			theta_Y_Up = currentRayPair.get_theta_y_up();
			X_Down = currentRayPair.get_x_down(z_Down);
			Y_Down = currentRayPair.get_y_down(z_Down);
			theta_X_Down = currentRayPair.get_theta_x_down();
			theta_Y_Down = currentRayPair.get_theta_y_down();
			Point PoCA = currentRayPair.get_PoCA();
			/*if(PoCA.get_Z()>600 && PoCA.get_Z()<800)*/ XY_correlation->Fill(PoCA.get_X(),PoCA.get_Y());
			pocaX->Fill(PoCA.get_X());
			pocaY->Fill(PoCA.get_Y());
			pocaZ->Fill(PoCA.get_Z());
			outTree->Fill();

		}
		delete currentCBEvent;
		/*if(jentry%500 == 0)*/ cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			c1->cd(1);
			thetaXUp->Draw();
			thetaYUp->Draw("SAME");
			c1->cd(2);
			thetaXDown->Draw();
			thetaYDown->Draw("SAME");
			c1->cd(3);
			doca->Draw();
			c1->cd(4);
			XY_correlation->Draw("COLZ");
			c1->Modified();
			c1->Update();
			c2->cd(1);
			pocaX->Draw();
			c2->cd(2);
			pocaY->Draw();
			c2->cd(3);
			pocaZ->Draw();
			upLim->SetY2(pocaZ->GetMaximum());
			downLim->SetY2(pocaZ->GetMaximum());
			upLim->Draw();
			downLim->Draw();
			c2->Modified();
			c2->Update();
			c3->cd(1);
			chi2_h->Draw();
			c3->cd(2);
			sizes->Draw();
			suitable_h->Draw("SAME");
			c3->cd(3);
			nclus_xh_h->Draw();
			c3->cd(4);
			nclus_yh_h->Draw();
			c3->cd(5);
			nclus_xb_h->Draw();
			c3->cd(6);
			nclus_yb_h->Draw();
			c3->Modified();
			c3->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	outFile->cd();
	outTree->Write();
	outFile->Close();
	c1->cd(1);
	thetaXUp->Draw();
	thetaYUp->Draw("SAME");
	c1->cd(2);
	thetaXDown->Draw();
	thetaYDown->Draw("SAME");
	c1->cd(3);
	doca->Draw();
	c1->cd(4);
	XY_correlation->Draw("COLZ");
	c1->Modified();
	c1->Update();
	c2->cd(1);
	pocaX->Draw();
	c2->cd(2);
	pocaY->Draw();
	c2->cd(3);
	pocaZ->Draw();
	upLim->SetY2(pocaZ->GetMaximum());
	downLim->SetY2(pocaZ->GetMaximum());
	upLim->Draw();
	downLim->Draw();
	c2->Modified();
	c2->Update();
	c3->cd(1);
	chi2_h->Draw();
	c3->cd(2);
	sizes->Draw();
	suitable_h->Draw("SAME");
	c3->cd(3);
	nclus_xh_h->Draw();
	c3->cd(4);
	nclus_yh_h->Draw();
	c3->cd(5);
	nclus_xb_h->Draw();
	c3->cd(6);
	nclus_yb_h->Draw();
	c3->Modified();
	c3->Update();
}

/*void Analyse::MultiGenDebug(int i){
	long eventReconstructed = 0;
	long eventSuitable = 0;

	TCanvas * c1 = new TCanvas("MGPos","MGPos");
	TH1D * mgPos = new TH1D("mgPos","mgPos",1024,0,1023);

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		for(vector<Event*>::iterator it=(currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if((*it)->get_type() == Tomography::MG && (*it)->get_n_in_tree() == i){
				MG_Event currentMGEvent(*dynamic_cast<MG_Event*>(*it));
				for(vector<MG_Cluster>::iterator jt = currentMGEvent.clusters.begin();jt!=currentMGEvent.clusters.end();++jt){
					mgPos->Fill(jt->pos);
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			c1->cd();
			mgPos->Draw();
			c1->Modified();
			c1->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	c1->cd();
	mgPos->Draw();
	c1->Modified();
	c1->Update();
}*/

void Analyse::bugtest(){
	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	int limit = 10;
	if(limit<nentries) nentries = limit;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		cout << evn << " : " << endl;
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,-1);
		for(vector<Event*>::iterator it = (CBEvent->events).begin();it!=(CBEvent->events).end();++it){
			if((*it)->get_type() == Tomography::MG){
				vector<Cluster*> current_clusters = (*it)->get_clusters();
				if(current_clusters.size()>0) cout << setw(10) << (current_clusters.front())->get_pos() << " | ";
				else cout << "No Cluster" << " | ";
				for(vector<Cluster*>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
					delete *clus_it;
				}
			}
		}
		cout << endl;
		delete CBEvent;
		for(int i=0;i<det_N[Tomography::MG]-1;i++){
			cout << setw(10) << reinterpret_cast<Double_t(*)[MG_Detector::MaxNClus]>(ClusPos[Tomography::MG])[i][0] << " | ";
		}
		cout << setw(10) << reinterpret_cast<Double_t(*)[MG_Detector::MaxNClus]>(ClusPos[Tomography::MG])[det_N[Tomography::MG]-1][0] << endl;
	}
}

void Analyse::CalcStripResponseFunction(int bin_nb){

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	gStyle->SetOptFit(0111);

	int det_num = get_det_N_tot();
	cout << setw(20) << "detector proccessed" << "|" << setw(20) << "event processed" << endl;
	TProfile * SRH;
	//TGraph * SRH2D[det_N];
	TH2D * SRH2D;
	TCanvas * c;
	TCanvas * d;
	TF1 * SRF;

	TCanvas ** c_coord = new TCanvas*[bin_nb];
	TProfile ** SRH_coord = new TProfile*[bin_nb];
	TF1 ** SRF_coord = new TF1*[bin_nb];
	TCanvas * c_param;
	TGraphErrors * gauss_width_graph;
	TGraphErrors * lorentz_width_graph;
	TGraphErrors * offset_graph;

	double chisquare_threshold = 10;

	TFile * signal_file = new TFile(signal_file_name.c_str(),"READ");
	TTree * signal_tree = (TTree*)(signal_file->Get("T"));

	if(fChain->GetEntriesFast() != signal_tree->GetEntriesFast()){
		cout << "total number of event in signal and analyse tree does not match" << endl;
		return;
	}

	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	float StripAmpl_MG_corr[det_N[Tomography::MG]][MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()];
	float StripAmpl_CM_corr[det_N[Tomography::CM]][CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()];
	int signal_evn;
	signal_tree->SetBranchAddress("Nevent",&signal_evn);
	if(det_N[Tomography::CM]>0) signal_tree->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
	if(det_N[Tomography::MG]>0) signal_tree->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);

	int det_x_n = 0;
	int det_y_n = 0;
	int nref_nb = 0;
	for(int j=0;j<det_num;j++){
		if(detectors[j]->get_is_X()) det_x_n++;
		else det_y_n++;
		if(!(detectors[j]->get_is_ref())) nref_nb++;
	}
	if(nref_nb>1){
		cout << "there is more than 1 non ref det, exiting" << endl;
		return;
	}

	for(int i=0;i<det_num;i++){
		if(detectors[i]->get_is_ref()) continue;
		int det_in_nref_dir = (detectors[i]->get_is_X()) ? det_x_n : det_y_n;

		ostringstream c_name;
		c_name << "c_" << i;
		ostringstream d_name;
		d_name << "d_" << i;
		ostringstream SRH_name;
		SRH_name << "SRH_" << i;
		ostringstream SRH2D_name;
		SRH2D_name << "SRH2D_" << i;
		ostringstream SRF_name;
		SRF_name << "SRF_" << i;

		double limit = 6;
		int bin_n = 200;
		if(detectors[i]->get_type() == Tomography::MG && detectors[i]->get_is_X() == false) limit *= 2.;
		c = new TCanvas(c_name.str().c_str(),c_name.str().c_str());
		d = new TCanvas(d_name.str().c_str(),d_name.str().c_str());
		SRH = new TProfile(SRH_name.str().c_str(),SRH_name.str().c_str(),bin_n,-limit,limit,0,1000);
		//SRH2D[i] = new TGraph();
		//int graph_point_nb = 0;
		SRH2D = new TH2D(SRH2D_name.str().c_str(),SRH2D_name.str().c_str(),200,-limit,limit,100,0,10);
		
		SRF = new TF1(SRF_name.str().c_str(),"[3] + ((1-[3])*(exp(-4*log(2)*(1-[0])*(x-[4])*(x-[4])/([1]*[1])))/(1+(4*[0]*(x-[4])*(x-[4])/([2]*[2]))))",-limit,limit);
		//([3] + (exp(-4*log(2)*(1-[0])*(x-[4])*(x-[4])/([1]*[1]))/(1+(4*[0]*(x-[4])*(x-[4])/([2]*[2]))))/(1+[3]))
		//[3] + ((1-[3])*(exp(-4*log(2)*(1-[0])*(x-[4])*(x-[4])/([1]*[1])))/(1+(4*[0]*(x-[4])*(x-[4])/([2]*[2]))))
		SRF->SetParameters(0.5,0.5,0.5,0.1,0);
		SRF->SetParLimits(0,0,1);
		SRF->SetParLimits(1,0,10000);
		SRF->SetParLimits(2,0,10000);
		SRF->SetParLimits(3,0,0.5);
		SRF->SetParLimits(4,-4,4);

		for(int j=0;j<bin_nb;j++){
			ostringstream c_graph_name;
			c_graph_name << "c_" << j*1024/bin_nb << "_"<< (j+1)*1024/bin_nb;
			ostringstream SRH_graph_name;
			SRH_graph_name << "SRH_" << j*1024/bin_nb << "_"<< (j+1)*1024/bin_nb;
			ostringstream SRF_graph_name;
			SRF_graph_name << "SRF_" << j*1024/bin_nb << "_"<< (j+1)*1024/bin_nb;
			c_coord[j] = new TCanvas(c_graph_name.str().c_str(),c_graph_name.str().c_str());
			SRH_coord[j] = new TProfile(SRH_graph_name.str().c_str(),SRH_graph_name.str().c_str(),bin_n,-limit,limit,0,1000);
			SRF_coord[j] = new TF1(SRF_graph_name.str().c_str(),"[3] + ((1-[3])*(exp(-4*log(2)*(1-[0])*(x-[4])*(x-[4])/([1]*[1])))/(1+(4*[0]*(x-[4])*(x-[4])/([2]*[2]))))",-limit,limit);
			SRF_coord[j]->SetParameters(0.5,0.5,0.5,0.1,0);
			SRF_coord[j]->SetParLimits(0,0,1);
			SRF_coord[j]->SetParLimits(1,0,10000);
			SRF_coord[j]->SetParLimits(2,0,10000);
			SRF_coord[j]->SetParLimits(3,0,0.5);
			SRF_coord[j]->SetParLimits(4,-4,4);
		}

		gauss_width_graph = new TGraphErrors();
		lorentz_width_graph = new TGraphErrors();
		offset_graph = new TGraphErrors();

		if (fChain == 0) return;
		for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			fChain->GetEntry(jentry);
			signal_tree->LoadTree(jentry);
			signal_tree->GetEntry(jentry);
			if(signal_evn!=(evn+0)){
				cout << "event numbers are different in analyse and signal trees" << endl;
				return;
			}
			CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
			vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);

			for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
				if(!((*it)->get_is_ref())){
					if((*it)->get_type() == Tomography::MG){
						vector<Cluster*> current_clusters = (*it)->get_clusters();
						for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){

							double chiSquare_in_nref_dir = (detectors[i]->get_is_X()) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
							int clus_in_nref_dir = (detectors[i]->get_is_X()) ? jt->get_clus_x_n() : jt->get_clus_y_n();
							if(clus_in_nref_dir<(det_in_nref_dir-1)) continue;
							if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;

							if((jt->get_chiSquare_X()+jt->get_chiSquare_Y()) > chisquare_threshold) continue;
							double residu = numeric_limits<double>::max();
							vector<Cluster*>::iterator matching_cluster = current_clusters.end();
							for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
								(*kt)->set_perp_pos_mm(*jt);
								double current_residu = jt->get_residu_ref(*kt);
								if(current_residu<residu){
									residu = current_residu;
									matching_cluster = kt;
								}
							}
							if(matching_cluster == current_clusters.end()) continue;
							
							double normalization = (*matching_cluster)->get_ampl()/(*matching_cluster)->get_size();
							double matching_position = ((*matching_cluster)->get_is_X()) ? jt->eval_X((*it)->get_z()) : jt->eval_Y((*it)->get_z());
							double matching_position_perp = ((*matching_cluster)->get_is_X()) ? jt->eval_Y((*it)->get_z()) : jt->eval_X((*it)->get_z());
							for(int strip_nb=0;strip_nb<1024;strip_nb++){
								int channel = MG_Detector::StripToChannel_a[strip_nb];
								if(Abs(residu)<50.){
									SRH->Fill(matching_position - (*matching_cluster)->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::get_instance()->get_Nsample()))/normalization);
									if(bin_nb>0) SRH_coord[Min(Max(FloorNint((bin_nb/2.)+matching_position_perp*bin_nb/Tomography::get_instance()->get_XY_size()),0),bin_nb-1)]->Fill(matching_position - (*matching_cluster)->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::get_instance()->get_Nsample()))/normalization);
								}
								//SRH2D[i]->SetPoint(graph_point_nb, matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::get_instance()->get_Nsample()))/normalization);
								//graph_point_nb++;
								SRH2D->Fill(matching_position - (*matching_cluster)->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::get_instance()->get_Nsample()))/normalization);
							}
							delete *matching_cluster;
							current_clusters.erase(matching_cluster);
						}
						for(vector<Cluster*>::iterator clus_it=current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
							delete *clus_it;
						}
					}
				}
			}

			delete currentCBEvent;

			if(jentry%500 == 0) cout << "\r" << setw(20) << i << "|" << setw(20) << jentry << flush;
			if(jentry%10000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
				c->cd();
				SRH->Draw();
				c->Modified();
				c->Update();
				d->cd();
				//SRH2D[i]->GetXaxis()->SetLimits(-limit,limit);
				SRH2D->Draw("colz");
				d->Modified();
				d->Update();
				for(int j=0;j<bin_nb;j++){
					c_coord[j]->cd();
					SRH_coord[j]->Draw();
					c_coord[j]->Modified();
					c_coord[j]->Update();
				}
			}
		}
		cout << "\r" << setw(20) << i << "|" << setw(20) << nentries << endl;
		/*
		double scaling_factor = 0;
		for(int n=0;n<bin_n;n++){
			if(SRH[i]->GetBinCenter(n)<-2) continue;
			if(SRH[i]->GetBinCenter(n)>2) continue;
			if(SRH[i]->GetBinContent(n)>scaling_factor) scaling_factor = SRH[i]->GetBinContent(n);
		}
		SRH[i]->Scale(1./scaling_factor);
		*/
		//TF1 * scaling_function = new TF1("scaling_function","[0]+[1]*x+[2]*x*x",-1,1);
		TF1 * scaling_function = new TF1("scaling_function","[0]*exp(-[1]*(x-[2])*(x-[2]))",-1,1);
		SRH->Fit(scaling_function,"QNR");
		//double scaling_factor = scaling_function->GetParameter(0) - 0.25*scaling_function->GetParameter(1)*scaling_function->GetParameter(1)/scaling_function->GetParameter(2);
		double scaling_factor = scaling_function->GetParameter(0);
		SRH->Scale(1./scaling_factor);
		delete scaling_function;
		
		SRH->Fit(SRF,"QN");
		c->cd();
		SRH->Draw();
		SRF->Draw("SAME");
		c->Modified();
		c->Update();
		d->cd();
		SRH2D->Draw("colz");
		d->Modified();
		d->Update();
		cout << "global fit " << endl;
		cout << setw(5) << " " << setw(30) << left << "offset : " << setw(10) << right << SRF->GetParameter(3) << endl;
		cout << setw(5) << " " << setw(30) << left << "misalignement : " << setw(10) << right << SRF->GetParameter(4) << endl;
		cout << setw(5) << " " << setw(30) << left << "gaussian width : " << setw(10) << right << SRF->GetParameter(1) << endl;
		cout << setw(5) << " " << setw(30) << left << "lorentzian width : " << setw(10) << right << SRF->GetParameter(2) << endl;
		cout << setw(5) << " " << setw(30) << left << "gauss/lorentz ratio : " << setw(10) << right << SRF->GetParameter(0) << endl;
		ostringstream fitQuality;
		fitQuality << SRF->GetChisquare() << "/" << SRF->GetNDF();
		cout << setw(5) << " " << setw(30) << left << "Chi2/NDF : " << setw(10) << right << fitQuality.str() << endl;

		for(int i_coord=0;i_coord<bin_nb;i_coord++){
			scaling_function = new TF1("scaling_function","[0]*exp(-[1]*(x-[2])*(x-[2]))",-1,1);
			SRH_coord[i_coord]->Fit(scaling_function,"QNR");
			SRH_coord[i_coord]->Scale(1./scaling_function->GetParameter(0));
			delete scaling_function;
			SRH_coord[i_coord]->Fit(SRF_coord[i_coord],"QN");
			gauss_width_graph->SetPoint(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.))+(i_coord*Tomography::get_instance()->get_XY_size()/bin_nb),SRF_coord[i_coord]->GetParameter(1));
			gauss_width_graph->SetPointError(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(1));
			lorentz_width_graph->SetPoint(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.))+(i_coord*Tomography::get_instance()->get_XY_size()/bin_nb),SRF_coord[i_coord]->GetParameter(2));
			lorentz_width_graph->SetPointError(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(2));
			offset_graph->SetPoint(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.))+(i_coord*Tomography::get_instance()->get_XY_size()/bin_nb),SRF_coord[i_coord]->GetParameter(3));
			offset_graph->SetPointError(i_coord,(Tomography::get_instance()->get_XY_size()/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(3));
			c_coord[i_coord]->cd();
			SRH_coord[i_coord]->Draw();
			SRF_coord[i_coord]->Draw("SAME");
			c_coord[i_coord]->Modified();
			c_coord[i_coord]->Update();
			cout << i_coord+1 << " fit " << endl;
			cout << setw(5) << " " << setw(30) << left << "offset : " << setw(10) << right << SRF_coord[i_coord]->GetParameter(3) << endl;
			cout << setw(5) << " " << setw(30) << left << "misalignement : " << setw(10) << right << SRF_coord[i_coord]->GetParameter(4) << endl;
			cout << setw(5) << " " << setw(30) << left << "gaussian width : " << setw(10) << right << SRF_coord[i_coord]->GetParameter(1) << endl;
			cout << setw(5) << " " << setw(30) << left << "lorentzian width : " << setw(10) << right << SRF_coord[i_coord]->GetParameter(2) << endl;
			cout << setw(5) << " " << setw(30) << left << "gauss/lorentz ratio : " << setw(10) << right << SRF_coord[i_coord]->GetParameter(0) << endl;
			ostringstream fitQuality_coord;
			fitQuality_coord << SRF_coord[i_coord]->GetChisquare() << "/" << SRF_coord[i_coord]->GetNDF();
			cout << setw(5) << " " << setw(30) << left << "Chi2/NDF : " << setw(10) << right << fitQuality_coord.str() << endl;
		}
		if(bin_nb>0){
			c_param = new TCanvas("param_evolution","param_evolution");
			c_param->Divide(3);
			gauss_width_graph->SetTitle("gauss_width");
			c_param->cd(1);
			gauss_width_graph->Draw("AP");
			lorentz_width_graph->SetTitle("lorentz_width");
			c_param->cd(2);
			lorentz_width_graph->Draw("AP");
			offset_graph->SetTitle("offset");
			c_param->cd(3);
			offset_graph->Draw("AP");
			c_param->Modified();
			c_param->Update();
		}
	}
	//cout << "\r" << setw(20) << det_N << "|" << setw(20) << nentries << endl;
}

void Analyse::EventDisplay(long event_nb, TCanvas * c1){
	long nentries = fChain->GetEntriesFast();
	if(event_nb<0 || event_nb>nentries){
		cout << "invalid event number" << endl;
		return;
	}
	TFile * signal_file = new TFile(signal_file_name.c_str(),"READ");
	TTree * signal_tree;// = (TTree*)(signal_file->Get("T"));
	Tsignal_R * signalT;// = new Tsignal_R(signal_tree,det_N);
	if(signal_file->IsOpen()){
		signal_tree = (TTree*)(signal_file->Get("T"));
		signalT = new Tsignal_R(signal_tree,det_N);
		if(nentries != signal_tree->GetEntriesFast()){
			cout << "total number of event in signal and analyse tree does not match" << endl;
			return;
		}
	}
	LoadTree(event_nb);
	GetEntry(event_nb);
	CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,-1);
	if(signal_file->IsOpen()){
		signalT->LoadTree(event_nb);
		signalT->GetEntry(event_nb);
		for(vector<Event*>::iterator ev_it = (CBEvent->events).begin();ev_it!=(CBEvent->events).end();++ev_it){
			(*ev_it)->set_strip_ampl(signalT->get_ampl<double>((*ev_it)->get_type(),(*ev_it)->get_n_in_tree()));
		}
		delete signal_file;
		delete signalT;
	}
	CBEvent->EventDisplay(c1);
	delete CBEvent;
	/*
	long nentries = fChain->GetEntriesFast();
	if(event_nb<0 || event_nb>nentries){
		cout << "invalid event number" << endl;
		return;
	}
	TFile * signal_file = new TFile(signal_file_name.c_str(),"READ");
	TTree * signal_tree = (TTree*)(signal_file->Get("T"));

	if(nentries != signal_tree->GetEntriesFast()){
		cout << "total number of event in signal and analyse tree does not match" << endl;
		return;
	}

	float StripAmpl_MG_corr[MG_N][MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()];
	float StripAmpl_CM_corr[CM_N][CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()];
	int signal_evn;
	signal_tree->SetBranchAddress("Nevent",&signal_evn);
	if(CM_N>0) signal_tree->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
	if(MG_N>0) signal_tree->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);

	double chisquare_threshold = 10000;
	int detector_color = 1;
	int ray_color = 2;
	int pos_color = 4;

	LoadTree(event_nb);
	GetEntry(event_nb);
	signal_tree->LoadTree(event_nb);
	signal_tree->GetEntry(event_nb);
	CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,false,-1);
	vector<Ray> eventRays = CBEvent->get_absorption_rays(chisquare_threshold);
	vector<Ray>::iterator rays_it = eventRays.begin();
	//delete rays with too big chi²
	
	while(rays_it!=eventRays.end()){
		if((rays_it->get_chiSquare_X()+rays_it->get_chiSquare_Y())>chisquare_threshold){
			eventRays.erase(rays_it);
			rays_it = eventRays.begin();
		}
		else{
			++rays_it;
		}
	}
	
	TGraph2D * ampl_graph_x = new TGraph2D();
	ampl_graph_x->SetTitle("XZ plane");
	int ampl_graph_x_point_nb = 0;
	TGraph2D * ampl_graph_y = new TGraph2D();
	ampl_graph_y->SetTitle("YZ plane");
	int ampl_graph_y_point_nb = 0;
	vector<pair<TPolyLine3D*, TPolyLine3D*> > det_pos;
	vector<pair<TPolyLine3D*, TPolyLine3D*> > ray_lines;
	vector<TPolyLine3D*> pos_lines_x;
	vector<TPolyLine3D*> pos_lines_y;
	double min_z = numeric_limits<double>::max();
	double max_z = numeric_limits<double>::min();
	//construct line representing det
	for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
		TPolyLine3D det_x(2);
		TPolyLine3D det_y(2);
		double det_z = (*det_it)->get_z();
		if(det_z>max_z) max_z = det_z;
		if(det_z<min_z) min_z = det_z;
		det_x.SetNextPoint(0,det_z,0);
		det_y.SetNextPoint(0,det_z,0);
		det_x.SetNextPoint(500,det_z,0);
		det_y.SetNextPoint(500,det_z,0);
		det_x.SetLineColor(detector_color);
		det_y.SetLineColor(detector_color);
		det_pos.push_back(pair<TPolyLine3D*, TPolyLine3D*>(new TPolyLine3D(det_x),new TPolyLine3D(det_y)));
	}
	min_z = min_z - 0.1*(max_z -min_z);
	max_z = max_z + 0.1*(max_z -min_z);
	for(rays_it = eventRays.begin();rays_it!=eventRays.end();++rays_it){
		//construct line representing ray
		TPolyLine3D ray_x(2);
		TPolyLine3D ray_y(2);
		ray_x.SetNextPoint(rays_it->eval_X(min_z),min_z,0);
		ray_y.SetNextPoint(rays_it->eval_Y(min_z),min_z,0);
		ray_x.SetNextPoint(rays_it->eval_X(max_z),max_z,0);
		ray_y.SetNextPoint(rays_it->eval_Y(max_z),max_z,0);
		ray_x.SetLineColor(ray_color);
		ray_y.SetLineColor(ray_color);
		ray_lines.push_back(pair<TPolyLine3D*, TPolyLine3D*>(new TPolyLine3D(ray_x),new TPolyLine3D(ray_y)));
		vector<Cluster*> ray_clusters = rays_it->get_clus();
	}
	for(vector<Event*>::iterator ev_it = (CBEvent->events).begin();ev_it!=(CBEvent->events).end();++ev_it){
		if((*ev_it)->get_type() == Tomography::MG){
			int det_n_in_tree = (*ev_it)->get_n_in_tree();
			vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*ev_it))->get_clusters();
			for(vector<MG_Cluster>::iterator clus_it = current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
				double det_z = clus_it->get_z();
				double current_clus_size = clus_it->get_size();
				double current_clus_pos = clus_it->get_pos_mm();
				double min_ampl = 0;
				double max_ampl = 0;
				for(int i_strip=0;i_strip<1024;i_strip++){
					int channel = MG_Detector::StripToChannel(i_strip);
					double current_strip_pos = clus_it->correct_strip_nb(i_strip);
					if(current_strip_pos>current_clus_pos+(1.*current_clus_size)) continue;
					if(current_strip_pos<current_clus_pos-(1.*current_clus_size)) continue;
					double current_ampl = *max_element(StripAmpl_MG_corr[det_n_in_tree][channel],StripAmpl_MG_corr[det_n_in_tree][channel]+Tomography::get_instance()->get_Nsample());
					if(current_ampl>max_ampl) max_ampl = current_ampl;
					if(current_ampl<min_ampl) min_ampl = current_ampl;
					if(clus_it->get_is_X()){
						ampl_graph_x->SetPoint(ampl_graph_x_point_nb,current_strip_pos,det_z,current_ampl);
						ampl_graph_x_point_nb++;
					}
					else{
						ampl_graph_y->SetPoint(ampl_graph_y_point_nb,current_strip_pos,det_z,current_ampl);
						ampl_graph_y_point_nb++;
					}
				}
				min_ampl = min_ampl - 0.1*(max_ampl - min_ampl);
				max_ampl = max_ampl + 0.1*(max_ampl - min_ampl);
				TPolyLine3D clus_line(2);
				clus_line.SetNextPoint(current_clus_pos,det_z,min_ampl);
				clus_line.SetNextPoint(current_clus_pos,det_z,max_ampl);
				clus_line.SetLineColor(pos_color);
				clus_line.SetLineStyle(2);
				if(clus_it->get_is_X()) pos_lines_x.push_back(new TPolyLine3D(clus_line));
				else pos_lines_y.push_back(new TPolyLine3D(clus_line));
			}

		}
		else if((*ev_it)->get_type() == Tomography::CM_Demux){
			//TODO : implement for CM
		}		
	}
	delete CBEvent;
	delete signal_file;
	if(ampl_graph_x_point_nb == 0 || ampl_graph_y_point_nb == 0){
		cout << "no cluster in this event, exiting" << endl;
		return;
	}
	ampl_graph_x->GetYaxis()->SetLimits(min_z,max_z);
	ampl_graph_y->GetYaxis()->SetLimits(min_z,max_z);
	ampl_graph_x->GetXaxis()->SetLimits(-50,550);
	ampl_graph_y->GetXaxis()->SetLimits(-50,550);
	ampl_graph_x->SetMarkerStyle(20);
	ampl_graph_y->SetMarkerStyle(20);
	TCanvas * c_visu = new TCanvas();
	c_visu->Divide(2);
	c_visu->cd(1);
	ampl_graph_x->Draw("pcol");
	for(vector<TPolyLine3D*>::iterator line_it = pos_lines_x.begin();line_it!=pos_lines_x.end();++line_it){
		(*line_it)->Draw();
	}
	c_visu->cd(2);
	ampl_graph_y->Draw("pcol");
	for(vector<TPolyLine3D*>::iterator line_it = pos_lines_y.begin();line_it!=pos_lines_y.end();++line_it){
		(*line_it)->Draw();
	}
	for(vector<pair<TPolyLine3D*, TPolyLine3D*> >::iterator pair_it = det_pos.begin();pair_it!=det_pos.end();++pair_it){
		c_visu->cd(1);
		(pair_it->first)->Draw();
		c_visu->cd(2);
		(pair_it->second)->Draw();
	}
	for(vector<pair<TPolyLine3D*, TPolyLine3D*> >::iterator pair_it = ray_lines.begin();pair_it!=ray_lines.end();++pair_it){
		c_visu->cd(1);
		(pair_it->first)->Draw();
		c_visu->cd(2);
		(pair_it->second)->Draw();
	}
	c_visu->Modified();
	c_visu->Update();
	*/
}

void Analyse::Correlation(){
	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
	TCanvas * cDisplay = new TCanvas();
	cDisplay->Divide(4,3);
	TH2D * correlation_X_ampl = new TH2D("correlation_X_ampl","correlation_X_ampl",1000,0,5000,1000,0,5000);
	TH2D * correlation_X_t = new TH2D("correlation_X_t","correlation_X_t",100,0,Tomography::get_instance()->get_Nsample(),100,0,Tomography::get_instance()->get_Nsample());
	TH2D * correlation_Y_ampl = new TH2D("correlation_Y_ampl","correlation_Y_ampl",1000,0,50000,1000,0,50000);
	TH2D * correlation_Y_t = new TH2D("correlation_Y_t","correlation_Y_t",100,0,Tomography::get_instance()->get_Nsample(),100,0,Tomography::get_instance()->get_Nsample());
	TH1D * sigma_X_ampl = new TH1D("sigma_X_ampl","sigma_X_ampl",100,0,5000);
	TH1D * sigma_X_t = new TH1D("sigma_X_t","sigma_X_t",100,0,Tomography::get_instance()->get_Nsample());
	TH1D * sigma_Y_ampl = new TH1D("sigma_Y_ampl","sigma_Y_ampl",100,0,50000);
	TH1D * sigma_Y_t = new TH1D("sigma_Y_t","sigma_Y_t",100,0,Tomography::get_instance()->get_Nsample());
	TH1D * corr_X_ampl = new TH1D("corr_X_ampl","corr_X_ampl",100,-1,1);
	TH1D * corr_X_t = new TH1D("corr_X_t","corr_X_t",100,-1,1);
	TH1D * corr_Y_ampl = new TH1D("corr_Y_ampl","corr_Y_ampl",100,-1,1);
	TH1D * corr_Y_t = new TH1D("corr_Y_t","corr_Y_t",100,-1,1);
	TCanvas * cCorrXY = new TCanvas();
	cCorrXY->Divide(2,2);
	TH2D * correlation_XY_ampl = new TH2D("correlation_XY_ampl","correlation_XY_ampl",1000,0,50000,1000,0,5000);
	TH2D * correlation_XY_t = new TH2D("correlation_XY_t","correlation_XY_t",100,0,Tomography::get_instance()->get_Nsample(),100,0,Tomography::get_instance()->get_Nsample());
	TH1D * corr_XY_ampl = new TH1D("corr_XY_ampl","corr_XY_ampl",100,-1,1);
	TH1D * corr_XY_t = new TH1D("corr_XY_t","corr_XY_t",100,-1,1);
	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	long eventSuitable = 0;
	cout << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	double global_meanXY_ampl = 0;
	double global_meanXY_t = 0;
	double global_sigmaXY_ampl = 0;
	double global_sigmaXY_t = 0;
	double global_meanXY_ampl_perp = 0;
	double global_meanXY_t_perp = 0;
	double global_sigmaXY_ampl_perp = 0;
	double global_sigmaXY_t_perp = 0;
	double global_corrXY_ampl = 0;
	double global_corrXY_t = 0;
	int global_sizeXY = 0;
	for (Long64_t jentry=0; jentry<nentries && Tomography::get_instance()->get_can_continue();jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,-1);
		CBEvent->do_cuts();
		eventSuitable+=CBEvent->get_clus_N()/(get_det_N_tot());
		/*
		bool is_single_event = true;
		for(vector<Detector*>::iterator it=detectors.begin();it!=detectors.end();++it){
			if((CBEvent->get_clus_N_by_det(*it))!=1) is_single_event = false;
		}
		if(is_single_event){
			for(vector<Event*>::iterator it=(CBEvent->events).begin();it!=(CBEvent->events).end();++it){
				double ampl_1 = 0;
				double t_1 = 0;
				if((*it)->get_type() == Tomography::MG){
					MG_Cluster current_cluster = ((dynamic_cast<MG_Event*>(*it))->get_clusters()).front();
					ampl_1 = current_cluster.get_ampl();
					t_1 = current_cluster.get_t();
				}
				else if((*it)->get_type() == Tomography::CM_Demux){
					CM_Demux_Cluster current_cluster = ((dynamic_cast<CM_Demux_Event*>(*it))->get_clusters()).front();
					ampl_1 = current_cluster.get_ampl();
					t_1 = current_cluster.get_t();
				}
				bool is_X = (*it)->get_is_X();
				for(vector<Event*>::iterator jt=(it+1);jt!=(CBEvent->events).end();++jt){
					if(is_X != (*jt)->get_is_X()) continue;
					double ampl_2 = 0;
					double t_2 = 0;
					if((*it)->get_type() == Tomography::MG){
						MG_Cluster current_cluster = ((dynamic_cast<MG_Event*>(*jt))->get_clusters()).front();
						ampl_2 = current_cluster.get_ampl();
						t_2 = current_cluster.get_t();
					}
					else if((*it)->get_type() == Tomography::CM_Demux){
						CM_Demux_Cluster current_cluster = ((dynamic_cast<CM_Demux_Event*>(*jt))->get_clusters()).front();
						ampl_2 = current_cluster.get_ampl();
						t_2 = current_cluster.get_t();
					}
					if(is_X){
						correlation_X_ampl->Fill(ampl_1,ampl_2);
						correlation_X_t->Fill(t_1,t_2);
					}
					else{
						correlation_Y_ampl->Fill(ampl_1,ampl_2);
						correlation_Y_t->Fill(t_1,t_2);
					}
				}
			}
		}
		*/

		vector<Ray> currentRays = CBEvent->get_absorption_rays();
		for(vector<Ray>::iterator ray_it = currentRays.begin();ray_it!=currentRays.end();++ray_it){
			vector<Cluster*> currentClusters = ray_it->get_clus();
			double sigmaX_ampl = 0;
			double sigmaY_ampl = 0;
			double meanX_ampl = 0;
			double meanY_ampl = 0;
			double sigmaX_t = 0;
			double sigmaY_t = 0;
			double meanX_t = 0;
			double meanY_t = 0;
			double corrX_ampl = 0;
			double corrY_ampl = 0;
			double corrX_t = 0;
			double corrY_t = 0;
			int sizeX_corr = 0;
			int sizeY_corr = 0;
			int sizeX = 0;
			int sizeY = 0;
			for(vector<Cluster*>::iterator clus_it = currentClusters.begin();clus_it!=currentClusters.end();++clus_it){
				double ampl_1 = (*clus_it)->get_ampl();
				double t_1 = (*clus_it)->get_t();
				//double t_1 = (*clus_it)->get_maxSample();
				bool is_X = (*clus_it)->get_is_X();
				if(is_X){
					meanX_ampl += ampl_1;
					sigmaX_ampl += ampl_1*ampl_1;
					meanX_t += t_1; 
					sigmaX_t += t_1*t_1;
					sizeX++;
				}
				else{
					meanY_ampl += ampl_1;
					sigmaY_ampl += ampl_1*ampl_1;
					meanY_t += t_1; 
					sigmaY_t += t_1*t_1;
					sizeY++;
				}
				for(vector<Cluster*>::iterator clus_jt = (clus_it+1);clus_jt!=currentClusters.end();++clus_jt){
					if(is_X != (*clus_jt)->get_is_X()) continue;
					if(is_X){
						correlation_X_ampl->Fill(ampl_1,(*clus_jt)->get_ampl());
						correlation_X_t->Fill(t_1,(*clus_jt)->get_t());
						//correlation_X_t->Fill(t_1,(*clus_jt)->get_maxSample());
						corrX_ampl += ampl_1*((*clus_jt)->get_ampl());
						corrX_t += t_1*((*clus_jt)->get_t());
						//corrX_t += t_1*((*clus_jt)->get_maxSample());
						sizeX_corr++;
					}
					else{
						correlation_Y_ampl->Fill(ampl_1,(*clus_jt)->get_ampl());
						correlation_Y_t->Fill(t_1,(*clus_jt)->get_t());
						//correlation_Y_t->Fill(t_1,(*clus_jt)->get_maxSample());
						corrY_ampl += ampl_1*((*clus_jt)->get_ampl());
						corrY_t += t_1*((*clus_jt)->get_t());
						//corrY_t += t_1*((*clus_jt)->get_maxSample());
						sizeY_corr++;
					}
				}
			}
			meanX_ampl /= sizeX;
			meanX_t /= sizeX;
			meanY_ampl /= sizeY;
			meanY_t /= sizeY;
			sigmaX_ampl = Sqrt((sigmaX_ampl/sizeX) - (meanX_ampl*meanX_ampl));
			sigmaX_t = Sqrt((sigmaX_t/sizeX) - (meanX_t*meanX_t));
			sigmaY_ampl = Sqrt((sigmaY_ampl/sizeY) - (meanY_ampl*meanY_ampl));
			sigmaY_t = Sqrt((sigmaY_t/sizeY) - (meanY_t*meanY_t));
			sigma_X_ampl->Fill(sigmaX_ampl);
			sigma_X_t->Fill(sigmaX_t);
			sigma_Y_ampl->Fill(sigmaY_ampl);
			sigma_Y_t->Fill(sigmaY_t);
			corrX_ampl = (corrX_ampl/sizeX_corr - meanX_ampl*meanX_ampl)/(sigmaX_ampl*sigmaX_ampl);
			corrX_t = (corrX_t/sizeX_corr - meanX_t*meanX_t)/(sigmaX_t*sigmaX_t);
			corrY_ampl = (corrY_ampl/sizeY_corr - meanY_ampl*meanY_ampl)/(sigmaY_ampl*sigmaY_ampl);
			corrY_t = (corrY_t/sizeY_corr - meanY_t*meanY_t)/(sigmaY_t*sigmaY_t);
			corr_X_ampl->Fill(corrX_ampl);
			corr_X_t->Fill(corrX_t);
			corr_Y_ampl->Fill(corrY_ampl);
			corr_Y_t->Fill(corrY_t);
			double meanXY_ampl = 0;
			double meanXY_t = 0;
			double sigmaXY_ampl = 0;
			double sigmaXY_t = 0;
			double meanXY_ampl_perp = 0;
			double meanXY_t_perp = 0;
			double sigmaXY_ampl_perp = 0;
			double sigmaXY_t_perp = 0;
			double corrXY_ampl = 0;
			double corrXY_t = 0;
			int sizeXY = 0;
			for(vector<Detector*>::iterator det_it = detectors.begin();det_it!=detectors.end();++det_it){
				if((*det_it)->get_perp_n()<0) continue;
				if((*det_it)->get_is_X()) continue;
				int n = (*det_it)->get_n_in_tree();
				vector<Detector*>::iterator det_jt = detectors.begin();
				while((*det_jt)->get_perp_n() != n && det_jt != detectors.end()){
					++det_jt;
				}
				double amplXY = 0;
				double amplXY_perp = 0;
				double tXY = 0;
				double tXY_perp = 0;
				int check = 0;
				for(vector<Cluster*>::iterator clus_it = currentClusters.begin();clus_it!=currentClusters.end();++clus_it){
					if((*clus_it)->is_in_det(*det_it)){
						amplXY = (*clus_it)->get_ampl();
						tXY = (*clus_it)->get_t();
						//tXY = (*clus_it)->get_maxSample();
						check++;
					}
					else if((*clus_it)->is_in_det(*det_jt)){
						amplXY_perp = (*clus_it)->get_ampl();
						tXY_perp = (*clus_it)->get_t();
						//tXY_perp = (*clus_it)->get_maxSample();
						check++;
					}
				}
				if(check!=2) continue;
				meanXY_ampl += amplXY;
				meanXY_t += tXY;
				sigmaXY_ampl += amplXY*amplXY;
				sigmaXY_t += tXY*tXY;
				meanXY_ampl_perp += amplXY_perp;
				meanXY_t_perp += tXY_perp;
				sigmaXY_ampl_perp += amplXY_perp*amplXY_perp;
				sigmaXY_t_perp += tXY_perp*tXY_perp;
				corrXY_ampl += amplXY*amplXY_perp;
				corrXY_t += tXY*tXY_perp;
				sizeXY++;
				correlation_XY_ampl->Fill(amplXY,amplXY_perp);
				correlation_XY_t->Fill(tXY,tXY_perp);
			}
			if(sizeXY>0){
				global_meanXY_ampl += meanXY_ampl;
				global_meanXY_t += meanXY_t;
				global_sigmaXY_ampl += sigmaXY_ampl;
				global_sigmaXY_t += sigmaXY_t;
				global_meanXY_ampl_perp += meanXY_ampl_perp;
				global_meanXY_t_perp += meanXY_t_perp;
				global_sigmaXY_ampl_perp += sigmaXY_ampl_perp;
				global_sigmaXY_t_perp += sigmaXY_t_perp;
				global_corrXY_ampl += corrXY_ampl;
				global_corrXY_t += corrXY_t;
				global_sizeXY += sizeXY;

				meanXY_ampl /= sizeXY;
				meanXY_t /= sizeXY;
				sigmaXY_ampl = Sqrt((sigmaXY_ampl/sizeXY) - (meanXY_ampl*meanXY_ampl));
				sigmaXY_t = Sqrt((sigmaXY_t/sizeXY) - (meanXY_t*meanXY_t));
				meanXY_ampl_perp /= sizeXY;
				meanXY_t_perp /= sizeXY;
				sigmaXY_ampl_perp = Sqrt((sigmaXY_ampl_perp/sizeXY) - (meanXY_ampl_perp*meanXY_ampl_perp));
				sigmaXY_t_perp = Sqrt((sigmaXY_t_perp/sizeXY) - (meanXY_t_perp*meanXY_t_perp));
				corrXY_ampl = ((corrXY_ampl/sizeXY) - (meanXY_ampl*meanXY_ampl_perp))/(sigmaXY_ampl*sigmaXY_ampl_perp);
				corrXY_t = ((corrXY_t/sizeXY) - (meanXY_t*meanXY_t_perp))/(sigmaXY_t*sigmaXY_t_perp);
				corr_XY_ampl->Fill(corrXY_ampl);
				corr_XY_t->Fill(corrXY_t);
			}
			for(vector<Cluster*>::iterator clus_it = currentClusters.begin();clus_it!=currentClusters.end();++clus_it){
				delete *clus_it;
			}
		}

		if(jentry%500 == 0) cout << "\r" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			cDisplay->cd(1);
			correlation_X_ampl->Draw("colz");
			cDisplay->cd(2);
			correlation_Y_ampl->Draw("colz");
			cDisplay->cd(3);
			correlation_X_t->Draw("colz");
			cDisplay->cd(4);
			correlation_Y_t->Draw("colz");
			cDisplay->cd(5);
			sigma_X_ampl->Draw();
			cDisplay->cd(6);
			sigma_Y_ampl->Draw();
			cDisplay->cd(7);
			sigma_X_t->Draw();
			cDisplay->cd(8);
			sigma_Y_t->Draw();
			cDisplay->cd(9);
			corr_X_ampl->Draw();
			cDisplay->cd(10);
			corr_Y_ampl->Draw();
			cDisplay->cd(11);
			corr_X_t->Draw();
			cDisplay->cd(12);
			corr_Y_t->Draw();
			cDisplay->Modified();
			cDisplay->Update();
			cCorrXY->cd(1);
			correlation_XY_ampl->Draw("colz");
			cCorrXY->cd(2);
			correlation_XY_t->Draw("colz");
			cCorrXY->cd(3);
			corr_XY_ampl->Draw();
			cCorrXY->cd(4);
			corr_XY_t->Draw();
			cCorrXY->Modified();
			cCorrXY->Update();
		}
		delete CBEvent;
	}
	cout << "\r" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	cDisplay->cd(1);
	correlation_X_ampl->Draw("colz");
	cDisplay->cd(2);
	correlation_Y_ampl->Draw("colz");
	cDisplay->cd(3);
	correlation_X_t->Draw("colz");
	cDisplay->cd(4);
	correlation_Y_t->Draw("colz");
	cDisplay->cd(5);
	sigma_X_ampl->Draw();
	cDisplay->cd(6);
	sigma_Y_ampl->Draw();
	cDisplay->cd(7);
	sigma_X_t->Draw();
	cDisplay->cd(8);
	sigma_Y_t->Draw();
	cDisplay->cd(9);
	corr_X_ampl->Draw();
	cDisplay->cd(10);
	corr_Y_ampl->Draw();
	cDisplay->cd(11);
	corr_X_t->Draw();
	cDisplay->cd(12);
	corr_Y_t->Draw();
	cDisplay->Modified();
	cDisplay->Update();
	cCorrXY->cd(1);
	correlation_XY_ampl->Draw("colz");
	cCorrXY->cd(2);
	correlation_XY_t->Draw("colz");
	cCorrXY->cd(3);
	corr_XY_ampl->Draw();
	cCorrXY->cd(4);
	corr_XY_t->Draw();
	cCorrXY->Modified();
	cCorrXY->Update();
	global_meanXY_ampl /= global_sizeXY;
	global_meanXY_t /= global_sizeXY;
	global_sigmaXY_ampl = Sqrt((global_sigmaXY_ampl/global_sizeXY) - (global_meanXY_ampl*global_meanXY_ampl));
	global_sigmaXY_t = Sqrt((global_sigmaXY_t/global_sizeXY) - (global_meanXY_t*global_meanXY_t));
	global_meanXY_ampl_perp /= global_sizeXY;
	global_meanXY_t_perp /= global_sizeXY;
	global_sigmaXY_ampl_perp = Sqrt((global_sigmaXY_ampl_perp/global_sizeXY) - (global_meanXY_ampl_perp*global_meanXY_ampl_perp));
	global_sigmaXY_t_perp = Sqrt((global_sigmaXY_t_perp/global_sizeXY) - (global_meanXY_t_perp*global_meanXY_t_perp));
	global_corrXY_ampl = ((global_corrXY_ampl/global_sizeXY) - (global_meanXY_ampl*global_meanXY_ampl_perp))/(global_sigmaXY_ampl*global_sigmaXY_ampl_perp);
	global_corrXY_t = ((global_corrXY_t/global_sizeXY) - (global_meanXY_t*global_meanXY_t_perp))/(global_sigmaXY_t*global_sigmaXY_t_perp);
	cout << "correlation XY en amplitude : " <<  global_corrXY_ampl << endl;
	cout << "correlation XY en temps : " << global_corrXY_t << endl;
}

void Analyse::SignalOverNoise(){
	map<string,TCanvas*> cDisplay;
	map<string,TH1D*> global_signal;
	map<string,TH1D*> global_noise;
	map<string,TProfile*> global_signal_over_noise;
	for(vector<Detector*>::iterator it=detectors.begin();it!=detectors.end();++it){
		ostringstream name;
		name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
		cDisplay[name.str()] = new TCanvas(name.str().c_str());
		cDisplay[name.str()]->Divide(3);
		//name << "_";
		global_signal[name.str()] = new TH1D((name.str() + "_signal").c_str(),(name.str() + "_signal").c_str(),500,0,4000);
		global_noise[name.str()] = new TH1D((name.str() + "_noise").c_str(),(name.str() + "_noise").c_str(),100,0,100);
		for(int i=0;i<(*it)->get_Nchannel();i++){
			global_noise[name.str()]->Fill((*it)->get_RMS(i));
		}
		global_signal_over_noise[name.str()] = new TProfile((name.str() + "_SoB").c_str(),(name.str() + "_SoB").c_str(),(*it)->get_Nchannel(),0,(*it)->get_Nchannel());
	}
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries && Tomography::get_instance()->get_can_continue();i++){
		LoadTree(i);
		GetEntry(i);

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			ostringstream name;
			name << (*it)->get_type() << "_" << (*it)->get_n_in_tree();
			Event * current_event = (*it)->build_event(this,evn);
			vector<Cluster*> current_cluster = current_event->get_clusters();
			for(vector<Cluster*>::iterator jt=current_cluster.begin();jt!=current_cluster.end();++jt){
				double current_signal = (*jt)->get_maxStripAmpl();
				global_signal[name.str()]->Fill(current_signal);
				int channel = (*it)->StripToChannel((*jt)->get_maxStrip());
				double current_noise = (*it)->get_RMS(channel);
				global_signal_over_noise[name.str()]->Fill(channel,current_signal/current_noise);
				delete *jt;
			}
			delete current_event;
		}
		if(i%100 == 0) cout << "\r" << i << "/" << nentries << flush;
		if(i%5000 == 0 && Tomography::get_instance()->get_live_graphic_display()){
			for(map<string,TCanvas*>::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
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
	
	for(map<string,TCanvas*>::iterator it = cDisplay.begin();it!=cDisplay.end();++it){
		it->second->cd(1);
		TFitResultPtr res_signal = global_signal[it->first]->Fit("landau","SQ");
		global_signal[it->first]->Draw();
		it->second->cd(2);
		TFitResultPtr res_noise = global_noise[it->first]->Fit("gaus","SQ");
		global_noise[it->first]->Draw();
		it->second->cd(3);
		TLine * average_SoN = new TLine(0,(res_signal->Parameter(1))/(res_noise->Parameter(1)),MG_Detector::Nchannel,(res_signal->Parameter(1))/(res_noise->Parameter(1)));
		average_SoN->SetLineStyle(2);
		average_SoN->SetLineColor(2);
		TLine * mean_SoN = new TLine(0,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()),MG_Detector::Nchannel,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()));
		mean_SoN->SetLineStyle(2);
		mean_SoN->SetLineColor(4);
		global_signal_over_noise[it->first]->Draw();
		average_SoN->Draw();
		mean_SoN->Draw();
		it->second->Modified();
		it->second->Update();
		cout << it->first << endl;
		cout << "    mean S/B : " << global_signal_over_noise[it->first]->GetMean(2) << endl;
		cout << "    mean S/mean B : " << (global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()) << endl;
		cout << "    sigma S/B : " << global_signal_over_noise[it->first]->GetRMS(2) << endl;
		cout << "    delta mean S/B : " << global_signal_over_noise[it->first]->GetMean(11) << endl;
	}
}
