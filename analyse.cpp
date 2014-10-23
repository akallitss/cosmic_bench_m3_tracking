#define analyse_cpp
//my class
#include "analyse.h"
#include "T.h"
#include "detector.h"
#include "ray.h"
#include "event.h"
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
//Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
//std
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include <utility>

//std
using std::string;
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
//boost
using boost::property_tree::ptree;
//tmath
using TMath::Sqrt;
using TMath::ATan;
using TMath::MaxElement;
using TMath::Abs;
using TMath::FloorNint;
using TMath::CeilNint;

Analyse::Analyse(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	TFile *f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	CM_N = 0;
	MG_N = 0;
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.CosMultis")){
		detectors.push_back(new CM_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("is_up"),child.second.get<int>("cm_n"),child.second.get<bool>("use_thin_strip"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle")));
		detectors.back()->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		(dynamic_cast<CM_Detector*>(detectors.back()))->set_ClusMaxStripAmplCut_Min_Wide(child.second.get<double>("ClusMaxStripAmplCut_Min_Wide"));
		(dynamic_cast<CM_Detector*>(detectors.back()))->set_ClusSizeCut_Max_Wide(child.second.get<double>("ClusSizeCut_Max_Wide"));
		CM_N++;
	}
	BOOST_FOREACH(const ptree::value_type& child, config_tree.get_child("CosmicBench.MultiGens")){
		detectors.push_back(new MG_Detector(child.second.get<double>("z"),child.second.get<bool>("is_X"),child.second.get<bool>("is_up"),child.second.get<int>("mg_n"),child.second.get<bool>("is_ref"),child.second.get<double>("offset"),child.second.get<bool>("direction"),child.second.get<double>("angle")));
		detectors.back()->set_ClusTOTCut_Min(child.second.get<double>("ClusTOTCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Min(child.second.get<double>("ClusMaxSampleCut_Min"));
		detectors.back()->set_ClusMaxSampleCut_Max(child.second.get<double>("ClusMaxSampleCut_Max"));
		(dynamic_cast<MG_Detector*>(detectors.back()))->set_ClusSizeCut_Min(child.second.get<double>("ClusSizeCut_Min"));
		MG_N++;
	}
	int total_CM_N = config_tree.get<int>("total_CM_N");
	int total_MG_N = config_tree.get<int>("total_MG_N");
	if((total_CM_N!=CM_N) || (total_MG_N!=MG_N)){
		cout << "problem in detectors number" << endl;
		return;
	}
	Init(tree,total_CM_N,total_MG_N,config_tree.get<bool>("use_SRF"));
	signal_file_name = config_tree.get<string>("signal_file");
}
Analyse::~Analyse(){

}
void Analyse::Residus(){
	TCanvas * c_CM[CM_N];
	TCanvas * c_MG[MG_N];
	TH1D * CM_residus[CM_N];
	TH1D * MG_residus[MG_N];
	int nbins = 200;
	int lim = 500;
	int eventReconstructed = 0;
	int eventSuitable = 0;
	for(int i=0;i<CM_N;i++){
		ostringstream name;
		name << "Cosmulti_" << i;
		c_CM[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		name << "_residu";
		CM_residus[i] = new TH1D(name.str().c_str(),name.str().c_str(),nbins,-lim,lim);
	}
	for(int i=0;i<MG_N;i++){
		ostringstream name;
		name << "MultiGen_" << i;
		c_MG[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		name << "_residu";
		MG_residus[i] = new TH1D(name.str().c_str(),name.str().c_str(),nbins,-lim,lim);
	}
	TCanvas * c0 = new TCanvas("chiSquare","chiSquare");
	TH1D * chiSquareH = new TH1D("chiSquares","chiSquares",nbins,0,3*lim);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
				chiSquareH->Fill(it->get_chiSquare_X()+it->get_chiSquare_Y());
			}
			int i_MG = 0;
			int i_CM = 0;
			for(vector<Detector*>::iterator jt=detectors.begin();jt!=detectors.end();++jt){
				double residu = it->get_residu(*jt);
				if(residu != numeric_limits<double>::min()){
					if((*jt)->get_type() == "CM"){
						CM_Detector * currentDet = dynamic_cast<CM_Detector*>(*jt);
						CM_residus[i_CM]->Fill(it->get_residu(currentDet));
						i_CM++;
					}
					if((*jt)->get_type() == "MG"){
						MG_Detector * currentDet = dynamic_cast<MG_Detector*>(*jt);
						MG_residus[i_MG]->Fill(it->get_residu(currentDet));
						i_MG++;
					}
				}
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
			for(int i=0;i<CM_N;i++){
				c_CM[i]->cd();
				CM_residus[i]->Draw();
				c_CM[i]->Modified();
				c_CM[i]->Update();
			}
			for(int i=0;i<MG_N;i++){
				c_MG[i]->cd();
				MG_residus[i]->Draw();
				c_MG[i]->Modified();
				c_MG[i]->Update();
			}
			c0->cd();
			chiSquareH->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	for(int i=0;i<CM_N;i++){
		c_CM[i]->cd();
		CM_residus[i]->Draw();
		c_CM[i]->Modified();
		c_CM[i]->Update();
	}
	for(int i=0;i<MG_N;i++){
		c_MG[i]->cd();
		MG_residus[i]->Draw();
		c_MG[i]->Modified();
		c_MG[i]->Update();
	}
	c0->cd();
	chiSquareH->Draw();
	c0->Modified();
	c0->Update();
}
void Analyse::Residus_ref(){
	int non_ref_n = 0;
	double chisquare_threshold = 10;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	gStyle->SetPalette(55,0);
	map<string,TCanvas*> c_MM;
	map<string,TH1D*> MM_residus;
	map<string,TH2D*> muon_seen;
	map<string,TH2D*> muon_total;
	map<string,TH2D*> efficacity_2D;
	map<string,TGraph*> correlation;
	map<string,TProfile*> angle_alignment;
	map<string,int> point_nb;
	map<string,double> efficacity;
	map<string,bool> is_seen;
	int nbins = 200;
	int lim = 500;
	double marge = 2./5.;
	int nbins_2D = 50*(1+2*marge);
	int eventReconstructed = 0;
	double eventSuitable = 0;
	map<string, unsigned int> det_in_nref_dir;
	map<string, bool> nref_is_X;
	unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			if((*it)->get_type() == "CM"){
				name << "Cosmulti_" << (dynamic_cast<CM_Detector*>(*it))->get_cm_n_in_tree();
				c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
				c_MM[name.str()]->Divide(2,2);
				MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
				muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				efficacity_2D[name.str()]->SetStats(false);
				correlation[name.str()] = new TGraph();
				angle_alignment[name.str()] = new TProfile((name.str()+"_angle").c_str(),(name.str()+"_angle").c_str(),500,0,500,-5,5);
				point_nb[name.str()] = 0;
				efficacity[name.str()] = 0;
				is_seen [name.str()] = false;
			}
			else if((*it)->get_type() == "MG"){
				name << "Multigen_" << (dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree();
				c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
				c_MM[name.str()]->Divide(2,2);
				MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
				muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
				efficacity_2D[name.str()]->SetStats(false);
				correlation[name.str()] = new TGraph();
				angle_alignment[name.str()] = new TProfile((name.str()+"_angle").c_str(),(name.str()+"_angle").c_str(),500,0,500,-5,5);
				point_nb[name.str()] = 0;
				efficacity[name.str()] = 0;
				is_seen [name.str()] = false;
			}
			nref_is_X[name.str()] = (*it)->get_is_X();
			if((*it)->get_is_X()) nref_x_n++;
		}
		if((*it)->get_is_X()) det_x_n++;
	}
	for(map<string,bool>::iterator it = nref_is_X.begin();it!=nref_is_X.end();++it){
		if(it->second) det_in_nref_dir[it->first] = det_x_n - nref_x_n;
		else det_in_nref_dir[it->first] = (MG_N + CM_N - det_x_n) - (nref_is_X.size() - nref_x_n);
	}
	/*
	c_MM["Multigen_2D_0"] = new TCanvas("Multigen_2D_0","Multigen_2D_0",1200,1000);
	c_MM["Multigen_2D_0"]->Divide(2);
	MM_residus["Multigen_2D_0"] = new TH1D("Multigen_2D_0_residu","Multigen_2D_0_residu",nbins,-lim,lim);
	muon_total["Multigen_2D_0"] = new TH2D("Multigen_2D_0_total","Multigen_2D_0_total",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	muon_seen["Multigen_2D_0"] = new TH2D("Multigen_2D_0_seen","Multigen_2D_0_seen",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	efficacity_2D["Multigen_2D_0"] = new TH2D("Multigen_2D_0_efficacity","Multigen_2D_0_efficacity",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	efficacity["Multigen_2D_0"] = 0;
	point_nb["Multigen_2D_0"] = 0;
	correlation["Multigen_2D_0"] = new TGraph();
	double z_MG2D_0 = 685;

	c_MM["Multigen_2D_1"] = new TCanvas("Multigen_2D_1","Multigen_2D_1",1200,1000);
	c_MM["Multigen_2D_1"]->Divide(2);
	MM_residus["Multigen_2D_1"] = new TH1D("Multigen_2D_1_residu","Multigen_2D_1_residu",nbins,-lim,lim);
	muon_total["Multigen_2D_1"] = new TH2D("Multigen_2D_1_total","Multigen_2D_1_total",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	muon_seen["Multigen_2D_1"] = new TH2D("Multigen_2D_1_seen","Multigen_2D_1_seen",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	efficacity_2D["Multigen_2D_1"] = new TH2D("Multigen_2D_1_efficacity","Multigen_2D_1_efficacity",nbins_2D,-marge*lim,(1+marge)*lim,nbins_2D,-marge*lim,(1+marge)*lim);
	efficacity["Multigen_2D_1"] = 0;
	point_nb["Multigen_2D_1"] = 0;
	correlation["Multigen_2D_1"] = new TGraph();
	double z_MG2D_1 = 705;
	*/
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",MG_N + CM_N + 2,0,MG_N + CM_N + 2);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
			if((jt->get_chiSquare_X()+jt->get_chiSquare_Y()) < chisquare_threshold){
				chisquares->Fill(jt->get_chiSquare_X()+jt->get_chiSquare_Y());
				ray_clus_n->Fill(jt->get_clus_n());
			}
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(CM_N+MG_N);
		
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				if((*it)->get_type() == "CM_Demux"){
					ostringstream name;
					name << "Cosmulti_" << (*it)->get_n_in_tree();
					vector<CM_Demux_Cluster> current_clusters = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						for(map<string,bool>::iterator nt = is_seen.begin();nt!=is_seen.end();++nt){
							nt->second = false;
						}
						double residu = numeric_limits<double>::max();
						vector<CM_Demux_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<CM_Demux_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							if(kt->get_is_X()){
								kt->set_perp_pos_mm(jt->eval_Y((*it)->get_z()));
							}
							else{
								kt->set_perp_pos_mm(jt->eval_X((*it)->get_z()));
							}
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						if(jt->eval_X((*it)->get_z())<500 && jt->eval_X((*it)->get_z())>0 && jt->eval_Y((*it)->get_z())<500 && jt->eval_Y((*it)->get_z())>0) muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						if(matching_cluster == current_clusters.end()) continue;
						if((*matching_cluster).get_is_X()){
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),(*matching_cluster).get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						}
						else{
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),(*matching_cluster).get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						}
						point_nb[name.str()]++;
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
						if(residu<chisquare_threshold){
							muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
							is_seen[name.str()] = true;
						}
					}
											/*
						muon_total["Multigen_2D_0"]->Fill(jt->eval_X(z_MG2D_0),jt->eval_Y(z_MG2D_0));
						muon_total["Multigen_2D_1"]->Fill(jt->eval_X(z_MG2D_1),jt->eval_Y(z_MG2D_1));
						if(is_seen["Multigen_0"] && is_seen["Multigen_1"]){
							muon_seen["Multigen_2D_0"]->Fill(jt->eval_X(z_MG2D_0),jt->eval_Y(z_MG2D_0));
						}
						if(is_seen["Multigen_2"] && is_seen["Multigen_3"]){
							muon_seen["Multigen_2D_1"]->Fill(jt->eval_X(z_MG2D_1),jt->eval_Y(z_MG2D_1));
						}
						*/
				}
				else if((*it)->get_type() == "MG"){
					ostringstream name;
					name << "Multigen_" << (*it)->get_n_in_tree();
					vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						for(map<string,bool>::iterator nt = is_seen.begin();nt!=is_seen.end();++nt){
							nt->second = false;
						}
						double residu = numeric_limits<double>::max();
						vector<MG_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<MG_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							if(kt->get_is_X()){
								kt->set_perp_pos_mm(jt->eval_Y((*it)->get_z()));
							}
							else{
								kt->set_perp_pos_mm(jt->eval_X((*it)->get_z()));
							}
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						if(jt->eval_X((*it)->get_z())<500 && jt->eval_X((*it)->get_z())>0 && jt->eval_Y((*it)->get_z())<500 && jt->eval_Y((*it)->get_z())>0) muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						if(matching_cluster == current_clusters.end()) continue;
						if((*matching_cluster).get_is_X()){
							(*matching_cluster).set_perp_pos_mm(jt->eval_Y((*it)->get_z()));
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),(*matching_cluster).get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
						}
						else{
							(*matching_cluster).set_perp_pos_mm(jt->eval_X((*it)->get_z()));
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),(*matching_cluster).get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
						}
						point_nb[name.str()]++;
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
						if(residu<chisquare_threshold){
							muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
							is_seen[name.str()] = true;
						}
						/*
						muon_total["Multigen_2D_0"]->Fill(jt->eval_X(z_MG2D_0),jt->eval_Y(z_MG2D_0));
						muon_total["Multigen_2D_1"]->Fill(jt->eval_X(z_MG2D_1),jt->eval_Y(z_MG2D_1));
						if(is_seen["Multigen_0"] && is_seen["Multigen_1"]){
							muon_seen["Multigen_2D_0"]->Fill(jt->eval_X(z_MG2D_0),jt->eval_Y(z_MG2D_0));
						}
						if(is_seen["Multigen_2"] && is_seen["Multigen_3"]){
							muon_seen["Multigen_2D_1"]->Fill(jt->eval_X(z_MG2D_1),jt->eval_Y(z_MG2D_1));
						}
						*/
					}
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<int>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
			for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
				it->second->cd(1);
				MM_residus[it->first]->Draw();
				for(int i=1;i<=nbins;i++){
					for(int j=1;j<=nbins;j++){
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
				it->second->Modified();
				it->second->Update();
			}
			c0->cd(1);
			chisquares->Draw();
			c0->cd(2);
			ray_clus_n->Draw();
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<int>(eventSuitable) << "|" << setw(20) << nentries << endl;
	for(map<string,TCanvas*>::iterator it = c_MM.begin();it!=c_MM.end();++it){
		it->second->cd(1);
		MM_residus[it->first]->Draw();
		double total_seen = 0;
		double total_passed = 0;
		for(int i=1;i<=nbins;i++){
			for(int j=1;j<=nbins;j++){
				int binN = muon_total[it->first]->GetBin(i,j);
				double binContent = 0;
				if(muon_total[it->first]->GetBinContent(binN) > 0) binContent = (muon_seen[it->first]->GetBinContent(binN))/(muon_total[it->first]->GetBinContent(binN));
				efficacity_2D[it->first]->SetBinContent(binN,binContent);
				double pos_X = muon_total[it->first]->GetXaxis()->GetBinCenter(i);
				double pos_Y = muon_total[it->first]->GetYaxis()->GetBinCenter(j);
				if(pos_X<=400 && pos_X>=100 && pos_Y<=400 && pos_Y>=100){
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
		angle_alignment[it->first]->Draw();
		it->second->Modified();
		it->second->Update();
		cout << it->first << " efficacity : " << 100.*efficacity[it->first] << "%" << endl;
	}
	//cout << "correlation MG2D_0 : " << (efficacity["Multigen_2D_0"] - efficacity["Multigen_0"]*efficacity["Multigen_1"])/Sqrt(efficacity["Multigen_0"]*(1-efficacity["Multigen_0"])*efficacity["Multigen_1"]*(1-efficacity["Multigen_1"])) << endl;
	//cout << "correlation MG2D_1 : " << (efficacity["Multigen_2D_1"] - efficacity["Multigen_2"]*efficacity["Multigen_3"])/Sqrt(efficacity["Multigen_2"]*(1-efficacity["Multigen_2"])*efficacity["Multigen_3"]*(1-efficacity["Multigen_3"])) << endl;
	c0->cd(1);
	chisquares->Draw();
	c0->cd(2);
	ray_clus_n->Draw();
	c0->Modified();
	c0->Update();
}
void Analyse::Efficacity(){
	TCanvas * c_CM[CM_N];
	TCanvas * c_MG[MG_N];
	TProfile * CM_residus[CM_N];
	TProfile * CM_spark_h[CM_N];
	TProfile * MG_residus[MG_N];
	TProfile * MG_spark_h[MG_N];
	int nbins = 200;
	//int lim = 100;
	int eventReconstructed = 0;
	int eventSuitable = 0;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	LoadTree(0);
	fChain->GetEntry(0);
	double beginTime = evttime;
	LoadTree(nentries-1);
	fChain->GetEntry(nentries-1);
	double endTime = evttime;

	for(int i=0;i<CM_N;i++){
		ostringstream name;
		name << "Cosmulti_" << i;
		c_CM[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		CM_residus[i] = new TProfile((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,beginTime,endTime,0,2);
		CM_residus[i]->SetLineColor(4);
		CM_spark_h[i] = new TProfile((name.str()+"_spark").c_str(),(name.str()+"_spark").c_str(),nbins,beginTime,endTime,0,2);
		CM_spark_h[i]->SetLineColor(2);
	}
	for(int i=0;i<MG_N;i++){
		ostringstream name;
		name << "MultiGen_" << i;
		c_MG[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		MG_residus[i] = new TProfile((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,beginTime,endTime,0,2);
		MG_residus[i]->SetLineColor(4);
		MG_spark_h[i] = new TProfile((name.str()+"_spark").c_str(),(name.str()+"_spark").c_str(),nbins,beginTime,endTime,0,2);
		MG_spark_h[i]->SetLineColor(2);
	}
	TCanvas * c0 = new TCanvas("chiSquare","chiSquare");
	TProfile * chiSquareH = new TProfile("chiSquares","chiSquares",nbins,0,nentries,-10,200);
	TCanvas * c_t = new TCanvas("diff_T","diff_T");
	//TH1D * diff_t = new TH1D("diff_t","diff_t",nbins,-20,20);
	TH1D * tLO = new TH1D("tLO","tLO",nbins,0,40);
	TH1D * tNLO = new TH1D("tNLO","tNLO",nbins,0,40);
	tLO->SetLineColor(2);
	tNLO->SetLineColor(4);
	TCanvas * c_sigma = new TCanvas("sigma_T","sigma_T");
	TH1D * sigmaLO = new TH1D("sigmaLO","sigmaLO",nbins,0,20);
	TH1D * sigmaNLO = new TH1D("sigmaNLO","sigmaNLO",nbins,0,20);
	sigmaLO->SetLineColor(2);
	sigmaNLO->SetLineColor(4);
	//TH1D * sigma_t = new TH1D("sigma_t","sigma_t",nbins,0,40);
	TCanvas * c_test = new TCanvas("test","test");
	TH1D * test_order = new TH1D("test_order","test_order",nbins,-10,200);

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		if(currentRays.size()>0){
			Ray first = currentRays[0];
			for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
				if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
					chiSquareH->Fill(evn,it->get_chiSquare_X()+it->get_chiSquare_Y());
					//sigma_t->Fill(it->get_t_sigma());
					if(it->get_t_mean()<first.get_t_mean()){
						first = Ray(*it);
					}
				}
			}
			test_order->Fill(first.get_chiSquare_X() + first.get_chiSquare_Y());
		}
		if(currentRays.size()>1){
			//diff_t->Fill(currentRays[0].get_t_mean() - currentRays[1].get_t_mean());
			double sigma0 = currentRays[0].get_t_sigma();
			double sigma1 = currentRays[1].get_t_sigma();
			double t0 = currentRays[0].get_t_mean();
			double t1 = currentRays[1].get_t_mean();
			if(t0<t1){
				tLO->Fill(t0);
				tNLO->Fill(t1);
				sigmaLO->Fill(sigma0);
				sigmaNLO->Fill(sigma1);
			}
			else{
				tLO->Fill(t1);
				tNLO->Fill(t0);
				sigmaLO->Fill(sigma1);
				sigmaNLO->Fill(sigma0);
			}
		}
		for(vector<Detector*>::iterator jt=detectors.begin();jt!=detectors.end();++jt){
			if((*jt)->get_type() == "CM"){
				CM_Detector * currentDet = dynamic_cast<CM_Detector*>(*jt);
				CM_residus[currentDet->get_cm_n_in_tree()]->Fill(evttime,currentCBEvent->get_clus_N_by_det(currentDet));
				CM_spark_h[currentDet->get_cm_n_in_tree()]->Fill(evttime,CM_Spark[currentDet->get_cm_n_in_tree()]);
			}
			if((*jt)->get_type() == "MG"){
				MG_Detector * currentDet = dynamic_cast<MG_Detector*>(*jt);
				MG_residus[currentDet->get_mg_n_in_tree()]->Fill(evttime,currentCBEvent->get_clus_N_by_det(currentDet));
				MG_spark_h[currentDet->get_mg_n_in_tree()]->Fill(evttime,MG_Spark[currentDet->get_mg_n_in_tree()]);
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
			for(int i=0;i<CM_N;i++){
				c_CM[i]->cd();
				CM_residus[i]->Draw();
				CM_spark_h[i]->Draw("SAME");
				c_CM[i]->Modified();
				c_CM[i]->Update();
			}
			for(int i=0;i<MG_N;i++){
				c_MG[i]->cd();
				MG_residus[i]->Draw();
				MG_spark_h[i]->Draw("SAME");
				c_MG[i]->Modified();
				c_MG[i]->Update();
			}
			c0->cd();
			chiSquareH->Draw();
			c0->Modified();
			c0->Update();
			c_t->cd();
			//diff_t->Draw();
			tLO->Draw();
			tNLO->Draw("SAME");
			c_t->Modified();
			c_t->Update();
			c_sigma->cd();
			//sigma_t->Draw();
			sigmaLO->Draw();
			sigmaNLO->Draw("SAME");
			c_sigma->Modified();
			c_sigma->Update();
			c_test->cd();
			test_order->Draw();
			c_test->Modified();
			c_test->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	for(int i=0;i<CM_N;i++){
		c_CM[i]->cd();
		CM_residus[i]->Draw();
		CM_spark_h[i]->Draw("SAME");
		c_CM[i]->Modified();
		c_CM[i]->Update();
	}
	for(int i=0;i<MG_N;i++){
		c_MG[i]->cd();
		MG_residus[i]->Draw();
		MG_spark_h[i]->Draw("SAME");
		c_MG[i]->Modified();
		c_MG[i]->Update();
	}
	c0->cd();
	chiSquareH->Draw();
	c0->Modified();
	c0->Update();
	c_t->cd();
	//diff_t->Draw();
	tLO->Draw();
	tNLO->Draw("SAME");
	c_t->Modified();
	c_t->Update();
	c_sigma->cd();
	//sigma_t->Draw();
	sigmaLO->Draw();
	sigmaNLO->Draw("SAME");
	c_sigma->Modified();
	c_sigma->Update();
	c_test->cd();
	test_order->Draw();
	c_test->Modified();
	c_test->Update();
}
TH2D * Analyse::AbsorptionFluxMap(double z, int nbins, TCanvas * c1){
	int eventReconstructed = 0;
	int eventSuitable = 0;

	//gStyle->SetPalette(1);
	//double z_Pb = 1553;
	if(c1 == 0){
		c1 = new TCanvas("fluxMap","fluxMap");
	}
	TH2D * fluxMapZ = new TH2D("fluxMapZ","fluxMapZ",nbins,-100,600,nbins,-100,600);
	fluxMapZ->SetStats(0);

	if (fChain == 0) return fluxMapZ;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
				fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
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
void Analyse::AbsorptionFluxMapNorm(double z,TH2D * background, int nbins, TCanvas * c1, TCanvas * c2, TCanvas * c3){
	int eventReconstructed = 0;
	int eventSuitable = 0;

	gStyle->SetPalette(55,0);
	//gStyle->SetPalette(1);
	//double z_Pb = 1553;
	if(c1==0) c1 = new TCanvas("fluxMap","fluxMap");
	TH2D * fluxMapZ = new TH2D("fluxMapSignal","fluxMapSignal",nbins,-100,600,nbins,-100,600);
	fluxMapZ->SetStats(0);
	if(c2 == 0) c2 = new TCanvas("fluxMapNorm","fluxMapNorm");
	TH2D * fluxMapSigma = new TH2D("fluxMapSigma","fluxMapSigma",nbins,-100,600,nbins,-100,600);
	fluxMapSigma->SetStats(0);
	if(c3 == 0) c3 = new TCanvas("fluxMap_Sigma","fluxMap_Sigma");

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays();
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			if(it->get_chiSquare_X()>-1 && it->get_chiSquare_Y()>-1){
				fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%10000 == 0){
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

	TCanvas * c1 = new TCanvas("thetaUp","thetaUp");
	TCanvas * c2 = new TCanvas("thetaDown","thetaDown");
	TH1D * thetaXUp = new TH1D("thetaXUp","thetaXUp",50,-2,2);
	thetaXUp->SetLineColor(2);
	TH1D * thetaYUp = new TH1D("thetaYUp","thetaYUp",50,-2,2);
	thetaYUp->SetLineColor(4);
	TH1D * thetaXDown = new TH1D("thetaXDown","thetaXDown",50,-2,2);
	thetaXDown->SetLineColor(2);
	TH1D * thetaYDown = new TH1D("thetaYDown","thetaYDown",50,-2,2);
	thetaYDown->SetLineColor(4);

	TCanvas * c4 = new TCanvas("correlation","correlation");
	TH2D * XY_correlation = new TH2D("correlation","correlation",70,-100,600,70,-100,600);

	TCanvas * c3 = new TCanvas("docas","docas");
	TH1D * doca = new TH1D("doca","doca",100,0,500);

	int eventReconstructed = 0;
	int eventSuitable = 0;

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

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
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
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		eventReconstructed+=currentCBEvent->get_rayPairs_N();
		for(unsigned int i=0;i<currentCBEvent->get_rayPairs_N();i++){
			RayPair currentRayPair = currentCBEvent->get_rayPair(i);
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
			if(PoCA.get_Z()>600 && PoCA.get_Z()<800) XY_correlation->Fill(PoCA.get_X(),PoCA.get_Y());
			outTree->Fill();

		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
			c1->cd();
			thetaXUp->Draw();
			thetaYUp->Draw("SAME");
			c1->Modified();
			c1->Update();
			c2->cd();
			thetaXDown->Draw();
			thetaYDown->Draw("SAME");
			c2->Modified();
			c2->Update();
			c3->cd();
			doca->Draw();
			c3->Modified();
			c3->Update();
			c4->cd();
			XY_correlation->Draw("COLZ");
			c4->Modified();
			c4->Update();
		}
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	outTree->Write();
	outFile->Close();
	c1->cd();
	thetaXUp->Draw();
	thetaYUp->Draw("SAME");
	c1->Modified();
	c1->Update();
	c2->cd();
	thetaXDown->Draw();
	thetaYDown->Draw("SAME");
	c2->Modified();
	c2->Update();
	c3->cd();
	doca->Draw();
	c3->Modified();
	c3->Update();
	c4->cd();
	XY_correlation->Draw("COLZ");
	c4->Modified();
	c4->Update();
}

/*void Analyse::MultiGenDebug(int i){
	int eventReconstructed = 0;
	int eventSuitable = 0;

	TCanvas * c1 = new TCanvas("MGPos","MGPos");
	TH1D * mgPos = new TH1D("mgPos","mgPos",1024,0,1023);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,-1);
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		for(vector<Event*>::iterator it=(currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if((*it)->get_type() == "MG" && (*it)->get_n_in_tree() == i){
				MG_Event currentMGEvent(*dynamic_cast<MG_Event*>(*it));
				for(vector<MG_Cluster>::iterator jt = currentMGEvent.clusters.begin();jt!=currentMGEvent.clusters.end();++jt){
					mgPos->Fill(jt->pos);
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0){
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

double Analyse::get_z_Up() const{
	double z_Up = numeric_limits<double>::max();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if((*it)->get_is_up() && (*it)->get_z()<z_Up){
			z_Up = (*it)->get_z();
		}
	}
	return z_Up;
}
double Analyse::get_z_Down() const{
	double z_Down = numeric_limits<double>::min();
	for(vector<Detector*>::const_iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_up()) && (*it)->get_z()>z_Down){
			z_Down = (*it)->get_z();
		}
	}
	return z_Down;
}

void Analyse::bugtest(){
	int before = 0;
	int afterAbsorption = 0;
	int afterDeviation = 0;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<  setw(20) << "before" <<  "|" << setw(20) << "after absorption" <<  "|" << setw(20) << "after deviation" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent1 = new CosmicBenchEvent(this,this,-1);
		currentCBEvent1->Demux_CM();
		before+=currentCBEvent1->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent1;
		CosmicBenchEvent * currentCBEvent2 = new CosmicBenchEvent(this,this,-1);
		currentCBEvent2->get_absorption_rays();
		afterAbsorption+=currentCBEvent2->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent2;
		CosmicBenchEvent * currentCBEvent3 = new CosmicBenchEvent(this,this,-1);
		currentCBEvent3->createPairs();
		afterDeviation+=currentCBEvent3->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent3;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << before << "|" << setw(20) << afterAbsorption << "|" << setw(20) << afterDeviation << flush;
	}
	cout << "\r"<< setw(20) << before << "|" << setw(20) << afterAbsorption << "|" << setw(20) << afterDeviation << endl;
}

void Analyse::CalcStripResponseFunction(int bin_nb){

	gStyle->SetPalette(55,0);
	gStyle->SetOptFit(0111);

	int det_N = CM_N + MG_N;
	cout << setw(20) << "detector proccessed" << "|" << setw(20) << "event processed" << endl;
	TProfile * SRH;
	//TGraph * SRH2D[det_N];
	TH2D * SRH2D;
	TCanvas * c;
	TCanvas * d;
	TF1 * SRF;

	TCanvas * c_coord[bin_nb];
	TProfile * SRH_coord[bin_nb];
	TF1 * SRF_coord[bin_nb];
	TCanvas * c_param;
	TGraphErrors * gauss_width_graph;
	TGraphErrors * lorentz_width_graph;
	TGraphErrors * offset_graph;

	double chisquare_threshold = 10;
	Long64_t nentries = fChain->GetEntriesFast();

	TFile * signal_file = new TFile(signal_file_name.c_str(),"READ");
	TTree * signal_tree = (TTree*)(signal_file->Get("T"));

	if(nentries != signal_tree->GetEntriesFast()){
		cout << "total number of event in signal and analyse tree does not match" << endl;
		return;
	}

	float StripAmpl_MG_corr[MG_N][61][32];
	float StripAmpl_CM_corr[CM_N][64][32];
	int signal_evn;
	signal_tree->SetBranchAddress("Nevent",&signal_evn);
	if(CM_N>0) signal_tree->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
	if(MG_N>0) signal_tree->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);

	int det_x_n = 0;
	int det_y_n = 0;
	int nref_nb = 0;
	for(int j=0;j<det_N;j++){
		if(detectors[j]->get_is_X()) det_x_n++;
		else det_y_n++;
		if(!(detectors[j]->is_ref)) nref_nb++;
	}
	if(nref_nb>1){
		cout << "there is more than 1 non ref det, exiting" << endl;
		return;
	}

	for(int i=0;i<det_N;i++){
		if(detectors[i]->is_ref) continue;
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
		if(detectors[i]->get_type() == "MG" && detectors[i]->get_is_X() == false) limit *= 2.;
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
		for (Long64_t jentry=0; jentry<nentries;jentry++){
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
			vector<Ray> currentRays = currentCBEvent->get_absorption_rays();

			for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
				if(!((*it)->get_is_ref())){
					if((*it)->get_type() == "CM_Demux"){
						continue;
					}
					else if((*it)->get_type() == "MG"){
						vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*it))->get_clusters();
						for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){

							double chiSquare_in_nref_dir = (detectors[i]->get_is_X()) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
							int clus_in_nref_dir = (detectors[i]->get_is_X()) ? jt->get_clus_x_n() : jt->get_clus_y_n();
							if(clus_in_nref_dir<(det_in_nref_dir-1)) continue;
							if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;

							if((jt->get_chiSquare_X()+jt->get_chiSquare_Y()) > chisquare_threshold) continue;
							double residu = numeric_limits<double>::max();
							vector<MG_Cluster>::iterator matching_cluster = current_clusters.end();
							for(vector<MG_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
								if(kt->get_is_X()){
									kt->set_perp_pos_mm(jt->eval_Y((*it)->get_z()));
								}
								else{
									kt->set_perp_pos_mm(jt->eval_X((*it)->get_z()));
								}
								double current_residu = jt->get_residu_ref(&(*kt));
								if(current_residu<residu){
									residu = current_residu;
									matching_cluster = kt;
								}
							}
							if(matching_cluster == current_clusters.end()) continue;
							
							double normalization = matching_cluster->get_ampl()/matching_cluster->get_size();
							double matching_position = (matching_cluster->get_is_X()) ? jt->eval_X((*it)->get_z()) : jt->eval_Y((*it)->get_z());
							//double matching_position = matching_cluster->get_pos_mm();
							/*
							for(int strip_nb = matching_cluster->get_pos()-1;strip_nb<(matching_cluster->get_pos()+2);strip_nb++){
								int channel = MG_Detector::StripToChannel(strip_nb);
								double current_max_ampl = *max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+32);
								if(current_max_ampl>normalization) normalization = current_max_ampl;
							}
							*/
							for(int strip_nb=0;strip_nb<1024;strip_nb++){
								int channel = MG_Detector::StripToChannel(strip_nb);
								if(Abs(residu)<50.){
									SRH->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+32))/normalization);
									if(bin_nb>0) SRH_coord[FloorNint(matching_position*bin_nb/500.)]->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+32))/normalization);
								}
								//SRH2D[i]->SetPoint(graph_point_nb, matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+32))/normalization);
								//graph_point_nb++;
								SRH2D->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+32))/normalization);
							}
							current_clusters.erase(matching_cluster);
						}
					}
				}
			}

			delete currentCBEvent;

			if(jentry%500 == 0) cout << "\r" << setw(20) << i << "|" << setw(20) << jentry << flush;
			if(jentry%10000 == 0){
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
			gauss_width_graph->SetPoint(i_coord,(250./bin_nb)+(i_coord*500./bin_nb),SRF_coord[i_coord]->GetParameter(1));
			gauss_width_graph->SetPointError(i_coord,250./bin_nb,SRF_coord[i_coord]->GetParError(1));
			lorentz_width_graph->SetPoint(i_coord,(250./bin_nb)+(i_coord*500./bin_nb),SRF_coord[i_coord]->GetParameter(2));
			lorentz_width_graph->SetPointError(i_coord,250./bin_nb,SRF_coord[i_coord]->GetParError(2));
			offset_graph->SetPoint(i_coord,(250./bin_nb)+(i_coord*500./bin_nb),SRF_coord[i_coord]->GetParameter(3));
			offset_graph->SetPointError(i_coord,250./bin_nb,SRF_coord[i_coord]->GetParError(3));
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

void Analyse::EventDisplay(int event_nb){
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

	float StripAmpl_MG_corr[MG_N][61][32];
	float StripAmpl_CM_corr[CM_N][64][32];
	int signal_evn;
	signal_tree->SetBranchAddress("Nevent",&signal_evn);
	if(CM_N>0) signal_tree->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
	if(MG_N>0) signal_tree->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);

	double chisquare_threshold = 50;
	int detector_color = 1;
	int ray_color = 2;
	int pos_color = 4;

	LoadTree(event_nb);
	GetEntry(event_nb);
	signal_tree->LoadTree(event_nb);
	signal_tree->GetEntry(event_nb);
	CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,-1);
	vector<Ray> eventRays = CBEvent->get_absorption_rays();
	delete CBEvent;
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
	if(eventRays.size() == 0){
		cout << "no ray in this event, exiting" << endl;
		return;
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
		cout << rays_it->eval_X(min_z) << endl;
		cout << rays_it->eval_Y(min_z) << endl;
		ray_x.SetNextPoint(rays_it->eval_X(min_z),min_z,0);
		ray_y.SetNextPoint(rays_it->eval_Y(min_z),min_z,0);
		ray_x.SetNextPoint(rays_it->eval_X(max_z),max_z,0);
		ray_y.SetNextPoint(rays_it->eval_Y(max_z),max_z,0);
		ray_x.SetLineColor(ray_color);
		ray_y.SetLineColor(ray_color);
		ray_lines.push_back(pair<TPolyLine3D*, TPolyLine3D*>(new TPolyLine3D(ray_x),new TPolyLine3D(ray_y)));
		vector<Cluster*> ray_clusters = rays_it->get_clus();
		for(vector<Cluster*>::iterator clus_it = ray_clusters.begin();clus_it!=ray_clusters.end();++clus_it){
			//add clus_pos line
			double det_z = (*clus_it)->get_z();
			double current_clus_size = (*clus_it)->get_size();
			double current_clus_pos = (*clus_it)->get_pos_mm();
			double min_ampl = 0;
			double max_ampl = 0;
			if((*clus_it)->get_type() == "MG"){
				int det_n_in_tree = dynamic_cast<MG_Detector*>(detectors[(*clus_it)->find_det(detectors)])->get_mg_n_in_tree();
				for(int i_strip=0;i_strip<1024;i_strip++){
					int channel = MG_Detector::StripToChannel(i_strip);
					double current_strip_pos = (*clus_it)->correct_strip_nb(i_strip);
					if(current_strip_pos>current_clus_pos+(2*current_clus_size)) continue;
					if(current_strip_pos<current_clus_pos-(2*current_clus_size)) continue;
					double current_ampl = *max_element(StripAmpl_MG_corr[det_n_in_tree][channel],StripAmpl_MG_corr[det_n_in_tree][channel]+32);
					if(current_ampl>max_ampl) max_ampl = current_ampl;
					if(current_ampl<min_ampl) min_ampl = current_ampl;
					if((*clus_it)->get_is_X()){
						ampl_graph_x->SetPoint(ampl_graph_x_point_nb,current_strip_pos,det_z,current_ampl);
						ampl_graph_x_point_nb++;
					}
					else{
						ampl_graph_y->SetPoint(ampl_graph_y_point_nb,current_strip_pos,det_z,current_ampl);
						ampl_graph_y_point_nb++;
					}

				}
			}
			else if((*clus_it)->get_type() == "CM" || (*clus_it)->get_type() == "CM_Demux"){
				//TODO : implement for CM
			}
			min_ampl = min_ampl - 0.1*(max_ampl - min_ampl);
			max_ampl = max_ampl + 0.1*(max_ampl - min_ampl);
			TPolyLine3D clus_line(2);
			clus_line.SetNextPoint(current_clus_pos,det_z,min_ampl);
			clus_line.SetNextPoint(current_clus_pos,det_z,max_ampl);
			clus_line.SetLineColor(pos_color);
			clus_line.SetLineStyle(2);
			if((*clus_it)->get_is_X()) pos_lines_x.push_back(new TPolyLine3D(clus_line));
			else pos_lines_y.push_back(new TPolyLine3D(clus_line));
		}
	}
	delete signal_file;
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
		cout << "blah" << endl;
		c_visu->cd(1);
		(pair_it->first)->Draw();
		c_visu->cd(2);
		(pair_it->second)->Draw();
	}
	c_visu->Modified();
	c_visu->Update();
}