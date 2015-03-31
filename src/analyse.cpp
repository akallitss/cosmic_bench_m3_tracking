#define analyse_cpp
//my class
#include "analyse.h"
#include "T.h"
#include "detector.h"
#include "ray.h"
#include "event.h"
#include "Tsignal.h"
#include "tomography.h"
#include "acceptanceFunction.h"
#include "Tray.h"
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
using TMath::Max;
using TMath::Min;
using TMath::Pi;

Analyse::Analyse(string configFilePath){
	ptree config_tree;
	read_json(configFilePath, config_tree);
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CosmicBench::Init(config_tree);
	T::Init(tree,CM_N,MG_N);
	signal_file_name = config_tree.get<string>("signal_file");
}
Analyse::Analyse(ptree config_tree){
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CosmicBench::Init(config_tree);
	T::Init(tree,CM_N,MG_N);
	signal_file_name = config_tree.get<string>("signal_file");
}
Analyse::~Analyse(){
	//delete f;
}
void Analyse::Residus(){
	TCanvas * c_CM[CM_N];
	TCanvas * c_MG[MG_N];
	TH1D * CM_residus[CM_N];
	TH1D * MG_residus[MG_N];
	int nbins = 200;
	long eventReconstructed = 0;
	long eventSuitable = 0;
	for(int i=0;i<CM_N;i++){
		ostringstream name;
		name << "Cosmulti_" << i;
		c_CM[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		name << "_residu";
		CM_residus[i] = new TH1D(name.str().c_str(),name.str().c_str(),nbins,-Tomography::XY_size,Tomography::XY_size);
	}
	for(int i=0;i<MG_N;i++){
		ostringstream name;
		name << "MultiGen_" << i;
		c_MG[i] = new TCanvas(name.str().c_str(),name.str().c_str());
		name << "_residu";
		MG_residus[i] = new TH1D(name.str().c_str(),name.str().c_str(),nbins,-Tomography::XY_size,Tomography::XY_size);
	}
	TCanvas * c0 = new TCanvas("chiSquare","chiSquare");
	TH1D * chiSquareH = new TH1D("chiSquares","chiSquares",nbins,0,3*Tomography::XY_size);

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
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
					if((*jt)->get_type() == Tomography::CM){
						CM_Detector * currentDet = dynamic_cast<CM_Detector*>(*jt);
						CM_residus[i_CM]->Fill(it->get_residu(currentDet));
						i_CM++;
					}
					if((*jt)->get_type() == Tomography::MG){
						MG_Detector * currentDet = dynamic_cast<MG_Detector*>(*jt);
						MG_residus[i_MG]->Fill(it->get_residu(currentDet));
						i_MG++;
					}
				}
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
	gStyle->SetNumberContours(512);
	map<string,TCanvas*> c_MM;
	map<string,TH1D*> MM_residus;
	map<string,TF1*> offset_fit;
	map<string,TH2D*> muon_seen;
	map<string,TH2D*> muon_total;
	map<string,TH2D*> efficacity_2D;
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
	map<int,int> perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			if((*it)->get_type() == Tomography::CM){
				name << "Cosmulti_" << (dynamic_cast<CM_Detector*>(*it))->get_cm_n_in_tree();
			}
			else if((*it)->get_type() == Tomography::MG){
				name << "Multigen_" << (dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree();
				if((*it)->get_perp_n()>-1) perp_pairs[(dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree()] = (*it)->get_perp_n();
			}
			c_MM[name.str()] = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
			c_MM[name.str()]->Divide(4,4);
			MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
			muon_seen[name.str()] = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
			muon_total[name.str()] = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
			efficacity_2D[name.str()] = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
			efficacity_2D[name.str()]->SetStats(false);
			correlation[name.str()] = new TGraph();
			angle_alignment[name.str()] = new TProfile((name.str()+"_angle").c_str(),(name.str()+"_angle").c_str(),500,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,-5,5);
			resVSpos[name.str()] = new TProfile((name.str()+"_resVSpos").c_str(),(name.str()+"_resVSpos").c_str(),500,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,-5,5);
			resVSampl[name.str()] = new TProfile((name.str()+"_resVSampl").c_str(),(name.str()+"_resVSampl").c_str(),500,-100,10000,-5,5);
			resVStime[name.str()] = new TProfile((name.str()+"_resVStime").c_str(),(name.str()+"_resVStime").c_str(),38,-2,34,-5,5);
			resVSangle[name.str()] = new TProfile((name.str()+"_resVSangle").c_str(),(name.str()+"_resVSangle").c_str(),50,-0.6,0.6,-5,5);
			resVSanglePerp[name.str()] = new TProfile((name.str()+"_resVSanglePerp").c_str(),(name.str()+"_resVSanglePerp").c_str(),50,-0.6,0.6,-5,5);
			resVStot[name.str()] = new TProfile((name.str()+"_resVStot").c_str(),(name.str()+"_resVStot").c_str(),26,0,25,-5,5);
			resVSsize[name.str()] = new TProfile((name.str()+"_resVSsize").c_str(),(name.str()+"_resVSsize").c_str(),50,0,50,-5,5);
			absResVSampl[name.str()] = new TProfile((name.str()+"_absResVSampl").c_str(),(name.str()+"_absResVSampl").c_str(),500,-100,10000,0,5);
			absResVStime[name.str()] = new TProfile((name.str()+"_absResVStime").c_str(),(name.str()+"_absResVStime").c_str(),38,-2,34,0,5);
			absResVSabsAngle[name.str()] = new TProfile((name.str()+"_absResVSangle").c_str(),(name.str()+"_absResVSabsAngle").c_str(),50,0,0.6,0,5);
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
	for(map<int,int>::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
		}
	}
	if(non_ref_n>2){
		cout << "too many non ref det" << endl;
	}
	if(nref_x_n>1){
		cout << "too many non ref det" << endl;
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2,2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",MG_N + CM_N + 2,0,MG_N + CM_N + 2);
	TH1D * ray_slope = new TH1D("slope","slope",100,0,1);
	TH1D * ray_phi = new TH1D("phi","phi",100,-Pi(),Pi());
	TH1D * ray_slope_X = new TH1D("slope_X","slope_X",100,0,1);
	TH1D * ray_slope_Y = new TH1D("slope_Y","slope_Y",100,0,1);
	ray_slope->SetLineColor(1);
	ray_slope_X->SetLineColor(2);
	ray_slope_Y->SetLineColor(3);
	if (fChain == 0) return;
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
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
		eventSuitable+=currentCBEvent->get_clus_N()*1./(CM_N+MG_N);
		
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				if((*it)->get_type() == Tomography::CM_Demux){
					ostringstream name;
					name << "Cosmulti_" << (*it)->get_n_in_tree();
					vector<CM_Demux_Cluster> current_clusters = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
						if(jt->get_clus_n()<static_cast<unsigned int>(CM_N+MG_N-non_ref_n)) continue;
						double residu = numeric_limits<double>::max();
						vector<CM_Demux_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<CM_Demux_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							kt->set_perp_pos_mm(*jt);
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						if(matching_cluster == current_clusters.end()) continue;
						matching_cluster->set_perp_pos_mm(*jt);
						if(matching_cluster->get_is_X()){
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),matching_cluster->get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSpos[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
							absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_X()),Abs(residu));
							resVSanglePerp[name.str()]->Fill(jt->get_slope_Y(),residu);
						}
						else{
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),matching_cluster->get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSpos[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
							absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_Y()),Abs(residu));
							resVSanglePerp[name.str()]->Fill(jt->get_slope_X(),residu);
						}
						point_nb[name.str()]++;
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
						resVStime[name.str()]->Fill(matching_cluster->get_t(),residu);
						resVSampl[name.str()]->Fill(matching_cluster->get_ampl(),residu);
						resVStot[name.str()]->Fill(matching_cluster->get_TOT(),residu);
						resVSsize[name.str()]->Fill(matching_cluster->get_size(),residu);
						absResVStime[name.str()]->Fill(matching_cluster->get_t(),Abs(residu));
						absResVSampl[name.str()]->Fill(matching_cluster->get_ampl(),Abs(residu));
						absResVStot[name.str()]->Fill(matching_cluster->get_TOT(),Abs(residu));
						absResVSsize[name.str()]->Fill(matching_cluster->get_size(),Abs(residu));
						if(residu<chisquare_threshold){
							muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						}
					}
				}
				else if((*it)->get_type() == Tomography::MG){
					ostringstream name;
					name << "Multigen_" << (*it)->get_n_in_tree();
					vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
						if(jt->get_clus_n()<static_cast<unsigned int>(CM_N+MG_N-non_ref_n)) continue;
						double residu = numeric_limits<double>::max();
						vector<MG_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<MG_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							kt->set_perp_pos_mm(*jt);
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						muon_total[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						if(matching_cluster == current_clusters.end()) continue;
						matching_cluster->set_perp_pos_mm(*jt);
						if(matching_cluster->get_is_X()){
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_X((*it)->get_z()),matching_cluster->get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSpos[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
							absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_X()),Abs(residu));
							resVSanglePerp[name.str()]->Fill(jt->get_slope_Y(),residu);
						}
						else{
							correlation[name.str()]->SetPoint(point_nb[name.str()],jt->eval_Y((*it)->get_z()),matching_cluster->get_pos_mm());
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSpos[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
							absResVSabsAngle[name.str()]->Fill(Abs(jt->get_slope_Y()),Abs(residu));
							resVSanglePerp[name.str()]->Fill(jt->get_slope_X(),residu);
						}
						point_nb[name.str()]++;
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
						resVStime[name.str()]->Fill(matching_cluster->get_t(),residu);
						resVSampl[name.str()]->Fill(matching_cluster->get_ampl(),residu);
						resVStot[name.str()]->Fill(matching_cluster->get_TOT(),residu);
						resVSsize[name.str()]->Fill(matching_cluster->get_size(),residu);
						absResVStime[name.str()]->Fill(matching_cluster->get_t(),Abs(residu));
						absResVSampl[name.str()]->Fill(matching_cluster->get_ampl(),Abs(residu));
						absResVStot[name.str()]->Fill(matching_cluster->get_TOT(),Abs(residu));
						absResVSsize[name.str()]->Fill(matching_cluster->get_size(),Abs(residu));
						if(residu<chisquare_threshold){
							muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						}
					}
				}
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
				double pos_X = muon_total[it->first]->GetXaxis()->GetBinCenter(i);
				double pos_Y = muon_total[it->first]->GetYaxis()->GetBinCenter(j);
				if(pos_X<=2*Tomography::XY_size/5. && pos_X>=-2*Tomography::XY_size/5. && pos_Y<=2*Tomography::XY_size/5. && pos_Y>=-2*Tomography::XY_size/5.){
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
double Analyse::Residus_ref_cost(){
	int non_ref_n = 0;
	double chisquare_threshold = 10;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
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
	map<int,int> perp_pairs;
	//map<string, unsigned int> det_in_nref_dir;
	//map<string, bool> nref_is_X;
	//unsigned int det_x_n = 0;
	unsigned int nref_x_n = 0;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			ostringstream name;
			if((*it)->get_type() == Tomography::CM){
				name << "Cosmulti_" << (dynamic_cast<CM_Detector*>(*it))->get_cm_n_in_tree();
			}
			else if((*it)->get_type() == Tomography::MG){
				name << "Multigen_" << (dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree();
				if((*it)->get_perp_n()>-1) perp_pairs[(dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree()] = (*it)->get_perp_n();
			}
			MM_residus[name.str()] = new TH1D((name.str()+"_residu").c_str(),(name.str()+"_residu").c_str(),nbins,-5,5);
			angle_alignment[name.str()] = new TProfile((name.str()+"_angle").c_str(),(name.str()+"_angle").c_str(),500,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,-5,5);
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
	for(map<int,int>::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
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
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
				++ray_it;
			}
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(CM_N+MG_N);
		
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				if((*it)->get_type() == Tomography::CM_Demux){
					ostringstream name;
					name << "Cosmulti_" << (*it)->get_n_in_tree();
					vector<CM_Demux_Cluster> current_clusters = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
						if(jt->get_clus_n()<static_cast<unsigned int>(CM_N+MG_N-non_ref_n)) continue;
						double residu = numeric_limits<double>::max();
						vector<CM_Demux_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<CM_Demux_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							kt->set_perp_pos_mm(*jt);
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						if(matching_cluster == current_clusters.end()) continue;
						matching_cluster->set_perp_pos_mm(*jt);
						if(matching_cluster->get_is_X()){
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
						}
						else{
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
						}
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
					}
				}
				else if((*it)->get_type() == Tomography::MG){
					ostringstream name;
					name << "Multigen_" << (*it)->get_n_in_tree();
					vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*it))->get_clusters();
					for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
						//double chiSquare_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_chiSquare_X() : jt->get_chiSquare_Y();
						//unsigned int clus_in_nref_dir = (nref_is_X[name.str()]) ? jt->get_clus_x_n() : jt->get_clus_y_n();
						//if(clus_in_nref_dir<det_in_nref_dir[name.str()]) continue;
						//if(chiSquare_in_nref_dir > chisquare_threshold/static_cast<double>(clus_in_nref_dir)) continue;
						if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(non_ref_n)) continue;
						if(jt->get_clus_n()<static_cast<unsigned int>(CM_N+MG_N-non_ref_n)) continue;
						double residu = numeric_limits<double>::max();
						vector<MG_Cluster>::iterator matching_cluster = current_clusters.end();
						for(vector<MG_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
							kt->set_perp_pos_mm(*jt);
							double current_residu = jt->get_residu_ref(&(*kt));
							if(current_residu<residu){
								residu = current_residu;
								matching_cluster = kt;
							}
						}
						if(matching_cluster == current_clusters.end()) continue;
						matching_cluster->set_perp_pos_mm(*jt);
						if(matching_cluster->get_is_X()){
							angle_alignment[name.str()]->Fill(jt->eval_Y((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_X(),residu);
						}
						else{
							angle_alignment[name.str()]->Fill(jt->eval_X((*it)->get_z()),residu);
							resVSangle[name.str()]->Fill(jt->get_slope_Y(),residu);
						}
						current_clusters.erase(matching_cluster);
						MM_residus[name.str()]->Fill(residu);
					}
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
	int non_ref_n = 0;
	double chisquare_threshold = 10;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())) non_ref_n++;
	}
	if(non_ref_n>2){
		cout << "too much non ref det" << endl;
		return;
	}
	if(non_ref_n<2){
		cout << "can't do 2D efficacity without at least 2 non ref det" << endl;
		return;
	}
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
	unsigned int det_n = CM_N + MG_N;
	vector<double> det_z;
	map<int,int> perp_pairs;
	ostringstream name;
	for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
		if(!((*it)->get_is_ref())){
			if((*it)->get_type() == Tomography::CM){
				if(name.str().size()>0) name << "_";
				name << "Cosmulti_" << (dynamic_cast<CM_Detector*>(*it))->get_cm_n_in_tree();
			}
			else if((*it)->get_type() == Tomography::MG){
				if(name.str().size()>0) name << "_";
				name << "Multigen_" << (dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree();
				if((*it)->get_perp_n()>-1) perp_pairs[(dynamic_cast<MG_Detector*>(*it))->get_mg_n_in_tree()] = (*it)->get_perp_n();
			}
			if((*it)->get_is_X()) nref_x_n++;
			det_z.push_back((*it)->get_z());
		}
		if((*it)->get_is_X()) det_x_n++;
	}
	if(det_z.size()!=2){
		cout << "problem in nref det number" << endl;
		return;
	}
	if(det_z[0]!=det_z[1]){
		cout << "you can only calculate 2D efficacity for 2D detector (which have same z)" << endl;
		return;
	}
	for(map<int,int>::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0){
			cout << "you can only calculate 2D efficacity for 2D detector" << endl;
			return;
		}
	}
	if(nref_x_n!=1){
		cout << "2D efficacity can only be done with one det in each direction" << endl;
		return;
	}

	c_MM = new TCanvas(name.str().c_str(),name.str().c_str(),1200,1000);
	muon_seen = new TH2D((name.str()+"_seen").c_str(),(name.str()+"_seen").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
	muon_total = new TH2D((name.str()+"_total").c_str(),(name.str()+"_total").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
	efficacity_2D = new TH2D((name.str()+"_efficacity").c_str(),(name.str()+"_efficacity").c_str(),nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2,nbins_2D,-(1+marge)*Tomography::XY_size/2,(1+marge)*Tomography::XY_size/2);
	efficacity_2D->SetStats(false);

	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2);
	TH1D * chisquares = new TH1D("chiSquares","chiSquares",nbins,0,chisquare_threshold);
	TH1D * ray_clus_n = new TH1D("clus_n","clus_n",MG_N + CM_N + 2,0,MG_N + CM_N + 2);

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
			if((jt->get_chiSquare_X()+jt->get_chiSquare_Y()) < chisquare_threshold){
				chisquares->Fill(jt->get_chiSquare_X()+jt->get_chiSquare_Y());
				ray_clus_n->Fill(jt->get_clus_n());
			}
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()*1./(CM_N+MG_N);
		
		vector<Event*> nref_event;
		for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
			if(!((*it)->get_is_ref())){
				if((*it)->get_type() == Tomography::CM_Demux){
					nref_event.push_back(new CM_Demux_Event(*dynamic_cast<CM_Demux_Event*>(*it)));
				}
				else if((*it)->get_type() == Tomography::MG){
					nref_event.push_back(new MG_Event(*dynamic_cast<MG_Event*>(*it)));
				}
			}
		}
		if(nref_event.size()!=2){
			cout << "problem in event size" << endl;
			return;
		}
		delete currentCBEvent;
		for(vector<Ray>::iterator jt=currentRays.begin();jt!=currentRays.end();++jt){
			if(jt->get_clus_n()<(det_n-2)){
				continue;
			}
			if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(det_n-2)){
				continue;
			}
			vector<unsigned int> seen_clus_in_array;
			for(vector<Event*>::iterator it = nref_event.begin();it!=nref_event.end();++it){
				if((*it)->get_type() == Tomography::CM_Demux){
					vector<CM_Demux_Cluster> current_clusters = (dynamic_cast<CM_Demux_Event*>(*it))->get_clusters();
					double residu = numeric_limits<double>::max();
					vector<CM_Demux_Cluster>::iterator matching_cluster = current_clusters.end();
					for(vector<CM_Demux_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						kt->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(&(*kt));
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					if(matching_cluster == current_clusters.end()) continue;
					matching_cluster->set_perp_pos_mm(*jt);
					if(residu<chisquare_threshold){
						seen_clus_in_array.push_back(matching_cluster - current_clusters.begin());
					}
				}
				else if((*it)->get_type() == Tomography::MG){
					vector<MG_Cluster> current_clusters = (dynamic_cast<MG_Event*>(*it))->get_clusters();
					double residu = numeric_limits<double>::max();
					vector<MG_Cluster>::iterator matching_cluster = current_clusters.end();
					for(vector<MG_Cluster>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						kt->set_perp_pos_mm(*jt);
						double current_residu = jt->get_residu_ref(&(*kt));
						if(current_residu<residu){
							residu = current_residu;
							matching_cluster = kt;
						}
					}
					if(matching_cluster == current_clusters.end()) continue;
					matching_cluster->set_perp_pos_mm(*jt);
					if(residu<chisquare_threshold){
						seen_clus_in_array.push_back(matching_cluster - current_clusters.begin());
					}
				}
			}
			muon_total->Fill(jt->eval_X(det_z[0]),jt->eval_Y(det_z[0]));
			if(seen_clus_in_array.size()==2){
				muon_seen->Fill(jt->eval_X(det_z[0]),jt->eval_Y(det_z[0]));
				for(int i_event=0;i_event<2;i_event++){
					if(nref_event[i_event]->get_type() == Tomography::CM_Demux){
						CM_Demux_Event * current_event = static_cast<CM_Demux_Event*>(nref_event[i_event]);
						(current_event->clusters).erase((current_event->clusters).begin()+seen_clus_in_array[i_event]);
					}
					else if(nref_event[i_event]->get_type() == Tomography::MG){
						MG_Event * current_event = static_cast<MG_Event*>(nref_event[i_event]);
						(current_event->clusters).erase((current_event->clusters).begin()+seen_clus_in_array[i_event]);
					}
				}
			}
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << static_cast<long>(eventSuitable) << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
			for(int i=1;i<=nbins_2D;i++){
				for(int j=1;j<=nbins_2D;j++){
					int binN = muon_total->GetBin(i,j);
					double binContent = 0;
					if(muon_total->GetBinContent(binN) > 0) binContent = (muon_seen->GetBinContent(binN))/(muon_total->GetBinContent(binN));
					efficacity_2D->SetBinContent(binN,binContent);
				}
			}
			c_MM->cd();
			efficacity_2D->Draw("COLZ");
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
			if(pos_X<=2*Tomography::XY_size/5. && pos_X>=-2*Tomography::XY_size/5. && pos_Y<=2*Tomography::XY_size/5. && pos_Y>=-2*Tomography::XY_size/5.){
				total_seen += muon_seen->GetBinContent(binN);
				total_passed += muon_total->GetBinContent(binN);
			}
		}
	}
	efficacity = total_seen/total_passed;
	c_MM->cd();
	efficacity_2D->Draw("COLZ");
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
	TCanvas * c_CM[CM_N];
	TCanvas * c_MG[MG_N];
	TProfile * CM_residus[CM_N];
	TProfile * CM_spark_h[CM_N];
	TProfile * MG_residus[MG_N];
	TProfile * MG_spark_h[MG_N];
	int nbins = 200;
	//int lim = 100;
	long eventReconstructed = 0;
	long eventSuitable = 0;
	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
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
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
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
			if((*jt)->get_type() == Tomography::CM){
				CM_Detector * currentDet = dynamic_cast<CM_Detector*>(*jt);
				CM_residus[currentDet->get_cm_n_in_tree()]->Fill(evttime,currentCBEvent->get_clus_N_by_det(currentDet));
				CM_spark_h[currentDet->get_cm_n_in_tree()]->Fill(evttime,CM_Spark[currentDet->get_cm_n_in_tree()]);
			}
			if((*jt)->get_type() == Tomography::MG){
				MG_Detector * currentDet = dynamic_cast<MG_Detector*>(*jt);
				MG_residus[currentDet->get_mg_n_in_tree()]->Fill(evttime,currentCBEvent->get_clus_N_by_det(currentDet));
				MG_spark_h[currentDet->get_mg_n_in_tree()]->Fill(evttime,MG_Spark[currentDet->get_mg_n_in_tree()]);
			}
		}
		delete currentCBEvent;
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		rayTree->fillTree(evn,evttime,currentRays,Z_Up,Z_Down);
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
	}
	cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << nentries << endl;
	rayTree->Write();
	rayTree->CloseFile();
}
TH2D * Analyse::AbsorptionFluxMap(double z, TCanvas * c1){
	long eventReconstructed = 0;
	long eventSuitable = 0;
	double chisquare_threshold = 100;

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
	double x_min = -Tomography::XY_size/2.;
	double x_max = Tomography::XY_size/2.;
	if(z>z_max){
		x_min -= Tomography::XY_size*(z - z_max)/(z_max - z_min);
		x_max = - x_min;
	}
	else if(z<z_min){
		x_min -= Tomography::XY_size*(z_min - z)/(z_max - z_min);
		x_max = - x_min;
	}
	double width = x_max - x_min;
	x_min -= 0.05*width;
	x_max += 0.05*width;

	//if (fChain == 0) return fluxMapZ;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	TH2D * fluxMapZ = new TH2D("fluxMapZ","fluxMapZ",Sqrt(0.02*nentries),x_min,x_max,Sqrt(0.02*nentries),x_min,x_max);
	fluxMapZ->SetStats(0);

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
void Analyse::AbsorptionFluxMapNormTheo(double z, TCanvas * c1, TCanvas * c2, TCanvas * c3, TCanvas * c4){
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
	double x_min = -Tomography::XY_size/2.;
	double x_max = Tomography::XY_size/2.;
	if(z>z_max){
		x_min -= Tomography::XY_size*(z - z_max)/(z_max - z_min);
		x_max = - x_min;
	}
	else if(z<z_min){
		x_min -= Tomography::XY_size*(z_min - z)/(z_max - z_min);
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
	acceptanceFunction acceptanceEstimation(0,500,0,500,z_max,z_min,0);
	TH2D * background = new TH2D(acceptanceEstimation.plot_XY(nbins,x_min,x_max,nbins,x_min,x_max,z));
	if(c4 == 0) c4 = new TCanvas("fluxMap_background","fluxMap_background");
	c4->cd();
	background->Draw("COLZ");
	c4->Modified();
	c4->Update();

	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
		vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
		vector<Ray>::iterator ray_it = currentRays.begin();
		while(ray_it!= currentRays.end()){
			if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold) ++ray_it;
			else ray_it = currentRays.erase(ray_it);
		}
		eventReconstructed+=currentRays.size();
		eventSuitable+=currentCBEvent->get_clus_N()/(CM_N+MG_N);
		delete currentCBEvent;
		for(vector<Ray>::iterator it=currentRays.begin();it!=currentRays.end();++it){
			fluxMapZ->Fill(it->eval_X(z),it->eval_Y(z));
		}
		if(jentry%500 == 0) cout << "\r"<< setw(20) << eventReconstructed << "|" << setw(20) << eventSuitable << "|" << setw(20) << jentry << flush;
		if(jentry%10000 == 0 && Tomography::live_graphic_display){
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
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
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
		if(jentry%10000 == 0 && Tomography::live_graphic_display){
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
	TH2D * XY_correlation = new TH2D("correlation","correlation",70,-6*Tomography::XY_size/10.,6*Tomography::XY_size/10.,70,-6*Tomography::XY_size/10.,6*Tomography::XY_size/10.);

	TCanvas * c3 = new TCanvas("docas","docas");
	TH1D * doca = new TH1D("doca","doca",100,0,Tomography::XY_size);

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

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
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

		CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
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
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
	outFile->cd();
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
	long eventReconstructed = 0;
	long eventSuitable = 0;

	TCanvas * c1 = new TCanvas("MGPos","MGPos");
	TH1D * mgPos = new TH1D("mgPos","mgPos",1024,0,1023);

	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	cout <<  setw(20) << "rays" <<  "|" << setw(20) << "suitable" <<  "|" << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
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
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
	if (fChain == 0) return;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	int limit = 10;
	if(limit<nentries) nentries = limit;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		cout << evn << " : " << endl;
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,false,-1);
		for(vector<Event*>::iterator it = (CBEvent->events).begin();it!=(CBEvent->events).end();++it){
			if((*it)->get_type() == Tomography::MG){
				vector<MG_Cluster> current_clusters = dynamic_cast<MG_Event*>(*it)->get_clusters();
				if(current_clusters.size()>0) cout << setw(10) << current_clusters.front().get_pos() << " | ";
				else cout << "No Cluster" << " | ";
			}
		}
		cout << endl;
		delete CBEvent;
		for(int i=0;i<MG_N-1;i++){
			cout << setw(10) << MG_ClusPos[i][0] << " | ";
		}
		cout << setw(10) << MG_ClusPos[MG_N-1][0] << endl;
	}
}

void Analyse::CalcStripResponseFunction(int bin_nb){

	gStyle->SetPalette(55,0);
	gStyle->SetNumberContours(512);
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

	TFile * signal_file = new TFile(signal_file_name.c_str(),"READ");
	TTree * signal_tree = (TTree*)(signal_file->Get("T"));

	if(fChain->GetEntriesFast() != signal_tree->GetEntriesFast()){
		cout << "total number of event in signal and analyse tree does not match" << endl;
		return;
	}

	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();

	float StripAmpl_MG_corr[MG_N][MG_Detector::Nstrip][Tomography::Nsample];
	float StripAmpl_CM_corr[CM_N][CM_Detector::Nstrip][Tomography::Nsample];
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
		if(!(detectors[j]->get_is_ref())) nref_nb++;
	}
	if(nref_nb>1){
		cout << "there is more than 1 non ref det, exiting" << endl;
		return;
	}

	for(int i=0;i<det_N;i++){
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
			CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(this,this,false,-1);
			vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);

			for(vector<Event*>::iterator it = (currentCBEvent->events).begin();it!=(currentCBEvent->events).end();++it){
				if(!((*it)->get_is_ref())){
					if((*it)->get_type() == Tomography::CM_Demux){
						continue;
					}
					else if((*it)->get_type() == Tomography::MG){
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
								kt->set_perp_pos_mm(*jt);
								double current_residu = jt->get_residu_ref(&(*kt));
								if(current_residu<residu){
									residu = current_residu;
									matching_cluster = kt;
								}
							}
							if(matching_cluster == current_clusters.end()) continue;
							
							double normalization = matching_cluster->get_ampl()/matching_cluster->get_size();
							double matching_position = (matching_cluster->get_is_X()) ? jt->eval_X((*it)->get_z()) : jt->eval_Y((*it)->get_z());
							double matching_position_perp = (matching_cluster->get_is_X()) ? jt->eval_Y((*it)->get_z()) : jt->eval_X((*it)->get_z());
							//double matching_position = matching_cluster->get_pos_mm();
							/*
							for(int strip_nb = matching_cluster->get_pos()-1;strip_nb<(matching_cluster->get_pos()+2);strip_nb++){
								int channel = MG_Detector::StripToChannel[strip_nb];
								double current_max_ampl = *max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::Nsample);
								if(current_max_ampl>normalization) normalization = current_max_ampl;
							}
							*/
							for(int strip_nb=0;strip_nb<1024;strip_nb++){
								int channel = MG_Detector::StripToChannel[strip_nb];
								if(Abs(residu)<50.){
									SRH->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::Nsample))/normalization);
									if(bin_nb>0) SRH_coord[Min(Max(FloorNint((bin_nb/2.)+matching_position_perp*bin_nb/Tomography::XY_size),0),bin_nb-1)]->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::Nsample))/normalization);
								}
								//SRH2D[i]->SetPoint(graph_point_nb, matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::Nsample))/normalization);
								//graph_point_nb++;
								SRH2D->Fill(matching_position - matching_cluster->correct_strip_nb(strip_nb), (*max_element(StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel],StripAmpl_MG_corr[(*it)->get_n_in_tree()][channel]+Tomography::Nsample))/normalization);
							}
							current_clusters.erase(matching_cluster);
						}
					}
				}
			}

			delete currentCBEvent;

			if(jentry%500 == 0) cout << "\r" << setw(20) << i << "|" << setw(20) << jentry << flush;
			if(jentry%10000 == 0 && Tomography::live_graphic_display){
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
			gauss_width_graph->SetPoint(i_coord,(Tomography::XY_size/(bin_nb*2.))+(i_coord*Tomography::XY_size/bin_nb),SRF_coord[i_coord]->GetParameter(1));
			gauss_width_graph->SetPointError(i_coord,(Tomography::XY_size/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(1));
			lorentz_width_graph->SetPoint(i_coord,(Tomography::XY_size/(bin_nb*2.))+(i_coord*Tomography::XY_size/bin_nb),SRF_coord[i_coord]->GetParameter(2));
			lorentz_width_graph->SetPointError(i_coord,(Tomography::XY_size/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(2));
			offset_graph->SetPoint(i_coord,(Tomography::XY_size/(bin_nb*2.))+(i_coord*Tomography::XY_size/bin_nb),SRF_coord[i_coord]->GetParameter(3));
			offset_graph->SetPointError(i_coord,(Tomography::XY_size/(bin_nb*2.)),SRF_coord[i_coord]->GetParError(3));
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
	TTree * signal_tree = (TTree*)(signal_file->Get("T"));
	Tsignal * signalT = new Tsignal(signal_tree,CM_N,MG_N);

	if(nentries != signal_tree->GetEntriesFast()){
		cout << "total number of event in signal and analyse tree does not match" << endl;
		return;
	}
	LoadTree(event_nb);
	GetEntry(event_nb);
	CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,false,-1);
	signalT->LoadTree(event_nb);
	signalT->GetEntry(event_nb);
	for(vector<Event*>::iterator ev_it = (CBEvent->events).begin();ev_it!=(CBEvent->events).end();++ev_it){
		if((*ev_it)->get_type() == Tomography::MG) (*ev_it)->set_strip_ampl(signalT->get_mg_ampl((*ev_it)->get_n_in_tree()));
		else if((*ev_it)->get_type() == Tomography::CM || (*ev_it)->get_type() == Tomography::CM_Demux) (*ev_it)->set_strip_ampl(signalT->get_cm_ampl((*ev_it)->get_n_in_tree()));
	}
	CBEvent->EventDisplay(c1);
	delete CBEvent;
	delete signal_file;
	delete signalT;
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

	float StripAmpl_MG_corr[MG_N][MG_Detector::Nstrip][Tomography::Nsample];
	float StripAmpl_CM_corr[CM_N][CM_Detector::Nstrip][Tomography::Nsample];
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
					double current_ampl = *max_element(StripAmpl_MG_corr[det_n_in_tree][channel],StripAmpl_MG_corr[det_n_in_tree][channel]+Tomography::Nsample);
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
	TH2D * correlation_X_t = new TH2D("correlation_X_t","correlation_X_t",100,0,Tomography::Nsample,100,0,Tomography::Nsample);
	TH2D * correlation_Y_ampl = new TH2D("correlation_Y_ampl","correlation_Y_ampl",1000,0,50000,1000,0,50000);
	TH2D * correlation_Y_t = new TH2D("correlation_Y_t","correlation_Y_t",100,0,Tomography::Nsample,100,0,Tomography::Nsample);
	TH1D * sigma_X_ampl = new TH1D("sigma_X_ampl","sigma_X_ampl",100,0,5000);
	TH1D * sigma_X_t = new TH1D("sigma_X_t","sigma_X_t",100,0,Tomography::Nsample);
	TH1D * sigma_Y_ampl = new TH1D("sigma_Y_ampl","sigma_Y_ampl",100,0,50000);
	TH1D * sigma_Y_t = new TH1D("sigma_Y_t","sigma_Y_t",100,0,Tomography::Nsample);
	TH1D * corr_X_ampl = new TH1D("corr_X_ampl","corr_X_ampl",100,-1,1);
	TH1D * corr_X_t = new TH1D("corr_X_t","corr_X_t",100,-1,1);
	TH1D * corr_Y_ampl = new TH1D("corr_Y_ampl","corr_Y_ampl",100,-1,1);
	TH1D * corr_Y_t = new TH1D("corr_Y_t","corr_Y_t",100,-1,1);
	TCanvas * cCorrXY = new TCanvas();
	cCorrXY->Divide(2,2);
	TH2D * correlation_XY_ampl = new TH2D("correlation_XY_ampl","correlation_XY_ampl",1000,0,50000,1000,0,5000);
	TH2D * correlation_XY_t = new TH2D("correlation_XY_t","correlation_XY_t",100,0,Tomography::Nsample,100,0,Tomography::Nsample);
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
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		CosmicBenchEvent * CBEvent = new CosmicBenchEvent(this,this,false,-1);
		CBEvent->do_cuts();
		eventSuitable+=CBEvent->get_clus_N()/(CM_N+MG_N);
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
				int n = -1;
				if((*det_it)->get_type() == Tomography::MG) n = dynamic_cast<MG_Detector*>(*det_it)->get_mg_n_in_tree();
				else if((*det_it)->get_type() == Tomography::CM) n = dynamic_cast<CM_Detector*>(*det_it)->get_cm_n_in_tree();
				else continue;
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
		if(jentry%5000 == 0 && Tomography::live_graphic_display){
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
	if(CM_N!=0){
		cout << "not implemented with CM" << endl;
		return;
	}
	map<int,TCanvas*> cDisplay;
	map<int,TH1D*> global_signal;
	map<int,TH1D*> global_noise;
	map<int,TProfile*> global_signal_over_noise;
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
			for(int i=0;i<MG_Detector::Nstrip;i++){
				global_noise[current_det->get_mg_n_in_tree()]->Fill(current_det->get_RMS(i));
			}
			global_signal_over_noise[current_det->get_mg_n_in_tree()] = new TProfile((name.str() + "SoB").c_str(),(name.str() + "SoB").c_str(),MG_Detector::Nstrip,0,MG_Detector::Nstrip);
		}
	}
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	for(long i=0;i<nentries;i++){
		LoadTree(i);
		GetEntry(i);

		for(vector<Detector*>::iterator it = detectors.begin();it!=detectors.end();++it){
			if((*it)->get_type() == Tomography::MG){
				MG_Detector * current_det = dynamic_cast<MG_Detector*>(*it);
				MG_Event current_event(this,current_det,false,evn);
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
		TLine * average_SoN = new TLine(0,(res_signal->Parameter(1))/(res_noise->Parameter(1)),MG_Detector::Nstrip,(res_signal->Parameter(1))/(res_noise->Parameter(1)));
		average_SoN->SetLineStyle(2);
		average_SoN->SetLineColor(2);
		TLine * mean_SoN = new TLine(0,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()),MG_Detector::Nstrip,(global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()));
		mean_SoN->SetLineStyle(2);
		mean_SoN->SetLineColor(4);
		global_signal_over_noise[it->first]->Draw();
		average_SoN->Draw();
		mean_SoN->Draw();
		it->second->Modified();
		it->second->Update();
		cout << "MG" << it->first << endl;
		cout << "    mean S/B : " << global_signal_over_noise[it->first]->GetMean(2) << endl;
		cout << "    mean S/mean B : " << (global_signal[it->first]->GetMean())/(global_noise[it->first]->GetMean()) << endl;
		cout << "    sigma S/B : " << global_signal_over_noise[it->first]->GetRMS(2) << endl;
		cout << "    delta mean S/B : " << global_signal_over_noise[it->first]->GetMean(11) << endl;
	}
}