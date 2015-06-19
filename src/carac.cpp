#define carac_cpp
#include "carac.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <ostream>
#include <map>
#include <vector>
#include <utility>
#include <limits>

#include "detector.h"
#include "ray.h"
#include "event.h"

#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TMath.h>

#include <boost/property_tree/json_parser.hpp>

using std::cout;
using std::endl;
using std::flush;
using std::ostringstream;
using std::map;
using std::vector;
using std::pair;
using std::setw;
using std::sort;
using std::numeric_limits;

using TMath::Min;
using TMath::Sqrt;
using TMath::Pi;
using TMath::Abs;
using TMath::ATan;

Carac::Carac(string configFilePath){
	read_json(configFilePath, config_tree);
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CM_N = config_tree.get<int>("total_CM_N");
	MG_N = config_tree.get<int>("total_MG_N");
	Tanalyse_R::Init(tree,CM_N,MG_N);
}
Carac::Carac(ptree config_tree_){
	config_tree = config_tree_;
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	CM_N = config_tree.get<int>("total_CM_N");
	MG_N = config_tree.get<int>("total_MG_N");
	Tanalyse_R::Init(tree,CM_N,MG_N);
}
Carac::~Carac(){
	//delete f;
}

void Carac::Residus_ref(){
	double chisquare_threshold = 10;
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
	vector< pair<int,int> > perp_pairs;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	CosmicBench temp_CB(config_tree);
	for(int i=0;i<(temp_CB.get_CM_N() + temp_CB.get_MG_N());i++){
		Detector * current_det = temp_CB.get_detector(i);
		ostringstream name;
		if(current_det->get_type() == Tomography::CM){
			name << "Cosmulti_" << (dynamic_cast<CM_Detector*>(current_det))->get_cm_n_in_tree();
		}
		else if(current_det->get_type() == Tomography::MG){
			name << "Multigen_" << (dynamic_cast<MG_Detector*>(current_det))->get_mg_n_in_tree();
			if(current_det->get_perp_n()>-1) perp_pairs.push_back(pair<int,int>((dynamic_cast<MG_Detector*>(current_det))->get_mg_n_in_tree(),current_det->get_perp_n()));
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
	}
	sort(perp_pairs.begin(),perp_pairs.end());
	for(vector<pair<int,int> >::iterator det_it = perp_pairs.begin();det_it!=perp_pairs.end();++det_it){
		vector<pair<int,int> >::iterator det_jt = (det_it+1);
		while(det_jt!=perp_pairs.end()){
			if(det_jt->first == det_it->first){
				if(det_jt->second != det_it->second){
					cout << "error in 2D detector pairing" << endl;
					return;
				}
				else{
					det_jt = perp_pairs.erase(det_jt);
				}
			}
			else if(det_jt->second == det_it->first){
				if(det_jt->first != det_it->second){
					cout << "error in 2D detector pairing" << endl;
					return;
				}
				else{
					det_jt = perp_pairs.erase(det_jt);
				}
			}
			else{
				det_jt++;
			}
		}
	}
	TCanvas * c0 = new TCanvas("stats","stats");
	c0->Divide(2,2);
	c0->GetPad(1)->SetLogy();
	c0->GetPad(2)->SetLogy();
	map<string,TH1D*> chisquares;
	map<string,TH1D*> ray_clus_n;
	map<string,TH1D*> ray_slope;
	map<string,TH1D*> ray_phi;
	int color_n = 1;
	for(vector<pair<int,int> >::iterator det_it = perp_pairs.begin();det_it!=perp_pairs.end();++det_it){
		ostringstream name_compl;
		name_compl << "_" << det_it->first << "_" << det_it->second;
		chisquares["MG"+name_compl.str()] = new TH1D(("chiSquares" + name_compl.str()).c_str(),("chiSquares" + name_compl.str()).c_str(),nbins,0,chisquare_threshold);
		chisquares["MG"+name_compl.str()]->SetLineColor(color_n);
		ray_clus_n["MG"+name_compl.str()] = new TH1D(("clus_n" + name_compl.str()).c_str(),("clus_n" + name_compl.str()).c_str(),MG_N + CM_N + 1,0,MG_N + CM_N + 1);
		ray_clus_n["MG"+name_compl.str()]->SetLineColor(color_n);
		ray_slope["MG"+name_compl.str()] = new TH1D(("slope" + name_compl.str()).c_str(),("slope" + name_compl.str()).c_str(),100,0,1);
		ray_slope["MG"+name_compl.str()]->SetLineColor(color_n);
		ray_phi["MG"+name_compl.str()] = new TH1D(("phi" + name_compl.str()).c_str(),("phi" + name_compl.str()).c_str(),100,-Pi(),Pi());
		ray_phi["MG"+name_compl.str()]->SetLineColor(color_n);
		color_n++;
	}
	if (fChain == 0) return;
	cout << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::can_continue;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		for(vector<pair<int,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
			ostringstream name_compl;
			name_compl << "_" << pair_it->first << "_" << pair_it->second;
			int non_ref_n = 0;
			for(ptree::iterator tree_it = config_tree.get_child("CosmicBench.MultiGens").begin(); tree_it != config_tree.get_child("CosmicBench.MultiGens").end(); tree_it++){
				if(tree_it->second.get<int>("mg_n") == pair_it->first || tree_it->second.get<int>("mg_n") == pair_it->second){
					tree_it->second.put<bool>("is_ref",false);
					non_ref_n++;
				}
				else{
					tree_it->second.put<bool>("is_ref",true);
				}
			}
			CosmicBench * current_CB = new CosmicBench(config_tree);
			CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(current_CB,this,false,-1);
			vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
			vector<Ray>::iterator ray_it = currentRays.begin();
			while(ray_it!= currentRays.end()){
				if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
					chisquares["MG"+name_compl.str()]->Fill(ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y());
					ray_clus_n["MG"+name_compl.str()]->Fill(ray_it->get_clus_n());
					double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
					ray_slope["MG"+name_compl.str()]->Fill(ATan(slope));
					double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
					if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
					ray_phi["MG"+name_compl.str()]->Fill(phi);
					++ray_it;
				}
				else ray_it = currentRays.erase(ray_it);
			}
			
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
			delete current_CB; delete currentCBEvent;
		}
		if(jentry%500 == 0) cout << "\r" << setw(20) << jentry << flush;
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
			for(vector<pair<int,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
				ostringstream name_compl;
				name_compl << "_" << pair_it->first << "_" << pair_it->second;
				string option = (pair_it == perp_pairs.begin()) ? "" : "SAME";
				c0->cd(1);
				ray_clus_n["MG"+name_compl.str()]->Draw(option.c_str());
				c0->cd(2);
				chisquares["MG"+name_compl.str()]->Draw(option.c_str());
				c0->cd(3);
				ray_slope["MG"+name_compl.str()]->Draw(option.c_str());
				c0->cd(4);
				ray_phi["MG"+name_compl.str()]->Draw(option.c_str());
			}
			c0->Modified();
			c0->Update();
		}
	}
	cout << "\r" << setw(20) << nentries << endl;
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
		it->second->Modified();
		it->second->Update();
		cout << it->first << " efficacity : " << 100.*efficacity[it->first] << "%" << endl;
	}
	for(vector<pair<int,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
		ostringstream name_compl;
		name_compl << "_" << pair_it->first << "_" << pair_it->second;
		string option = (pair_it == perp_pairs.begin()) ? "" : "SAME";
		c0->cd(1);
		ray_clus_n["MG"+name_compl.str()]->Draw(option.c_str());
		c0->cd(2);
		chisquares["MG"+name_compl.str()]->Draw(option.c_str());
		c0->cd(3);
		ray_slope["MG"+name_compl.str()]->Draw(option.c_str());
		c0->cd(4);
		ray_phi["MG"+name_compl.str()]->Draw(option.c_str());
	}
	c0->Modified();
	c0->Update();
}