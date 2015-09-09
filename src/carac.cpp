#define carac_cpp
#include "carac.h"
#include <iostream>
//#include <string>
#include <iomanip>
#include <ostream>
#include <map>
#include <vector>
#include <utility>
#include <limits>

#include "detector.h"
#include "ray.h"
#include "event.h"
#include "cluster.h"

#include <TFile.h>
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
	Init();
}
Carac::Carac(ptree config_tree_){
	config_tree = config_tree_;
	Init();
}
void Carac::Init(){
	f = new TFile((config_tree.get<string>("Tree")).c_str());
	cout << config_tree.get<string>("Tree") << endl;
	TTree * tree = (TTree*)(f->Get("T"));
	max_event = config_tree.get<long>("max_event");
	Tanalyse_R::Init(tree,CosmicBench(config_tree).get_det_N());
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
	map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> > perp_pairs;
	long nentries = (max_event>0) ? Min(static_cast<long>(fChain->GetEntriesFast()),max_event) : fChain->GetEntriesFast();
	CosmicBench temp_CB(config_tree);
	for(int i=0;i<(temp_CB.get_det_N_tot());i++){
		Detector * current_det = temp_CB.get_detector(i);
		ostringstream name;
		name << current_det->get_type() << "_" << current_det->get_n_in_tree();
		if(current_det->get_perp_n()>-1) perp_pairs[pair<Tomography::det_type,int>(current_det->get_type(),current_det->get_n_in_tree())] = pair<Tomography::det_type,int>(current_det->get_perp_type(),current_det->get_perp_n());
		else{
			perp_pairs[pair<Tomography::det_type,int>(current_det->get_type(),current_det->get_n_in_tree())] = pair<Tomography::det_type,int>(current_det->get_type(),-1);
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
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator map_it = perp_pairs.begin();map_it!=perp_pairs.end();++map_it){
		if(perp_pairs.count(map_it->second)==0 && (map_it->second).second > -1){
			cout << "2D detectors must be set to non ref in both direction" << endl;
			return;
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
	map<string,CosmicBench*> all_CB;
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
		ostringstream name_compl;
		name_compl << (pair_it->first).first << "_" << (pair_it->first).second;
		if((pair_it->second).second > -1) name_compl << "_" << (pair_it->second).first << "_" << (pair_it->second).second;
		chisquares[name_compl.str()] = new TH1D(("chiSquares_" + name_compl.str()).c_str(),("chiSquares_" + name_compl.str()).c_str(),nbins,0,chisquare_threshold);
		chisquares[name_compl.str()]->SetLineColor(color_n);
		ray_clus_n[name_compl.str()] = new TH1D(("clus_n_" + name_compl.str()).c_str(),("clus_n_" + name_compl.str()).c_str(),temp_CB.get_det_N_tot() + 1,0,temp_CB.get_det_N_tot() + 1);
		ray_clus_n[name_compl.str()]->SetLineColor(color_n);
		ray_slope[name_compl.str()] = new TH1D(("slope_" + name_compl.str()).c_str(),("slope_" + name_compl.str()).c_str(),100,0,1);
		ray_slope[name_compl.str()]->SetLineColor(color_n);
		ray_phi[name_compl.str()] = new TH1D(("phi_" + name_compl.str()).c_str(),("phi_" + name_compl.str()).c_str(),100,-Pi(),Pi());
		ray_phi[name_compl.str()]->SetLineColor(color_n);
		color_n++;
		for(map<const Tomography::det_type,const Detector* const>::const_iterator type_it = Tomography::Static_Detector.begin();type_it!=Tomography::Static_Detector.end();++type_it){
			ostringstream tree_child;
			tree_child << "CosmicBench." << type_it->first;
			if(config_tree.get_child_optional(tree_child.str().c_str())){
				for(ptree::iterator tree_it = config_tree.get_child(tree_child.str().c_str()).begin(); tree_it != config_tree.get_child(tree_child.str().c_str()).end(); tree_it++){
					Detector * current_det = (type_it->second)->build_det(*tree_it);
					if(type_it->first == (pair_it->first).first && current_det->get_n_in_tree() == (pair_it->first).second){
						tree_it->second.put<bool>("is_ref",false);
					}
					else if((pair_it->second).second > -1 && type_it->first == (pair_it->second).first && current_det->get_n_in_tree() == (pair_it->second).second){
						tree_it->second.put<bool>("is_ref",false);
					}
					else tree_it->second.put<bool>("is_ref",true);
					delete current_det;
				}
			}
		}
		all_CB[name_compl.str()] = new CosmicBench(config_tree);
	}
	if (fChain == 0) return;
	cout << setw(20) << "total processed" << endl;
	for (Long64_t jentry=0; jentry<nentries && Tomography::can_continue;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		for(map<string,CosmicBench*>::iterator CB_it=all_CB.begin();CB_it!=all_CB.end();++CB_it){
			CosmicBenchEvent * currentCBEvent = new CosmicBenchEvent(CB_it->second,this,-1);
			vector<Ray> currentRays = currentCBEvent->get_absorption_rays(chisquare_threshold);
			vector<Ray>::iterator ray_it = currentRays.begin();
			while(ray_it!= currentRays.end()){
				if(ray_it->get_chiSquare_X()>-1 && ray_it->get_chiSquare_Y()>-1 && ((ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y())/ray_it->get_clus_n())<chisquare_threshold){
					chisquares[CB_it->first]->Fill(ray_it->get_chiSquare_X()+ray_it->get_chiSquare_Y());
					ray_clus_n[CB_it->first]->Fill(ray_it->get_clus_n());
					double slope = Sqrt((ray_it->get_slope_Y()*ray_it->get_slope_Y()) + (ray_it->get_slope_X()*ray_it->get_slope_X()));
					ray_slope[CB_it->first]->Fill(ATan(slope));
					double phi = 2*ATan((ray_it->get_slope_Y())/(slope + ray_it->get_slope_X()));
					if(ray_it->get_slope_X()==0 && ray_it->get_slope_Y()<0) phi = Pi();
					ray_phi[CB_it->first]->Fill(phi);
					++ray_it;
				}
				else ray_it = currentRays.erase(ray_it);
			}
			
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
						if((jt->get_chiSquare_X() + jt->get_chiSquare_Y()) > chisquare_threshold/static_cast<double>(CB_it->second->get_non_ref_N())) continue;
						if(jt->get_clus_n()<static_cast<unsigned int>(CB_it->second->get_det_N_tot()-CB_it->second->get_non_ref_N())) continue;
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
						delete *matching_cluster;
						current_clusters.erase(matching_cluster);
						if(residu<chisquare_threshold){
							muon_seen[name.str()]->Fill(jt->eval_X((*it)->get_z()),jt->eval_Y((*it)->get_z()));
						}
					}
					for(vector<Cluster*>::iterator kt = current_clusters.begin();kt!=current_clusters.end();++kt){
						delete *kt;
					}
				}
			}
			delete currentCBEvent;
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
			for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
				ostringstream name_compl;
				name_compl << (pair_it->first).first << "_" << (pair_it->first).second;
				if((pair_it->second).second > -1) name_compl << "_" << (pair_it->second).first << "_" << (pair_it->second).second;
				string option = (pair_it == perp_pairs.begin()) ? "" : "SAME";
				c0->cd(1);
				ray_clus_n[name_compl.str()]->Draw(option.c_str());
				c0->cd(2);
				chisquares[name_compl.str()]->Draw(option.c_str());
				c0->cd(3);
				ray_slope[name_compl.str()]->Draw(option.c_str());
				c0->cd(4);
				ray_phi[name_compl.str()]->Draw(option.c_str());
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
	for(map<pair<Tomography::det_type,int>,pair<Tomography::det_type,int> >::iterator pair_it = perp_pairs.begin();pair_it!=perp_pairs.end();++pair_it){
		ostringstream name_compl;
		name_compl << (pair_it->first).first << "_" << (pair_it->first).second;
		if((pair_it->second).second > -1) name_compl << "_" << (pair_it->second).first << "_" << (pair_it->second).second;
		string option = (pair_it == perp_pairs.begin()) ? "" : "SAME";
		c0->cd(1);
		ray_clus_n[name_compl.str()]->Draw(option.c_str());
		c0->cd(2);
		chisquares[name_compl.str()]->Draw(option.c_str());
		c0->cd(3);
		ray_slope[name_compl.str()]->Draw(option.c_str());
		c0->cd(4);
		ray_phi[name_compl.str()]->Draw(option.c_str());
		delete all_CB[name_compl.str()];
	}
	c0->Modified();
	c0->Update();
}