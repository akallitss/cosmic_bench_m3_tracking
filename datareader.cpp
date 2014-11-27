#define datareader_cpp
#include "datareader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>
#include <TBranch.h>
#include <TH1F.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <TMath.h>
#include <sstream>
#include "dataline.h"
#include "header.h"

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::sort;
using std::setw;
using std::setfill;
using std::ostringstream;

using TMath::Min;

const int DataReader::Nsample = 32;
const int DataReader::Nstrip_MG = 61;
const int DataReader::Nstrip_CM = 64;

DataReader::DataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_,bool ped_done_,bool cns_done_, int max_event_){
	exists = exists_;
	ped_done = ped_done_;
	cns_done = cns_done_;
	max_event = max_event_;
	global_offset = 0;
	MG_N = 0;
	CM_N = 0;
	det_type_by_asic = det_type_by_asic_;
	det_n_by_asic = det_n_by_asic_;
	if(det_n_by_asic.size()!= det_type_by_asic.size()){
		cout << "problem in detector caracs" << endl;
		return;
	}
	for(unsigned int i=0;i<det_n_by_asic.size();i++){
		if(det_n_by_asic[i]>=0){
			if(det_type_by_asic[i] == "MG") MG_N++;
			else if(det_type_by_asic[i] == "CM") CM_N++;
			else{
				cout << "detector type unknown : " << det_type_by_asic[i] << endl;
				return;
			}
		}
	}
	if(MG_N>0){
		StripAmpl_MG = new float[MG_N][Nstrip_MG][Nsample];
		StripAmpl_MG_ped = new float[MG_N][Nstrip_MG][Nsample];
		StripAmpl_MG_corr = new float[MG_N][Nstrip_MG][Nsample];
		Pedestal_MG = new float[MG_N][Nstrip_MG];
	}
	if(CM_N>0){
		StripAmpl_CM = new float[CM_N][Nstrip_CM][Nsample];
		StripAmpl_CM_ped = new float[CM_N][Nstrip_CM][Nsample];
		StripAmpl_CM_corr = new float[CM_N][Nstrip_CM][Nsample];
		Pedestal_CM = new float[CM_N][Nstrip_CM];
	}
	outFileName = baseFileName + "_signal.root";
	PedFileName = baseFileName + "_Ped.dat";
	RMSPedFileName = baseFileName + "_RMSPed.dat";
	string fileOption = (exists) ? "UPDATE" : "RECREATE";
	outFile = new TFile(outFileName.c_str(),fileOption.c_str());
	outTree->SetMaxTreeSize(100000000000LL);
	Nevent = 0;
	if(!exists){
		outTree = new TTree("T","event");
		outTree->Branch("Nevent", &Nevent, "Nevent/I"); // event number
		char leefTsampleNum[100];
		sprintf(leefTsampleNum,"TsampleNum[%d]/I",Nsample);
		outTree->Branch("TsampleNum", TsampleNum, leefTsampleNum); // time sample number
		// For the MG
		if(MG_N>0){
			char leefStripAmpl_MG[100];
			sprintf(leefStripAmpl_MG,"StripAmpl_MG[%d][%d][%d]/F",MG_N,Nstrip_MG,Nsample);
			outTree->Branch("StripAmpl_MG", StripAmpl_MG, leefStripAmpl_MG); // raw amplitude
		}
		// For the CM
		if(CM_N>0){
			char leefStripAmpl_CM[100];
			sprintf(leefStripAmpl_CM,"StripAmpl_CM[%d][%d][%d]/F",CM_N,Nstrip_CM,Nsample);
			outTree->Branch("StripAmpl_CM", StripAmpl_CM, leefStripAmpl_CM); // raw amplitude
		}
	}
	else{
		outTree = (TTree*)(outFile->Get("T"));
		outTree->SetBranchAddress("Nevent",&Nevent);
		outTree->SetBranchAddress("TsampleNum",TsampleNum);
		if(CM_N>0) outTree->SetBranchAddress("StripAmpl_CM",StripAmpl_CM);
		if(MG_N>0) outTree->SetBranchAddress("StripAmpl_MG",StripAmpl_MG);
	}
	if(ped_done){
		if(CM_N>0) outTree->SetBranchAddress("StripAmpl_CM_ped",StripAmpl_CM_ped);
		if(MG_N>0) outTree->SetBranchAddress("StripAmpl_MG_ped",StripAmpl_MG_ped);
	}
	if(cns_done){
		if(CM_N>0) outTree->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
		if(MG_N>0) outTree->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);
	}
	for(int k=0;k<Nsample;k++){
		TsampleNum[k] = k;
	}
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			for(int k=0;k<Nsample;k++){
				StripAmpl_MG[i][j][k] = 0;
			}
		}
	}
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			for(int k=0;k<Nsample;k++){
				StripAmpl_CM[i][j][k] = 0;
			}
		}
	}
	is_first = true;
	DAQType = "";
	dumb_branch = new TBranch();
}
DataReader::~DataReader(){
	delete dumb_branch;
	outFile->Close();
	delete outFile;
}
void DataReader::add_file_to_process(string inFileName){
	file_names.push_back(inFileName);
}
void DataReader::Fill(){
	outTree->Fill();
}
void DataReader::Write(){
	outTree->Write();
}
void DataReader::reset_tree_leaf(){
	//Nevent = 0;
	/*
	for(int k=0;k<Nsample;k++){
		TsampleNum[k] = k;
	}
	*/
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			for(int k=0;k<Nsample;k++){
				StripAmpl_MG[i][j][k] = 0;
			}
		}
	}
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			for(int k=0;k<Nsample;k++){
				StripAmpl_CM[i][j][k] = 0;
			}
		}
	}
}
void DataReader::compute_ped(){
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			Pedestal_MG[i][j] = 0;
		}
	}
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			Pedestal_CM[i][j] = 0;
		}
	}
	int nentries = outTree->GetEntries();
	for(int n=0;n<nentries;n++){
		outTree->LoadTree(n);
		outTree->GetEntry(n);
		for(unsigned int i=0;i<MG_N;i++){
			for(int j=0;j<Nstrip_MG;j++){
				vector<float> current_strip(Nsample,0);
				for(int k=0;k<Nsample;k++){
					current_strip[k] = StripAmpl_MG[i][j][k];
				}
				sort(current_strip.begin(),current_strip.end());
				Pedestal_MG[i][j] += current_strip[Nsample/2];
			}
		}
		for(unsigned int i=0;i<CM_N;i++){
			for(int j=0;j<Nstrip_CM;j++){
				vector<float> current_strip(Nsample,0);
				for(int k=0;k<Nsample;k++){
					current_strip[k] = StripAmpl_CM[i][j][k];
				}
				sort(current_strip.begin(),current_strip.end());
				Pedestal_CM[i][j] += current_strip[Nsample/2];
			}
		}
		if((n%100) == 0) cout << "\r" << "computing pedestal (" << n << "/" << nentries << ")" << flush;
	}
	cout << "\r" << "computing pedestal (" << nentries << "/" << nentries << ")" << endl;
	cout << "writing it to file..." << flush;
	ofstream pedFile(PedFileName.c_str());
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			Pedestal_CM[i][j] /= nentries;
			pedFile << i << " " << j << " " << Pedestal_CM[i][j] << "\n";
		}
	}
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			Pedestal_MG[i][j] /= nentries;
			pedFile << i << " " << j << " " << Pedestal_MG[i][j] << "\n";
		}
	}
	pedFile.close();
	cout << "!" << endl;
}

void DataReader::read_ped(string ped_file){
	if(ped_file != "") PedFileName = ped_file;
	cout << "reading pedestal from file : " << PedFileName << "..." << flush;
	ifstream pedFile(PedFileName.c_str());
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			pedFile >> i >> j >> Pedestal_CM[i][j];
		}
	}
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			pedFile >> i >> j >> Pedestal_MG[i][j];
		}
	}
	pedFile.close();
	cout << "!" << endl;
}
void DataReader::do_ped_sub(string ped_file){
	if(ped_done){
		cout << "pedestal already substracted" << endl;
		return;
	}
	read_ped(ped_file);
	TBranch *newBranch_MG = dumb_branch;
	TBranch *newBranch_CM = dumb_branch;
	if(MG_N>0){
		char leefStripAmpl_MG[100];
		sprintf(leefStripAmpl_MG,"StripAmpl_MG_ped[%d][%d][%d]/F",MG_N,Nstrip_MG,Nsample);
		newBranch_MG = outTree->Branch("StripAmpl_MG_ped", StripAmpl_MG_ped,leefStripAmpl_MG);
	}
	if(CM_N>0){
		char leefStripAmpl_CM[100];
		sprintf(leefStripAmpl_CM,"StripAmpl_CM_ped[%d][%d][%d]/F",CM_N,Nstrip_CM,Nsample);
		newBranch_CM = outTree->Branch("StripAmpl_CM_ped", StripAmpl_CM_ped,leefStripAmpl_CM);
	}
	int nentries = outTree->GetEntries();
	for(int n=0;n<nentries;n++){
		outTree->LoadTree(n);
		outTree->GetEntry(n);
		for(unsigned int i=0;i<MG_N;i++){
			for(int j=0;j<Nstrip_MG;j++){
				for(int k=0;k<Nsample;k++){
					StripAmpl_MG_ped[i][j][k] = StripAmpl_MG[i][j][k] - Pedestal_MG[i][j];
				}
			}
		}
		for(unsigned int i=0;i<CM_N;i++){
			for(int j=0;j<Nstrip_CM;j++){
				for(int k=0;k<Nsample;k++){
					StripAmpl_CM_ped[i][j][k] = StripAmpl_CM[i][j][k] - Pedestal_CM[i][j];
				}
			}
		}
		if(CM_N>0) newBranch_CM->Fill();
		if(MG_N>0) newBranch_MG->Fill();
		if((n%100) == 0) cout << "\r" << "substracting pedestal (" << n << "/" << nentries << ")" << flush;
	}
	cout << "\r" << "substracting pedestal (" << nentries << "/" << nentries << ")" << endl;
	Write();
	ped_done = true;
}
void DataReader::do_common_noise_sub(){
	if(cns_done){
		cout << "common noise already substracted" << endl;
		return;
	}
	if(!ped_done){
		cout << "cannot substract common noise before the pedestal substraction" << endl;
		return;
	}
	TBranch *newBranch_MG = dumb_branch;
	TBranch *newBranch_CM = dumb_branch;
	if(MG_N>0){
		char leefStripAmpl_MG[100];
		sprintf(leefStripAmpl_MG,"StripAmpl_MG_corr[%d][%d][%d]/F",MG_N,Nstrip_MG,Nsample);
		newBranch_MG = outTree->Branch("StripAmpl_MG_corr", StripAmpl_MG_corr,leefStripAmpl_MG);
	}
	if(CM_N>0){
		char leefStripAmpl_CM[100];
		sprintf(leefStripAmpl_CM,"StripAmpl_CM_corr[%d][%d][%d]/F",CM_N,Nstrip_CM,Nsample);
		newBranch_CM = outTree->Branch("StripAmpl_CM_corr", StripAmpl_CM_corr,leefStripAmpl_CM);
	}
	int detector_div_MG = 2;
	int detector_div_CM = 2;
	int nentries = outTree->GetEntries();
	for(int n=0;n<nentries;n++){
		outTree->LoadTree(n);
		outTree->GetEntry(n);
		for(unsigned int i=0;i<MG_N;i++){
			for(int k=0;k<Nsample;k++){
				for(int det_div=0;det_div<detector_div_MG;det_div++){
					int strip_nb = (Nstrip_MG/detector_div_MG) + (Nstrip_MG%detector_div_MG);
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && j+strip_offset<Nstrip_MG);j++){
						current_sample[j] = StripAmpl_MG_ped[i][j+strip_offset][k];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && j+strip_offset<Nstrip_MG);j++){
						StripAmpl_MG_corr[i][j+strip_offset][k] = StripAmpl_MG_ped[i][j+strip_offset][k] - median;
					}
				}
			}
		}
		for(unsigned int i=0;i<CM_N;i++){
			for(int k=0;k<Nsample;k++){
				for(int det_div=0;det_div<detector_div_CM;det_div++){
					int strip_nb = Nstrip_CM/detector_div_CM + (Nstrip_CM%detector_div_CM);
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && j+strip_offset<Nstrip_CM);j++){
						current_sample[j] = StripAmpl_CM_ped[i][j+strip_offset][k];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && j+strip_offset<Nstrip_CM);j++){
						StripAmpl_CM_corr[i][j+strip_offset][k] = StripAmpl_CM_ped[i][j+strip_offset][k] - median;
					}
				}
			}
		}
		if(CM_N>0) newBranch_CM->Fill();
		if(MG_N>0) newBranch_MG->Fill();
		if((n%100) == 0) cout << "\r" << "substracting common noise (" << n << "/" << nentries << ")" << flush;
	}
	cout << "\r" << "substracting common noise (" << nentries << "/" << nentries << ")" << endl;
	Write();
	cns_done = true;
}
void DataReader::compute_RMSPed(){
	double Ymin=-500;
	double Ymax=500;
	int bin_n = 500;
	Long64_t max_event = 500;
	int nentries = Min(outTree->GetEntries(),max_event);
	ofstream RMSPedFile(RMSPedFileName.c_str());
	if(CM_N>0){
		outTree->SetBranchStatus("*",0);
		outTree->SetBranchStatus("StripAmpl_CM_corr",1);
		for(unsigned int i=0;i<CM_N;i++){
			cout << "\r" << "computing RMS Ped for CM_" << i << flush;
			vector<TH1F*> ampl_hist(Nstrip_CM);
			for(int j=0;j<Nstrip_CM;j++){
				ostringstream name;
				name << "ampl_hist_" << j;
				ampl_hist[j] = new TH1F(name.str().c_str(),name.str().c_str(),bin_n,Ymin,Ymax);
			}
			for(int n=0;n<nentries;n++){
				outTree->LoadTree(n);
				outTree->GetEntry(n);
				for(int j=0;j<Nstrip_CM;j++){
					for(int k=0;k<Nsample;k++){
						ampl_hist[j]->Fill(StripAmpl_CM_corr[i][j][k]);
					}
				}
			}
			for(int j=0;j<Nstrip_CM;j++){
				TFitResultPtr res = ampl_hist[j]->Fit("gaus","SQN");
				RMSPedFile << i << " " << j << " " << res->Parameter(2) << "\n";
				delete ampl_hist[j];
			}
			/*
			for(int j=0;j<Nstrip_CM;j++){
				cout << "\r" << "computing RMS Ped for CM_" << i << " and strip_" << setw(2) << setfill('0') << j << flush;
				TH1F * ampl_hist = new TH1F("ampl_hist","ampl_hist",bin_n,Ymin,Ymax);
				for(int n=0;n<nentries;n++){
					outTree->LoadTree(n);
					outTree->GetEntry(n);
					for(int k=0;k<Nsample;k++){
						ampl_hist->Fill(StripAmpl_CM_corr[i][j][k]);
					}
				}
				TFitResultPtr res = ampl_hist->Fit("gaus","SQN");
				RMSPedFile << i << " " << j << " " << res->Parameter(2) << "\n";
				delete ampl_hist;
			}
			*/
		}
		cout << "\r" << "RMS Ped for CMs computed !                             " << endl;
	}
	if(MG_N>0){
		outTree->SetBranchStatus("*",0);
		outTree->SetBranchStatus("StripAmpl_MG_corr",1);
		for(unsigned int i=0;i<MG_N;i++){
			cout << "\r" << "computing RMS Ped for MG_" << i << flush;
			vector<TH1F*> ampl_hist(Nstrip_MG);
			for(int j=0;j<Nstrip_MG;j++){
				ostringstream name;
				name << "ampl_hist_" << j;
				ampl_hist[j] = new TH1F(name.str().c_str(),name.str().c_str(),bin_n,Ymin,Ymax);
			}
			for(int n=0;n<nentries;n++){
				outTree->LoadTree(n);
				outTree->GetEntry(n);
				for(int j=0;j<Nstrip_MG;j++){
					for(int k=0;k<Nsample;k++){
						ampl_hist[j]->Fill(StripAmpl_MG_corr[i][j][k]);
					}
				}
			}
			for(int j=0;j<Nstrip_MG;j++){
				TFitResultPtr res = ampl_hist[j]->Fit("gaus","SQN");
				RMSPedFile << i << " " << j << " " << res->Parameter(2) << "\n";
				delete ampl_hist[j];
			}
			/*
			for(int j=0;j<Nstrip_MG;j++){
				cout << "\r" << "computing RMS Ped for MG_" << i << " and strip_" << setw(2) << setfill('0') << j << flush;
				TH1F * ampl_hist = new TH1F("ampl_hist","ampl_hist",bin_n,Ymin,Ymax);
				for(int n=0;n<nentries;n++){
					outTree->LoadTree(n);
					outTree->GetEntry(n);
					for(int k=0;k<Nsample;k++){
						ampl_hist->Fill(StripAmpl_MG_corr[i][j][k]);
					}
				}
				TFitResultPtr res = ampl_hist->Fit("gaus","SQN");
				RMSPedFile << i << " " << j << " " << res->Parameter(2) << "\n";
				delete ampl_hist;
			}
			*/
		}
		cout << "\r" << "RMS Ped for MGs computed !                             " << endl;
	}
	outTree->SetBranchStatus("*",1);
	RMSPedFile.close();
}
DreamDataReader::DreamDataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_,bool ped_done_,bool cns_done_, int max_event_): DataReader(baseFileName,det_type_by_asic_,det_n_by_asic_,exists_,ped_done_,cns_done_,max_event_){
	DAQType = "Dream";
}
DreamDataReader::~DreamDataReader(){

}
int DreamDataReader::mapping(string det_type, int channel){
	if(det_type == "MG"){
		return channel + 1 - (2*(channel%2));
	}
	return channel;
}
void DreamDataReader::process(){
	if(exists){
		cout << "tree already initiated" << endl;
		return;
	}
	for(vector<string>::iterator it=file_names.begin();it!=file_names.end();++it){
		int current_offset = outTree->GetEntries();
		read_file(*it,current_offset);
	}
	Write();
	exists = true;
}
void DreamDataReader::read_file(string file_name,int evn_offset){
	ifstream iFile(file_name.c_str(),ifstream::binary);
	if(!iFile.is_open()){
		cout << "file : " << file_name << " can't be opened" << endl;
		return;
	}
	HeaderC current_header;
	int evNinFile = 0;
	// Loop on event
	int isample=-1; int isample_prev=-2;
	int ichannel=0;
	int asicN=0;
	int detN=0;
	int channelN=0;
	DataLineDream current_data;
	while(iFile.good() && evNinFile<26000 && !(((evNinFile + evn_offset - global_offset)>max_event)*(max_event>0))){
		reset_tree_leaf();
		isample=-1; isample_prev=-2;
		while(isample<Nsample-1 && isample==isample_prev+1){
			iFile.ignore(2);
			if(iFile.read((char*)&current_header,sizeof(current_header)).good()){
				current_header.ntohs_();
				isample_prev = isample;
				isample = current_header.get_sampleIndex();
				if(isample!=isample_prev+1){
					cout << "problem in sample index" << endl;
					return;
				}
				if(!current_header.check_type()){
					cout << "problem with header type" << endl;
					return;
				}
				if(is_first){
					cout << "INFO: pedestal mode is: " << current_header.get_ped_mode() << endl;
					cout << "INFO: common noise mode is: " << current_header.get_cms_mode() << endl;
					cout << "INFO: zero suppress mode is: " << current_header.get_zs_mode() << endl;
					is_first = false;
				}
				ichannel=0;
				asicN=0;
				detN=0;
				channelN=0;
				current_data.data = 0;
				iFile.read((char*)&current_data,sizeof(current_data));
				current_data.ntohs_();
				while(!(current_data.is_final_trailer())){
					if(current_data.is_first_line()){
						iFile.ignore(2*sizeof(current_data));
						iFile.read((char*)&current_data,sizeof(current_data));
						current_data.ntohs_();
						asicN = current_data.get_dream_ID();
						detN = det_n_by_asic[asicN];
					}
					if(current_data.is_channel_ID() && current_header.get_zs_mode()){
						ichannel = current_data.get_channel_ID();
						channelN = mapping(det_type_by_asic[asicN],ichannel);
						iFile.read((char*)&current_data,sizeof(current_data));
						current_data.ntohs_();

						if(det_type_by_asic[asicN] == "MG"){
							if(channelN>-1 && channelN<Nstrip_MG) StripAmpl_MG[detN][channelN][isample] = current_data.get_data();
							TsampleNum[isample] = isample;
						}
						else if(det_type_by_asic[asicN] == "CM"){
							if(channelN>-1 && channelN<Nstrip_CM) StripAmpl_CM[detN][channelN][isample] = current_data.get_data();
							TsampleNum[isample] = isample;
						}
					}
					if(current_data.is_data() && !(current_header.get_zs_mode())){
						channelN = mapping(det_type_by_asic[asicN],ichannel);
						if(det_type_by_asic[asicN] == "MG"){
							if(channelN>-1 && channelN<Nstrip_MG) StripAmpl_MG[detN][channelN][isample] = current_data.get_data();
							TsampleNum[isample] = isample;
						}
						else if(det_type_by_asic[asicN] == "CM"){
							if(channelN>-1 && channelN<Nstrip_CM) StripAmpl_CM[detN][channelN][isample] = current_data.get_data();
							TsampleNum[isample] = isample;
						}
						ichannel++;
						ichannel = ichannel%64;
					}
					iFile.read((char*)&current_data,sizeof(current_data));
					current_data.ntohs_();
				}
				iFile.read((char*)&current_data,sizeof(current_data));
				current_data.ntohs_();
			}
			else{
				iFile.close();
				break;
			}
		}
		Nevent = evNinFile+evn_offset;
		if((evNinFile%100) == 0) cout << "\r" << "event processed in file : " << file_name << " : " << evNinFile << " (total number of event : " << evNinFile + evn_offset - global_offset << ")" << flush;
		evNinFile++;
		Fill();
	}
	cout << "\r" << "event processed in file : " << file_name << " : " << evNinFile << " (total number of event : " << evNinFile + evn_offset - global_offset << ")" << endl;
	iFile.close();
}
map<string,vector<vector<vector<double> > > > DreamDataReader::read_event(ifstream * file,int event_nb, bool fill_tree){
	map<string,vector<vector<vector<double> > > > event_ampl;
	return event_ampl;
}
FeminosDataReader::FeminosDataReader(string baseFileName, map<int,string> det_type_by_asic_, map<int,int> det_n_by_asic_, bool exists_,bool ped_done_,bool cns_done_, int max_event_): DataReader(baseFileName,det_type_by_asic_,det_n_by_asic_,exists_,ped_done_,cns_done_,max_event_){
	DAQType = "Feminos";
}
FeminosDataReader::~FeminosDataReader(){

}
int FeminosDataReader::mapping(string det_type, int channel){
	if(det_type == "MG"){
		int tmpchan = -1;
		if(channel<14) tmpchan = channel-2;
		else if(channel<25) tmpchan = channel-3;
		else if(channel<48) tmpchan = channel-4;
		else if(channel<59) tmpchan = channel-5;
		else tmpchan = channel-6;
		if(tmpchan<16 || tmpchan>47) tmpchan = 62 + (2*(tmpchan%2)) - tmpchan;
		if(tmpchan>15 && tmpchan<48) return tmpchan;
		else return 63 - tmpchan;
	}
	return channel;
}
void FeminosDataReader::process(){
	if(exists){
		cout << "tree already initiated" << endl;
		return;
	}
	global_offset = get_first_event_nb(file_names.front());
	Nevent = global_offset;
	for(vector<string>::iterator it=file_names.begin();it!=file_names.end();++it){
		//int current_offset = global_offset + outTree->GetEntries();
		read_file(*it,Nevent);
		Nevent++;
	}
	Write();
	exists = true;
}
void FeminosDataReader::read_file(string file_name,int evn_offset){
	ifstream iFile(file_name.c_str(),ifstream::binary);
	if(!iFile.is_open()){
		cout << "file : " << file_name << " can't be opened" << endl;
		return;
	}
	iFile.ignore(26); //You can read run UID here
	iFile.ignore(2); //Beginning of frame (28 bytes from the beginning)
	int evNinFile = 0;
	int card=0;
	int chip=0;
	int channel=0;
	int itime = 0;
	int channelN=0;
	int det = 0;
	int detN=0;
	bool inEvent = false;
	bool inFrame = false;
	int event_started = 0;
	DataLineFeminos current_data;
	iFile.read((char*)&current_data,sizeof(current_data));
	while(iFile.good() && !(((evNinFile + evn_offset - global_offset)>max_event)*(max_event>0))){


		if(inEvent){
			if(inFrame){
				if(current_data.is_event_start()){
					event_started++;
					iFile.ignore(3*sizeof(current_data)); //contain timestamp
					int current_event;
					iFile.read((char*)&current_event,sizeof(current_event));
					if(current_event != (evNinFile + evn_offset)){
						cout << "problem in event number" << endl;
						cout << std::hex << current_event << " " << evNinFile << endl;
						return;
					}
				}
				else if(current_data.is_end_of_event()){
					event_started--;
					iFile.ignore(sizeof(current_data)); //contain eventsize
				}
				else if(current_data.is_end_of_frame()){
					inFrame = false;
				}
				else if(current_data.is_info()){
					card = current_data.get_card_ID();
					chip = current_data.get_chip_ID();
					det = chip + (4*card);
					detN = det_n_by_asic[det];
					channel = current_data.get_channel_ID();
					channelN = mapping(det_type_by_asic[det],channel);
					itime = 0;
				}
				else if(current_data.is_time()){
					itime= current_data.get_time();
				}
				else if(current_data.is_data()){
					if(det_type_by_asic[det] == "MG" && channelN>-1 && channelN<61){
						StripAmpl_MG[detN][channelN][itime] = current_data.get_data();
						itime++;
					}
					else if(det_type_by_asic[det] == "CM" && channelN>-1 && channelN<64){
						StripAmpl_CM[detN][channelN][itime] = current_data.get_data();
						itime++;
					}
				}

			}
			else if(current_data.is_end_of_built_event()){
				if(event_started != 0){
					cout << "problem in fem number" << endl;
					return;
				}
				inEvent = false;
				Nevent = evNinFile + evn_offset;
				if((evNinFile%100) == 0) cout << "\r" << "event processed in file : " << file_name << " : " << evNinFile << " (total number of event : " << Nevent - global_offset << ")" << flush;
				evNinFile++;
				Fill();
			}
			else if(current_data.is_frame_start()){
				inFrame = true;
			}
		}
		else if(current_data.is_built_event_start()){
			inEvent = true;
			reset_tree_leaf();
			chip=0;
			channel=0;
			itime = 0;
			channelN=0;
			det=0;
			detN=0;
			event_started = 0;
		}
		iFile.read((char*)&current_data,sizeof(current_data));
	}
	cout << "\r" << "event processed in file : " << file_name << " : " << evNinFile << " (total number of event : " << evNinFile + evn_offset - global_offset << ")" << endl;
	iFile.close();
}

map<string,vector<vector<vector<double> > > > FeminosDataReader::read_event(ifstream * file,int event_nb, bool fill_tree){
	int card=0;
	int chip=0;
	int channel=0;
	int itime = 0;
	int channelN=0;
	int det = 0;
	int detN=0;
	bool inEvent = false;
	bool inFrame = false;
	int event_started = 0;
	DataLineFeminos current_data;
	file->read((char*)&current_data,sizeof(current_data));
	bool event_complete = false;
	map<string,vector<vector<vector<double> > > > event_ampl;
	while(file->good()){


		if(inEvent){
			if(inFrame){
				if(current_data.is_event_start()){
					event_started++;
					file->ignore(3*sizeof(current_data)); //contain timestamp
					int current_event;
					file->read((char*)&current_event,sizeof(current_event));
					if(current_event != event_nb){
						cout << "warning : event numbers does not match" << endl;
						cout << current_event << " " << event_nb << endl;
					}
				}
				else if(current_data.is_end_of_event()){
					event_started--;
					file->ignore(sizeof(current_data)); //contain eventsize
				}
				else if(current_data.is_end_of_frame()){
					inFrame = false;
				}
				else if(current_data.is_info()){
					card = current_data.get_card_ID();
					chip = current_data.get_chip_ID();
					det = chip + (4*card);
					detN = det_n_by_asic[det];
					channel = current_data.get_channel_ID();
					channelN = mapping(det_type_by_asic[det],channel);
					itime = 0;
				}
				else if(current_data.is_time()){
					itime= current_data.get_time();
				}
				else if(current_data.is_data()){
					if(det_type_by_asic[det] == "MG" && channelN>-1 && channelN<61){
						StripAmpl_MG[detN][channelN][itime] = current_data.get_data();
						itime++;
					}
					else if(det_type_by_asic[det] == "CM" && channelN>-1 && channelN<64){
						StripAmpl_CM[detN][channelN][itime] = current_data.get_data();
						itime++;
					}
				}

			}
			else if(current_data.is_end_of_built_event()){
				if(event_started != 0){
					cout << "problem in fem number" << endl;
					return event_ampl;
				}
				inEvent = false;
				event_complete = true;
				break;
			}
			else if(current_data.is_frame_start()){
				inFrame = true;
			}
		}
		else if(current_data.is_built_event_start()){
			inEvent = true;
			reset_tree_leaf();
			chip=0;
			channel=0;
			itime = 0;
			channelN=0;
			det=0;
			detN=0;
			event_started = 0;
		}
		file->read((char*)&current_data,sizeof(current_data));
	}
	if(!event_complete) return event_ampl;
	Nevent = event_nb;
	if(fill_tree && !exists) Fill();
	vector<vector<vector<double> > > MG_Ampl(MG_N,vector<vector<double> >(Nstrip_MG,vector<double>(Nsample,0)));
	vector<vector<vector<double> > > CM_Ampl(CM_N,vector<vector<double> >(Nstrip_CM,vector<double>(Nsample,0)));
	for(unsigned int i=0;i<MG_N;i++){
		for(int j=0;j<Nstrip_MG;j++){
			for(int k=0;k<Nsample;k++){
				MG_Ampl[i][j][k] = StripAmpl_MG[i][j][k];
			}
		}
	}
	for(unsigned int i=0;i<CM_N;i++){
		for(int j=0;j<Nstrip_CM;j++){
			for(int k=0;k<Nsample;k++){
				CM_Ampl[i][j][k] = StripAmpl_CM[i][j][k];
			}
		}
	}
	event_ampl["MG"] = MG_Ampl;
	event_ampl["CM"] = CM_Ampl;
	return event_ampl;
}

int FeminosDataReader::get_first_event_nb(string file_name){
	ifstream iFile(file_name.c_str(),ifstream::binary);
	if(!iFile.is_open()){
		cout << "file : " << file_name << " can't be opened" << endl;
		return -1;
	}
	iFile.ignore(26); //You can read run UID here
	iFile.ignore(2); //Beginning of frame (28 bytes from the beginning)
	DataLineFeminos current_data;
	while(iFile.good() /*&& evNinFile<6000*/){
		if(current_data.is_event_start()){
			iFile.ignore(3*sizeof(current_data)); //contain timestamp
			int current_event;
			iFile.read((char*)&current_event,sizeof(current_event));
			return current_event;
		}
		iFile.read((char*)&current_data,sizeof(current_data));
	}
	iFile.close();
	return -1;
}