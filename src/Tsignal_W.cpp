#define Tsignal_W_cxx
#include "Tsignal_W.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>

using std::vector;

#include <iostream>
using std::cout;
using std::endl;


Tsignal_W::Tsignal_W(string saveFileName, int CM_n, int MG_n)
{
   saveFile = new TFile(saveFileName.c_str(),"RECREATE");
   T = new TTree("T","event");
   T->SetMaxTreeSize(10000000000000LL);
   CMN = CM_n;
   MGN = MG_n;
   Init();
}

Tsignal_W::~Tsignal_W()
{
   if(MGN>0){
      delete StripAmpl_MG;
      delete StripAmpl_MG_ped;
      delete StripAmpl_MG_corr;
   }
   if(CMN>0){
      delete StripAmpl_CM;
      delete StripAmpl_CM_ped;
      delete StripAmpl_CM_corr;
   }
   delete saveFile;
}

void Tsignal_W::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   T->Branch("Nevent", &Nevent, "Nevent/I");
   T->Branch("evttime", &evttime, "evttime/D");
   if(MGN>0){
      StripAmpl_MG = new Float_t[MGN][61][Tomography::Nsample];
      StripAmpl_MG_ped = new Float_t[MGN][61][Tomography::Nsample];
      StripAmpl_MG_corr = new Float_t[MGN][61][Tomography::Nsample];

      char leefStripAmpl_MG[100];
      sprintf(leefStripAmpl_MG,"StripAmpl_MG[%d][%d][%d]/F",MGN,61,Tomography::Nsample);
      T->Branch("StripAmpl_MG", StripAmpl_MG, leefStripAmpl_MG);

      char leefStripAmpl_MG_ped[100];
      sprintf(leefStripAmpl_MG_ped,"StripAmpl_MG_ped[%d][%d][%d]/F",MGN,61,Tomography::Nsample);
      T->Branch("StripAmpl_MG_ped", StripAmpl_MG_ped, leefStripAmpl_MG_ped);

      char leefStripAmpl_MG_corr[100];
      sprintf(leefStripAmpl_MG_corr,"StripAmpl_MG_corr[%d][%d][%d]/F",MGN,61,Tomography::Nsample);
      T->Branch("StripAmpl_MG_corr", StripAmpl_MG_corr, leefStripAmpl_MG_corr);
   }
   if(CMN>0){
      StripAmpl_CM = new Float_t[CMN][64][Tomography::Nsample];
      StripAmpl_CM_ped = new Float_t[CMN][64][Tomography::Nsample];
      StripAmpl_CM_corr = new Float_t[CMN][64][Tomography::Nsample];

      char leefStripAmpl_CM[100];
      sprintf(leefStripAmpl_CM,"StripAmpl_CM[%d][%d][%d]/F",CMN,64,Tomography::Nsample);
      T->Branch("StripAmpl_CM", StripAmpl_CM, leefStripAmpl_CM);

      char leefStripAmpl_CM_ped[100];
      sprintf(leefStripAmpl_CM_ped,"StripAmpl_CM_ped[%d][%d][%d]/F",CMN,64,Tomography::Nsample);
      T->Branch("StripAmpl_CM_ped", StripAmpl_CM_ped, leefStripAmpl_CM_ped);

      char leefStripAmpl_CM_corr[100];
      sprintf(leefStripAmpl_CM_corr,"StripAmpl_CM_corr[%d][%d][%d]/F",CMN,64,Tomography::Nsample);
      T->Branch("StripAmpl_CM_corr", StripAmpl_CM_corr, leefStripAmpl_CM_corr);
   }

   // Set branch addresses and branch pointers
   /*
   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("evttime",&evttime, &b_evttime);
   fChain->SetBranchAddress("TsampleNum", TsampleNum, &b_TsampleNum);
   if(MGN>0){
	   fChain->SetBranchAddress("StripAmpl_MG", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_MG_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_MG_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
	if(CMN>0){
	   fChain->SetBranchAddress("StripAmpl_CM", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_CM_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_CM_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
   */
}

TTree * Tsignal_W::getTree() const{
   return T->CloneTree();
}
void Tsignal_W::Write(){
   saveFile->cd();
   T->Write();
}
void Tsignal_W::CloseFile(){
   saveFile->Close();
}
void Tsignal_W::fillTree_raw(int evn_, double evttime_, vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl){
   Nevent = evn_;
   T->GetBranch("Nevent")->Fill();
   evttime = evttime_;
   T->GetBranch("evttime")->Fill();
   if(mg_ampl.size()!=static_cast<unsigned int>(MGN)){
      cout << "problem in MG number" << endl;
      return;
   }
   for(int i=0;i<MGN;i++){
      if(mg_ampl[i].size()!=61){
         cout << "problem in MG strip number" << endl;
         return;
      }
      for(int j=0;j<61;j++){
         if(mg_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_MG[i][j][k] = mg_ampl[i][j][k];
         }
      }
   }
   if(MGN>0) T->GetBranch("StripAmpl_MG")->Fill();
   if(cm_ampl.size()!=static_cast<unsigned int>(CMN)){
      cout << "problem in CM number" << endl;
      return;
   }
   for(int i=0;i<CMN;i++){
      if(cm_ampl[i].size()!=64){
         cout << "problem in CM strip number" << endl;
         return;
      }
      for(int j=0;j<64;j++){
         if(cm_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_CM[i][j][k] = cm_ampl[i][j][k];
         }
      }
   }
   if(CMN>0) T->GetBranch("StripAmpl_CM")->Fill();
   T->SetEntries((T->GetEntries())+1);
}
void Tsignal_W::fillTree_ped(vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl){
   if(mg_ampl.size()!=static_cast<unsigned int>(MGN)){
      cout << "problem in MG number" << endl;
      return;
   }
   for(int i=0;i<MGN;i++){
      if(mg_ampl[i].size()!=61){
         cout << "problem in MG strip number" << endl;
         return;
      }
      for(int j=0;j<61;j++){
         if(mg_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_MG_ped[i][j][k] = mg_ampl[i][j][k];
         }
      }
   }
   if(MGN>0) T->GetBranch("StripAmpl_MG_ped")->Fill();
   if(cm_ampl.size()!=static_cast<unsigned int>(CMN)){
      cout << "problem in CM number" << endl;
      return;
   }
   for(int i=0;i<CMN;i++){
      if(cm_ampl[i].size()!=64){
         cout << "problem in CM strip number" << endl;
         return;
      }
      for(int j=0;j<64;j++){
         if(cm_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_CM_ped[i][j][k] = cm_ampl[i][j][k];
         }
      }
   }
   if(CMN>0) T->GetBranch("StripAmpl_CM_ped")->Fill();
}
void Tsignal_W::fillTree_corr(vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl){
   if(mg_ampl.size()!=static_cast<unsigned int>(MGN)){
      cout << "problem in MG number" << endl;
      return;
   }
   for(int i=0;i<MGN;i++){
      if(mg_ampl[i].size()!=61){
         cout << "problem in MG strip number" << endl;
         return;
      }
      for(int j=0;j<61;j++){
         if(mg_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_MG_corr[i][j][k] = mg_ampl[i][j][k];
         }
      }
   }
   if(MGN>0) T->GetBranch("StripAmpl_MG_corr")->Fill();
   if(cm_ampl.size()!=static_cast<unsigned int>(CMN)){
      cout << "problem in CM number" << endl;
      return;
   }
   for(int i=0;i<CMN;i++){
      if(cm_ampl[i].size()!=64){
         cout << "problem in CM strip number" << endl;
         return;
      }
      for(int j=0;j<64;j++){
         if(cm_ampl[i][j].size()!=static_cast<unsigned int>(Tomography::Nsample)){
            cout << "problem in sample number" << endl;
            return;
         }
         for(int k=0;k<Tomography::Nsample;k++){
            StripAmpl_CM_corr[i][j][k] = cm_ampl[i][j][k];
         }
      }
   }
   if(CMN>0) T->GetBranch("StripAmpl_CM_corr")->Fill();
}
void Tsignal_W::Reset_raw(){
   if(MGN>0) T->GetBranch("StripAmpl_MG")->Reset();
   if(CMN>0) T->GetBranch("StripAmpl_CM")->Reset();
}
void Tsignal_W::Reset_ped(){
   if(MGN>0) T->GetBranch("StripAmpl_MG_ped")->Reset();
   if(CMN>0) T->GetBranch("StripAmpl_CM_ped")->Reset();
}
void Tsignal_W::Reset_corr(){
   if(MGN>0) T->GetBranch("StripAmpl_MG_corr")->Reset();
   if(CMN>0) T->GetBranch("StripAmpl_CM_corr")->Reset();
}
map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_raw(long entry){
   map<Tomography::det_type,vector<vector<vector<float> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   if(MGN>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(MGN,vector<vector<float> >(61,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<MGN;i++){
      for(int j=0;j<61;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG[i][j][k];
         }
      }
   }
   if(CMN>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(CMN,vector<vector<float> >(64,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<CMN;i++){
      for(int j=0;j<64;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::CM][i][j][k] = StripAmpl_CM[i][j][k];
         }
      }
   }
   return return_map;
}
map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_ped(long entry){
   map<Tomography::det_type,vector<vector<vector<float> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   if(MGN>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(MGN,vector<vector<float> >(61,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<MGN;i++){
      for(int j=0;j<61;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG_ped[i][j][k];
         }
      }
   }
   if(CMN>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(CMN,vector<vector<float> >(64,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<CMN;i++){
      for(int j=0;j<64;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::CM][i][j][k] = StripAmpl_CM_ped[i][j][k];
         }
      }
   }
   return return_map;
}
map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_corr(long entry){
   map<Tomography::det_type,vector<vector<vector<float> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   if(MGN>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(MGN,vector<vector<float> >(61,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<MGN;i++){
      for(int j=0;j<61;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG_corr[i][j][k];
         }
      }
   }
   if(CMN>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(CMN,vector<vector<float> >(64,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<CMN;i++){
      for(int j=0;j<64;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::CM][i][j][k] = StripAmpl_CM_corr[i][j][k];
         }
      }
   }
   return return_map;
}