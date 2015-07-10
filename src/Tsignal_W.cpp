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


Tsignal_W::Tsignal_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_)
{
   saveFile = new TFile(saveFileName.c_str(),"RECREATE");
   T = new TTree("T","event");
   T->SetMaxTreeSize(10000000000000LL);
   det_N = det_N_;
   Init();
}

Tsignal_W::~Tsignal_W()
{
   if(det_N[Tomography::MG]>0){
      delete StripAmpl_MG;
      delete StripAmpl_MG_ped;
      delete StripAmpl_MG_corr;
   }
   if(det_N[Tomography::MGv2]>0){
      delete StripAmpl_MGv2;
      delete StripAmpl_MGv2_ped;
      delete StripAmpl_MGv2_corr;
   }
   if(det_N[Tomography::CM]>0){
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
   if(det_N[Tomography::MG]>0){
      StripAmpl_MG = new Float_t[det_N[Tomography::MG]][MG_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_MG_ped = new Float_t[det_N[Tomography::MG]][MG_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_MG_corr = new Float_t[det_N[Tomography::MG]][MG_Detector::Nchannel][Tomography::Nsample];

      char leefStripAmpl_MG[100];
      sprintf(leefStripAmpl_MG,"StripAmpl_MG[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MG", StripAmpl_MG, leefStripAmpl_MG);

      char leefStripAmpl_MG_ped[100];
      sprintf(leefStripAmpl_MG_ped,"StripAmpl_MG_ped[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MG_ped", StripAmpl_MG_ped, leefStripAmpl_MG_ped);

      char leefStripAmpl_MG_corr[100];
      sprintf(leefStripAmpl_MG_corr,"StripAmpl_MG_corr[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MG_corr", StripAmpl_MG_corr, leefStripAmpl_MG_corr);
   }
   if(det_N[Tomography::MGv2]>0){
      StripAmpl_MGv2 = new Float_t[det_N[Tomography::MGv2]][MGv2_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_MGv2_ped = new Float_t[det_N[Tomography::MGv2]][MGv2_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_MGv2_corr = new Float_t[det_N[Tomography::MGv2]][MGv2_Detector::Nchannel][Tomography::Nsample];

      char leefStripAmpl_MGv2[100];
      sprintf(leefStripAmpl_MGv2,"StripAmpl_MGv2[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MGv2", StripAmpl_MGv2, leefStripAmpl_MGv2);

      char leefStripAmpl_MGv2_ped[100];
      sprintf(leefStripAmpl_MGv2_ped,"StripAmpl_MGv2_ped[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MGv2_ped", StripAmpl_MGv2_ped, leefStripAmpl_MGv2_ped);

      char leefStripAmpl_MGv2_corr[100];
      sprintf(leefStripAmpl_MGv2_corr,"StripAmpl_MGv2_corr[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_MGv2_corr", StripAmpl_MGv2_corr, leefStripAmpl_MGv2_corr);
   }
   if(det_N[Tomography::CM]>0){
      StripAmpl_CM = new Float_t[det_N[Tomography::CM]][CM_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_CM_ped = new Float_t[det_N[Tomography::CM]][CM_Detector::Nchannel][Tomography::Nsample];
      StripAmpl_CM_corr = new Float_t[det_N[Tomography::CM]][CM_Detector::Nchannel][Tomography::Nsample];

      char leefStripAmpl_CM[100];
      sprintf(leefStripAmpl_CM,"StripAmpl_CM[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_CM", StripAmpl_CM, leefStripAmpl_CM);

      char leefStripAmpl_CM_ped[100];
      sprintf(leefStripAmpl_CM_ped,"StripAmpl_CM_ped[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::Nsample);
      T->Branch("StripAmpl_CM_ped", StripAmpl_CM_ped, leefStripAmpl_CM_ped);

      char leefStripAmpl_CM_corr[100];
      sprintf(leefStripAmpl_CM_corr,"StripAmpl_CM_corr[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::Nsample);
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
void Tsignal_W::fillTree_raw(int evn_, double evttime_, map<Tomography::det_type,vector<vector<vector<float> > > > ampl){
   Nevent = evn_;
   T->GetBranch("Nevent")->Fill();
   evttime = evttime_;
   T->GetBranch("evttime")->Fill();
   for(map<Tomography::det_type,vector<vector<vector<float> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MG[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MGv2[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_CM[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
   }
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG")->Fill();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2")->Fill();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM")->Fill();
   T->SetEntries((T->GetEntries())+1);
}
void Tsignal_W::fillTree_ped(map<Tomography::det_type,vector<vector<vector<float> > > > ampl){
   for(map<Tomography::det_type,vector<vector<vector<float> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MG_ped[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MGv2_ped[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_CM_ped[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
   }
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG_ped")->Fill();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2_ped")->Fill();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM_ped")->Fill();
}
void Tsignal_W::fillTree_corr(map<Tomography::det_type,vector<vector<vector<float> > > > ampl){
   for(map<Tomography::det_type,vector<vector<vector<float> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MG_corr[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_MGv2_corr[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::Nsample)){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::Nsample;k++){
                  StripAmpl_CM_corr[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
   }
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG_corr")->Fill();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2_corr")->Fill();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM_corr")->Fill();
}
void Tsignal_W::Reset_raw(){
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG")->Reset();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2")->Reset();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM")->Reset();
}
void Tsignal_W::Reset_ped(){
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG_ped")->Reset();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2_ped")->Reset();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM_ped")->Reset();
}
void Tsignal_W::Reset_corr(){
   if(det_N[Tomography::MG]>0) T->GetBranch("StripAmpl_MG_corr")->Reset();
   if(det_N[Tomography::MGv2]>0) T->GetBranch("StripAmpl_MGv2_corr")->Reset();
   if(det_N[Tomography::CM]>0) T->GetBranch("StripAmpl_CM_corr")->Reset();
}
map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_raw(long entry){
   map<Tomography::det_type,vector<vector<vector<float> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   if(det_N[Tomography::MG]>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(det_N[Tomography::MG],vector<vector<float> >(MG_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MG];i++){
      for(int j=0;j<MG_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG[i][j][k];
         }
      }
   }
   if(det_N[Tomography::MGv2]>0) return_map[Tomography::MGv2] = vector<vector<vector<float> > >(det_N[Tomography::MGv2],vector<vector<float> >(MGv2_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MGv2];i++){
      for(int j=0;j<MGv2_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MGv2][i][j][k] = StripAmpl_MGv2[i][j][k];
         }
      }
   }
   if(det_N[Tomography::CM]>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(det_N[Tomography::CM],vector<vector<float> >(CM_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::CM];i++){
      for(int j=0;j<CM_Detector::Nchannel;j++){
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
   if(det_N[Tomography::MG]>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(det_N[Tomography::MG],vector<vector<float> >(MG_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MG];i++){
      for(int j=0;j<MG_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG_ped[i][j][k];
         }
      }
   }
   if(det_N[Tomography::MGv2]>0) return_map[Tomography::MGv2] = vector<vector<vector<float> > >(det_N[Tomography::MGv2],vector<vector<float> >(MGv2_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MGv2];i++){
      for(int j=0;j<MGv2_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MGv2][i][j][k] = StripAmpl_MGv2_ped[i][j][k];
         }
      }
   }
   if(det_N[Tomography::CM]>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(det_N[Tomography::CM],vector<vector<float> >(CM_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::CM];i++){
      for(int j=0;j<CM_Detector::Nchannel;j++){
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
   if(det_N[Tomography::MG]>0) return_map[Tomography::MG] = vector<vector<vector<float> > >(det_N[Tomography::MG],vector<vector<float> >(MG_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MG];i++){
      for(int j=0;j<MG_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MG][i][j][k] = StripAmpl_MG_corr[i][j][k];
         }
      }
   }
   if(det_N[Tomography::MGv2]>0) return_map[Tomography::MGv2] = vector<vector<vector<float> > >(det_N[Tomography::MGv2],vector<vector<float> >(MGv2_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::MGv2];i++){
      for(int j=0;j<MGv2_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::MGv2][i][j][k] = StripAmpl_MGv2_corr[i][j][k];
         }
      }
   }
   if(det_N[Tomography::CM]>0) return_map[Tomography::CM] = vector<vector<vector<float> > >(det_N[Tomography::CM],vector<vector<float> >(CM_Detector::Nchannel,vector<float>(Tomography::Nsample,0)));
   for(int i=0;i<det_N[Tomography::CM];i++){
      for(int j=0;j<CM_Detector::Nchannel;j++){
         for(int k=0;k<Tomography::Nsample;k++){
            return_map[Tomography::CM][i][j][k] = StripAmpl_CM_corr[i][j][k];
         }
      }
   }
   return return_map;
}