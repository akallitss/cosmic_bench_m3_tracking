#define Tsignal_W_cxx
#include "Tsignal_W.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::ostringstream;


Tsignal_W::Tsignal_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_)
{
   saveFile = new TFile(saveFileName.c_str(),"UPDATE");
   saveFile->SetCompressionLevel(1);
   T = (TTree*)saveFile->Get("T");
   if(!T){
      saveFile->SetCompressionLevel(5);
      T = new TTree("T","event");
   }
   T->SetMaxTreeSize(10000000000000LL);
   det_N = det_N_;
   Init();
}

Tsignal_W::~Tsignal_W()
{
   /*
   if(det_N.count(Tomography::MG)>0){
      delete[] StripAmpl_MG;
      delete[] StripAmpl_MG_ped;
      delete[] StripAmpl_MG_corr;
   }
   if(det_N.count(Tomography::MGv2)>0){
      delete[] StripAmpl_MGv2;
      delete[] StripAmpl_MGv2_ped;
      delete[] StripAmpl_MGv2_corr;
   }
   if(det_N.count(Tomography::CM)>0){
      delete[] StripAmpl_CM;
      delete[] StripAmpl_CM_ped;
      delete[] StripAmpl_CM_corr;
   }
   */
   for(map<Tomography::det_type,Float_t*>::iterator type_it=StripAmpl.begin();type_it!=StripAmpl.end();++type_it){
      delete[] type_it->second;
   }
   StripAmpl.clear();
   for(map<Tomography::det_type,Float_t*>::iterator type_it=StripAmpl_ped.begin();type_it!=StripAmpl_ped.end();++type_it){
      delete[] type_it->second;
   }
   StripAmpl_ped.clear();
   for(map<Tomography::det_type,Float_t*>::iterator type_it=StripAmpl_corr.begin();type_it!=StripAmpl_corr.end();++type_it){
      delete[] type_it->second;
   }
   StripAmpl_corr.clear();
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
   if(T->GetBranch("Nevent")) T->SetBranchAddress("Nevent",&Nevent);
   else T->Branch("Nevent", &Nevent, "Nevent/I");
   if(T->GetBranch("evttime")) T->SetBranchAddress("evttime",&evttime);
   else T->Branch("evttime", &evttime, "evttime/D");
   /*
   if(det_N.count(Tomography::MG)>0){
      StripAmpl_MG = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MG_ped = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MG_corr = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];

      char leefStripAmpl_MG[100];
      sprintf(leefStripAmpl_MG,"StripAmpl_MG[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MG")) T->SetBranchAddress("StripAmpl_MG",StripAmpl_MG);
      else T->Branch("StripAmpl_MG", StripAmpl_MG, leefStripAmpl_MG);

      char leefStripAmpl_MG_ped[100];
      sprintf(leefStripAmpl_MG_ped,"StripAmpl_MG_ped[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MG_ped")) T->SetBranchAddress("StripAmpl_MG_ped",StripAmpl_MG_ped);
      else T->Branch("StripAmpl_MG_ped", StripAmpl_MG_ped, leefStripAmpl_MG_ped);

      char leefStripAmpl_MG_corr[100];
      sprintf(leefStripAmpl_MG_corr,"StripAmpl_MG_corr[%d][%d][%d]/F",det_N[Tomography::MG],MG_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MG_corr")) T->SetBranchAddress("StripAmpl_MG_corr",StripAmpl_MG_corr);
      else T->Branch("StripAmpl_MG_corr", StripAmpl_MG_corr, leefStripAmpl_MG_corr);
   }
   if(det_N.count(Tomography::MGv2)>0){
      StripAmpl_MGv2 = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MGv2_ped = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MGv2_corr = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];

      char leefStripAmpl_MGv2[100];
      sprintf(leefStripAmpl_MGv2,"StripAmpl_MGv2[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MGv2")) T->SetBranchAddress("StripAmpl_MGv2",StripAmpl_MGv2);
      else T->Branch("StripAmpl_MGv2", StripAmpl_MGv2, leefStripAmpl_MGv2);

      char leefStripAmpl_MGv2_ped[100];
      sprintf(leefStripAmpl_MGv2_ped,"StripAmpl_MGv2_ped[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MGv2_ped")) T->SetBranchAddress("StripAmpl_MGv2_ped",StripAmpl_MGv2_ped);
      else T->Branch("StripAmpl_MGv2_ped", StripAmpl_MGv2_ped, leefStripAmpl_MGv2_ped);

      char leefStripAmpl_MGv2_corr[100];
      sprintf(leefStripAmpl_MGv2_corr,"StripAmpl_MGv2_corr[%d][%d][%d]/F",det_N[Tomography::MGv2],MGv2_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_MGv2_corr")) T->SetBranchAddress("StripAmpl_MGv2_corr",StripAmpl_MGv2_corr);
      else T->Branch("StripAmpl_MGv2_corr", StripAmpl_MGv2_corr, leefStripAmpl_MGv2_corr);
   }
   if(det_N.count(Tomography::CM)>0){
      StripAmpl_CM = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_CM_ped = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_CM_corr = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];

      char leefStripAmpl_CM[100];
      sprintf(leefStripAmpl_CM,"StripAmpl_CM[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_CM")) T->SetBranchAddress("StripAmpl_CM",StripAmpl_CM);
      else T->Branch("StripAmpl_CM", StripAmpl_CM, leefStripAmpl_CM);

      char leefStripAmpl_CM_ped[100];
      sprintf(leefStripAmpl_CM_ped,"StripAmpl_CM_ped[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_CM_ped")) T->SetBranchAddress("StripAmpl_CM_ped",StripAmpl_CM_ped);
      else T->Branch("StripAmpl_CM_ped", StripAmpl_CM_ped, leefStripAmpl_CM_ped);

      char leefStripAmpl_CM_corr[100];
      sprintf(leefStripAmpl_CM_corr,"StripAmpl_CM_corr[%d][%d][%d]/F",det_N[Tomography::CM],CM_Detector::Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch("StripAmpl_CM_corr")) T->SetBranchAddress("StripAmpl_CM_corr",StripAmpl_CM_corr);
      else T->Branch("StripAmpl_CM_corr", StripAmpl_CM_corr, leefStripAmpl_CM_corr);
   }
   */
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      int current_Nchannel = Tomography::Static_Detector[type_it->first]->get_Nchannel();
      StripAmpl[type_it->first] = new Float_t[(type_it->second)*(Tomography::get_instance()->get_Nsample())*current_Nchannel];
      StripAmpl_ped[type_it->first] = new Float_t[(type_it->second)*(Tomography::get_instance()->get_Nsample())*current_Nchannel];
      StripAmpl_corr[type_it->first] = new Float_t[(type_it->second)*(Tomography::get_instance()->get_Nsample())*current_Nchannel];
      string current_name = Tomography::Static_Detector[type_it->first]->Name();
      char leefStripAmpl[100];
      sprintf(leefStripAmpl,"StripAmpl_%s[%d][%d][%d]/F",current_name.c_str(),type_it->second,current_Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch(("StripAmpl_"+current_name).c_str())) T->SetBranchAddress(("StripAmpl_"+current_name).c_str(),StripAmpl[type_it->first]);
      else T->Branch(("StripAmpl_"+current_name).c_str(),StripAmpl[type_it->first],leefStripAmpl);

      char leefStripAmpl_ped[100];
      sprintf(leefStripAmpl_ped,"StripAmpl_%s_ped[%d][%d][%d]/F",current_name.c_str(),type_it->second,current_Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch(("StripAmpl_"+current_name+"_ped").c_str())) T->SetBranchAddress(("StripAmpl_"+current_name+"_ped").c_str(),StripAmpl_ped[type_it->first]);
      else T->Branch(("StripAmpl_"+current_name+"_ped").c_str(),StripAmpl_ped[type_it->first],leefStripAmpl_ped);

      char leefStripAmpl_corr[100];
      sprintf(leefStripAmpl_corr,"StripAmpl_%s_corr[%d][%d][%d]/F",current_name.c_str(),type_it->second,current_Nchannel,Tomography::get_instance()->get_Nsample());
      if(T->GetBranch(("StripAmpl_"+current_name+"_corr").c_str())) T->SetBranchAddress(("StripAmpl_"+current_name+"_corr").c_str(),StripAmpl_corr[type_it->first]);
      else T->Branch(("StripAmpl_"+current_name+"_corr").c_str(),StripAmpl_corr[type_it->first],leefStripAmpl_corr);
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
template<typename A>
void Tsignal_W::fillTree_raw(int evn_, double evttime_, map<Tomography::det_type,vector<vector<vector<A> > > > ampl){
   Nevent = evn_;
   T->GetBranch("Nevent")->Fill();
   evttime = evttime_;
   T->GetBranch("evttime")->Fill();
   for(typename map<Tomography::det_type,vector<vector<vector<A> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if(det_N.count(type_it->first)==0){
         cout << "detector type mismatch : " << type_it->first << endl;
         return;
      }
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      for(int i=0;i<det_N[type_it->first];i++){
         if((type_it->second)[i].size() != static_cast<unsigned int>(Tomography::Static_Detector[type_it->first]->get_Nchannel())){
            cout << "problem in " << type_it->first << " channel/strip number" << endl;
            return;
         }
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
               cout << "problem in sample number" << endl;
               return;
            }
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl[type_it->first])[i][j][k] = (type_it->second)[i][j][k];
            }
         }
      }
      T->GetBranch(("StripAmpl_"+Tomography::Static_Detector[type_it->first]->Name()).c_str())->Fill();
      /*
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else{
         cout << "det type not implemented : " << type_it->first << endl;
	      return;
      }
      */
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG")->Fill();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2")->Fill();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM")->Fill();
   */
   T->SetEntries((T->GetEntries())+1);
}
template void Tsignal_W::fillTree_raw(int evn_, double evttime_, map<Tomography::det_type,vector<vector<vector<float> > > > ampl);
template void Tsignal_W::fillTree_raw(int evn_, double evttime_, map<Tomography::det_type,vector<vector<vector<double> > > > ampl);

template<typename A>
void Tsignal_W::fillTree_ped(map<Tomography::det_type,vector<vector<vector<A> > > > ampl){
   for(typename map<Tomography::det_type,vector<vector<vector<A> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if(det_N.count(type_it->first)==0){
         cout << "detector type mismatch : " << type_it->first << endl;
         return;
      }
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      for(int i=0;i<det_N[type_it->first];i++){
         if((type_it->second)[i].size() != static_cast<unsigned int>(Tomography::Static_Detector[type_it->first]->get_Nchannel())){
            cout << "problem in " << type_it->first << " channel/strip number" << endl;
            return;
         }
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
               cout << "problem in sample number" << endl;
               return;
            }
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl_ped[type_it->first])[i][j][k] = (type_it->second)[i][j][k];
            }
         }
      }
      T->GetBranch(("StripAmpl_"+Tomography::Static_Detector[type_it->first]->Name()+"_ped").c_str())->Fill();
      /*
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_ped)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_ped)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_ped)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else{
         cout << "det type not implemented : " << type_it->first << endl;
	      return;
      }
      */
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG_ped")->Fill();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2_ped")->Fill();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM_ped")->Fill();
   */
}
template void Tsignal_W::fillTree_ped(map<Tomography::det_type,vector<vector<vector<float> > > > ampl);
template void Tsignal_W::fillTree_ped(map<Tomography::det_type,vector<vector<vector<double> > > > ampl);

template<typename A>
void Tsignal_W::fillTree_corr(map<Tomography::det_type,vector<vector<vector<A> > > > ampl){
   for(typename map<Tomography::det_type,vector<vector<vector<A> > > >::const_iterator type_it=ampl.begin();type_it!=ampl.end();++type_it){
      if(det_N.count(type_it->first)==0){
         cout << "detector type mismatch : " << type_it->first << endl;
         return;
      }
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in " << type_it->first << " number" << endl;
         return;
      }
      for(int i=0;i<det_N[type_it->first];i++){
         if((type_it->second)[i].size() != static_cast<unsigned int>(Tomography::Static_Detector[type_it->first]->get_Nchannel())){
            cout << "problem in " << type_it->first << " channel/strip number" << endl;
            return;
         }
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
               cout << "problem in sample number" << endl;
               return;
            }
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl_corr[type_it->first])[i][j][k] = (type_it->second)[i][j][k];
            }
         }
      }
      T->GetBranch(("StripAmpl_"+Tomography::Static_Detector[type_it->first]->Name()+"_corr").c_str())->Fill();
      /*
      if(type_it->first == Tomography::MG){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MG_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MG_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_corr)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::MGv2){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != MGv2_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<MGv2_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_corr)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else if(type_it->first == Tomography::CM){
         for(int i=0;i<det_N[type_it->first];i++){
            if((type_it->second)[i].size() != CM_Detector::Nchannel){
               cout << "problem in " << type_it->first << " strip number" << endl;
               return;
            }
            for(int j=0;j<CM_Detector::Nchannel;j++){
               if((type_it->second)[i][j].size() != static_cast<unsigned int>(Tomography::get_instance()->get_Nsample())){
                  cout << "problem in sample number" << endl;
                  return;
               }
               for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
                  reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_corr)[i][j][k] = (type_it->second)[i][j][k];
               }
            }
         }
      }
      else{
          cout << "det type not implemented : " << type_it->first << endl;
	       return;
      }
      */
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG_corr")->Fill();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2_corr")->Fill();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM_corr")->Fill();
   */
}
template void Tsignal_W::fillTree_corr(map<Tomography::det_type,vector<vector<vector<float> > > > ampl);
template void Tsignal_W::fillTree_corr(map<Tomography::det_type,vector<vector<vector<double> > > > ampl);

void Tsignal_W::Reset_raw(){
   for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      if(type_it->second > 0){
         ostringstream name;
         name << "StripAmpl_" << type_it->first;
         T->GetBranch(name.str().c_str())->Reset();
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG")->Reset();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2")->Reset();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM")->Reset();
   */
}
void Tsignal_W::Reset_ped(){
   for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      if(type_it->second > 0){
         ostringstream name;
         name << "StripAmpl_" << type_it->first << "_ped";
         T->GetBranch(name.str().c_str())->Reset();
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG_ped")->Reset();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2_ped")->Reset();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM_ped")->Reset();
   */
}
void Tsignal_W::Reset_corr(){
   for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      if(type_it->second > 0){
         ostringstream name;
         name << "StripAmpl_" << type_it->first << "_corr";
         T->GetBranch(name.str().c_str())->Reset();
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0) T->GetBranch("StripAmpl_MG_corr")->Reset();
   if(det_N.count(Tomography::MGv2)>0) T->GetBranch("StripAmpl_MGv2_corr")->Reset();
   if(det_N.count(Tomography::CM)>0) T->GetBranch("StripAmpl_CM_corr")->Reset();
   */
}
void Tsignal_W::disable_data_branches(){
   T->SetBranchStatus("*",0);
   T->SetBranchStatus("Nevent",1);
   T->SetBranchStatus("evttime",1);
}
void Tsignal_W::enable_all_branches(){
   T->SetBranchStatus("*",1);
}
void Tsignal_W::enable_raw_branches(){
   for(map<Tomography::det_type,unsigned short>::const_iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      if(type_it->second > 0){
         ostringstream name;
         name << "StripAmpl_" << type_it->first;
         T->SetBranchStatus(name.str().c_str(),1);
      }
   }
}
void Tsignal_W::enable_ped_branches(){
   T->SetBranchStatus("StripAmpl_*_ped",1);
}
void Tsignal_W::enable_corr_branches(){
   T->SetBranchStatus("StripAmpl_*_corr",1);
}
template<typename A>
map<Tomography::det_type,vector<vector<vector<A> > > > Tsignal_W::read_raw(long entry){
   map<Tomography::det_type,vector<vector<vector<A> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      return_map[type_it->first] = vector<vector<vector<A> > >(det_N[type_it->first],vector<vector<A> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<type_it->second;i++){
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[type_it->first][i][j][k] = reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl[type_it->first])[i][j][k];
            }
         }
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0){
      return_map[Tomography::MG] = vector<vector<vector<A> > >(det_N[Tomography::MG],vector<vector<A> >(MG_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MG];i++){
         for(int j=0;j<MG_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MG][i][j][k] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::MGv2)>0){
      return_map[Tomography::MGv2] = vector<vector<vector<A> > >(det_N[Tomography::MGv2],vector<vector<A> >(MGv2_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MGv2];i++){
         for(int j=0;j<MGv2_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MGv2][i][j][k] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::CM)>0){
      return_map[Tomography::CM] = vector<vector<vector<A> > >(det_N[Tomography::CM],vector<vector<A> >(CM_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::CM];i++){
         for(int j=0;j<CM_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::CM][i][j][k] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM)[i][j][k];
            }
         }
      }
   }
   */
   return return_map;
}
template map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_raw(long entry);
template map<Tomography::det_type,vector<vector<vector<double> > > > Tsignal_W::read_raw(long entry);

template<typename A>
map<Tomography::det_type,vector<vector<vector<A> > > > Tsignal_W::read_ped(long entry){
   map<Tomography::det_type,vector<vector<vector<A> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      return_map[type_it->first] = vector<vector<vector<A> > >(det_N[type_it->first],vector<vector<A> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<type_it->second;i++){
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[type_it->first][i][j][k] = reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl_ped[type_it->first])[i][j][k];
            }
         }
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0){
      return_map[Tomography::MG] = vector<vector<vector<A> > >(det_N[Tomography::MG],vector<vector<A> >(MG_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MG];i++){
         for(int j=0;j<MG_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MG][i][j][k] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_ped)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::MGv2)>0){
      return_map[Tomography::MGv2] = vector<vector<vector<A> > >(det_N[Tomography::MGv2],vector<vector<A> >(MGv2_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MGv2];i++){
         for(int j=0;j<MGv2_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MGv2][i][j][k] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_ped)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::CM)>0){
      return_map[Tomography::CM] = vector<vector<vector<A> > >(det_N[Tomography::CM],vector<vector<A> >(CM_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::CM];i++){
         for(int j=0;j<CM_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::CM][i][j][k] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_ped)[i][j][k];
            }
         }
      }
   }
   */
   return return_map;
}
template map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_ped(long entry);
template map<Tomography::det_type,vector<vector<vector<double> > > > Tsignal_W::read_ped(long entry);

template<typename A>
map<Tomography::det_type,vector<vector<vector<A> > > > Tsignal_W::read_corr(long entry){
   map<Tomography::det_type,vector<vector<vector<A> > > > return_map;
   if(entry>= T->GetEntries()){
      return return_map;
   }
   T->LoadTree(entry);
   T->GetEntry(entry);
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      return_map[type_it->first] = vector<vector<vector<A> > >(det_N[type_it->first],vector<vector<A> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<type_it->second;i++){
         for(int j=0;j<Tomography::Static_Detector[type_it->first]->get_Nchannel();j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[type_it->first][i][j][k] = reinterpret_cast<Float_t(*)[Tomography::Static_Detector[type_it->first]->get_Nchannel()][Tomography::get_instance()->get_Nsample()]>(StripAmpl_corr[type_it->first])[i][j][k];
            }
         }
      }
   }
   /*
   if(det_N.count(Tomography::MG)>0){
      return_map[Tomography::MG] = vector<vector<vector<A> > >(det_N[Tomography::MG],vector<vector<A> >(MG_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MG];i++){
         for(int j=0;j<MG_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MG][i][j][k] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_corr)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::MGv2)>0){
      return_map[Tomography::MGv2] = vector<vector<vector<A> > >(det_N[Tomography::MGv2],vector<vector<A> >(MGv2_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::MGv2];i++){
         for(int j=0;j<MGv2_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::MGv2][i][j][k] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_corr)[i][j][k];
            }
         }
      }
   }
   if(det_N.count(Tomography::CM)>0){
      return_map[Tomography::CM] = vector<vector<vector<A> > >(det_N[Tomography::CM],vector<vector<A> >(CM_Detector::Nchannel,vector<A>(Tomography::get_instance()->get_Nsample(),0)));
      for(int i=0;i<det_N[Tomography::CM];i++){
         for(int j=0;j<CM_Detector::Nchannel;j++){
            for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
               return_map[Tomography::CM][i][j][k] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_corr)[i][j][k];
            }
         }
      }
   }
   */
   return return_map;
}
template map<Tomography::det_type,vector<vector<vector<float> > > > Tsignal_W::read_corr(long entry);
template map<Tomography::det_type,vector<vector<vector<double> > > > Tsignal_W::read_corr(long entry);
