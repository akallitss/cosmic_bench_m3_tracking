#define Tsignal_R_cxx
#include "Tsignal_R.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using std::cout;
using std::endl;

Tsignal_R::Tsignal_R(){
   
}

Tsignal_R::Tsignal_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_) : fChain(0) 
{
   Init(tree, det_N_);
}

Tsignal_R::~Tsignal_R()
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
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tsignal_R::GetEntry(Long64_t entry)
{
// Read contents of entry.
   current_entry = entry;
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tsignal_R::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Tsignal_R::Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   det_N = det_N_;
   /*
   if(det_N.count(Tomography::MG)>0){
      StripAmpl_MG = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MG_ped = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MG_corr = new Float_t[det_N[Tomography::MG]*MG_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
   }
   if(det_N.count(Tomography::MGv2)>0){
      StripAmpl_MGv2 = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MGv2_ped = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_MGv2_corr = new Float_t[det_N[Tomography::MGv2]*MGv2_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
   }
   if(det_N.count(Tomography::CM)>0){
      StripAmpl_CM = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_CM_ped = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
      StripAmpl_CM_corr = new Float_t[det_N[Tomography::CM]*CM_Detector::Nchannel*Tomography::get_instance()->get_Nsample()];
   }
   */
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   current_entry = -1;
   
   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("evttime",&evttime, &b_evttime);
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      int current_channel = Tomography::Static_Detector[type_it->first]->get_Nchannel();
      StripAmpl[type_it->first] = new Float_t[(type_it->second)*current_channel*(Tomography::get_instance()->get_Nsample())];
      StripAmpl_ped[type_it->first] = new Float_t[(type_it->second)*current_channel*(Tomography::get_instance()->get_Nsample())];
      StripAmpl_corr[type_it->first] = new Float_t[(type_it->second)*current_channel*(Tomography::get_instance()->get_Nsample())];
      string current_name = Tomography::Static_Detector[type_it->first]->Name();
      fChain->SetBranchAddress(("StripAmpl_"+current_name).c_str(),StripAmpl[type_it->first],&b_StripAmpl[type_it->first]);
      fChain->SetBranchAddress(("StripAmpl_"+current_name+"_ped").c_str(),StripAmpl_ped[type_it->first],&b_StripAmpl_ped[type_it->first]);
      fChain->SetBranchAddress(("StripAmpl_"+current_name+"_corr").c_str(),StripAmpl_corr[type_it->first],&b_StripAmpl_corr[type_it->first]);
   }
   /*
   if(det_N.count(Tomography::MG)>0){
	   fChain->SetBranchAddress("StripAmpl_MG", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_MG_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_MG_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
   if(det_N.count(Tomography::MGv2)>0){
      fChain->SetBranchAddress("StripAmpl_MGv2", StripAmpl_MGv2, &b_StripAmpl_MGv2);
      fChain->SetBranchAddress("StripAmpl_MGv2_ped", StripAmpl_MGv2_ped, &b_StripAmpl_MGv2_ped);
      fChain->SetBranchAddress("StripAmpl_MGv2_corr", StripAmpl_MGv2_corr, &b_StripAmpl_MGv2_corr);
   }
	if(det_N.count(Tomography::CM)>0){
	   fChain->SetBranchAddress("StripAmpl_CM", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_CM_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_CM_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
   */
 Notify();
}

Bool_t Tsignal_R::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tsignal_R::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tsignal_R::Cut(Long64_t /*entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
bool Tsignal_R::GetNext(){
   current_entry++;
   if(current_entry<0 || current_entry>=fChain->GetEntries()){
      current_entry = fChain->GetEntries();
      return false;
   }
   else{
      LoadTree(current_entry);
      GetEntry(current_entry);
      return true;
   }
}
long Tsignal_R::GetCurrentEntry() const{
   return current_entry;
}
template<typename T>
vector<vector<T> > Tsignal_R::get_ampl(Tomography::det_type type_, unsigned short det_n_){
   vector<vector<T> > return_array;
   int current_channel = Tomography::Static_Detector[type_]->get_Nchannel();
   return_array = vector<vector<T> >(current_channel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
   for(int i=0;i<current_channel;i++){
      for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
         return_array[i][j] = reinterpret_cast<Float_t(*)[current_channel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_corr[type_])[det_n_][i][j];
      }
   }
   /*
   if(type_ == Tomography::MG){
      return_array = vector<vector<T> >(MG_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MG_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_corr)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::MGv2){
      return_array = vector<vector<T> >(MGv2_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MGv2_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_corr)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::CM){
      return_array = vector<vector<T> >(CM_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
         for(int i=0;i<CM_Detector::Nchannel;i++){
            for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
               return_array[i][j] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_corr)[det_n_][i][j];
            }
         }
   }
   else{

   }
   */
   return return_array;
}
template vector<vector<float> > Tsignal_R::get_ampl(Tomography::det_type type_, unsigned short det_n_);
template vector<vector<double> > Tsignal_R::get_ampl(Tomography::det_type type_, unsigned short det_n_);

template<typename T>
vector<vector<T> > Tsignal_R::get_ampl_ped(Tomography::det_type type_, unsigned short det_n_){
   vector<vector<T> > return_array;
   int current_channel = Tomography::Static_Detector[type_]->get_Nchannel();
   return_array = vector<vector<T> >(current_channel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
   for(int i=0;i<current_channel;i++){
      for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
         return_array[i][j] = reinterpret_cast<Float_t(*)[current_channel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_ped[type_])[det_n_][i][j];
      }
   }
   /*
   if(type_ == Tomography::MG){
      return_array = vector<vector<T> >(MG_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MG_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG_ped)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::MGv2){
      return_array = vector<vector<T> >(MGv2_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MGv2_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2_ped)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::CM){
      return_array = vector<vector<T> >(CM_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
         for(int i=0;i<CM_Detector::Nchannel;i++){
            for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
               return_array[i][j] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM_ped)[det_n_][i][j];
            }
         }
   }
   else{

   }
   */
   return return_array;
}
template vector<vector<float> > Tsignal_R::get_ampl_ped(Tomography::det_type type_, unsigned short det_n_);
template vector<vector<double> > Tsignal_R::get_ampl_ped(Tomography::det_type type_, unsigned short det_n_);

template<typename T>
vector<vector<T> > Tsignal_R::get_ampl_raw(Tomography::det_type type_, unsigned short det_n_){
   vector<vector<T> > return_array;
   int current_channel = Tomography::Static_Detector[type_]->get_Nchannel();
   return_array = vector<vector<T> >(current_channel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
   for(int i=0;i<current_channel;i++){
      for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
         return_array[i][j] = reinterpret_cast<Float_t(*)[current_channel][Tomography::get_instance()->get_Nsample()]>(StripAmpl[type_])[det_n_][i][j];
      }
   }
   /*
   if(type_ == Tomography::MG){
      return_array = vector<vector<T> >(MG_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MG_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MG_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MG)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::MGv2){
      return_array = vector<vector<T> >(MGv2_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
      for(int i=0;i<MGv2_Detector::Nchannel;i++){
         for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
            return_array[i][j] = reinterpret_cast<Float_t(*)[MGv2_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_MGv2)[det_n_][i][j];
         }
      }
   }
   else if(type_ == Tomography::CM){
      return_array = vector<vector<T> >(CM_Detector::Nchannel,vector<T>(Tomography::get_instance()->get_Nsample(),0));
         for(int i=0;i<CM_Detector::Nchannel;i++){
            for(int j=0;j<Tomography::get_instance()->get_Nsample();j++){
               return_array[i][j] = reinterpret_cast<Float_t(*)[CM_Detector::Nchannel][Tomography::get_instance()->get_Nsample()]>(StripAmpl_CM)[det_n_][i][j];
            }
         }
   }
   else{

   }
   */
   return return_array;
}
template vector<vector<float> > Tsignal_R::get_ampl_raw(Tomography::det_type type_, unsigned short det_n_);
template vector<vector<double> > Tsignal_R::get_ampl_raw(Tomography::det_type type_, unsigned short det_n_);
