#define Tanalyse_R_cxx
#include "Tanalyse_R.h"
#include <iostream>

using std::cout;
using std::endl;

Tanalyse_R::Tanalyse_R(){
   
}

Tanalyse_R::Tanalyse_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/* if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run6_Pb_analyse.root");
      if (!f) {
         f = new TFile("run6_Pb_analyse.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }*/
   Init(tree, det_N_);
}

Tanalyse_R::~Tanalyse_R()
{
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      delete NClus[type_it->first];
      delete Spark[type_it->first];
      delete ClusAmpl[type_it->first];
      delete ClusSize[type_it->first];
      delete ClusPos[type_it->first];
      delete ClusTOT[type_it->first];
      delete ClusT[type_it->first];
      delete ClusMaxStrip[type_it->first];
      delete ClusMaxSample[type_it->first];
      delete ClusMaxStripAmpl[type_it->first];
      delete StripMaxAmpl[type_it->first];
   }
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tanalyse_R::GetEntry(Long64_t entry)
{
// Read contents of entry.
   current_entry = entry;
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tanalyse_R::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}
bool Tanalyse_R::GetNext(){
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
long Tanalyse_R::GetCurrentEntry() const{
   return current_entry;
}
void Tanalyse_R::Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   det_N = det_N_;
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   current_entry = -1;
   
   fChain->SetBranchAddress("evn", &evn, &b_evn);
   //evttime = 0;
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      NClus[type_it->first] = new int[type_it->second];
      string current_name = Tomography::Static_Detector[type_it->first]->Name();
      fChain->SetBranchAddress((current_name+"_NClus").c_str(), NClus[type_it->first], &b_NClus[type_it->first]);
      Spark[type_it->first] = new int[type_it->second];
      for(int i = 0;i<(type_it->second);i++){
         Spark[type_it->first][i] = 0;
      }
      //fChain->SetBranchAddress((current_name+"_Spark").c_str(), Spark[type_it->first], &b_Spark[type_it->first]);
      int current_MaxNClus = Tomography::Static_Detector[type_it->first]->get_MaxNClus();
      ClusAmpl[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusAmpl").c_str(), ClusAmpl[type_it->first], &b_ClusAmpl[type_it->first]);
      ClusSize[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusSize").c_str(), ClusSize[type_it->first], &b_ClusSize[type_it->first]);
      ClusPos[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusPos").c_str(), ClusPos[type_it->first], &b_ClusPos[type_it->first]);
      ClusMaxStripAmpl[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusMaxStripAmpl").c_str(), ClusMaxStripAmpl[type_it->first], &b_ClusMaxStripAmpl[type_it->first]);
      ClusMaxStrip[type_it->first] = new Int_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusMaxStrip").c_str(), ClusMaxStrip[type_it->first], &b_ClusMaxStrip[type_it->first]);
      ClusMaxSample[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusMaxSample").c_str(), ClusMaxSample[type_it->first], &b_ClusMaxSample[type_it->first]);
      ClusTOT[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusTOT").c_str(), ClusTOT[type_it->first], &b_ClusTOT[type_it->first]);
      ClusT[type_it->first] = new Double_t[(type_it->second)*current_MaxNClus];
      fChain->SetBranchAddress((current_name+"_ClusT").c_str(), ClusT[type_it->first], &b_ClusT[type_it->first]);
      StripMaxAmpl[type_it->first] = new Double_t[(type_it->second)*(Tomography::Static_Detector[type_it->first]->get_Nchannel())];
      fChain->SetBranchAddress((current_name+"_StripMaxAmpl").c_str(), StripMaxAmpl[type_it->first], &b_StripMaxAmpl[type_it->first]);
   }
   /*
   if(det_N.count(Tomography::MG)>0){
      MG_NClus = new int[det_N[Tomography::MG]];
      fChain->SetBranchAddress("MG_NClus", MG_NClus, &b_MG_NClus);
      MG_Spark = new int[det_N[Tomography::MG]];
      for(int i = 0;i<det_N[Tomography::MG];i++){
         MG_Spark[i] = 0;
      }
      //fChain->SetBranchAddress("MG_Spark", MG_Spark, &b_MG_Spark);
      MG_ClusAmpl = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusAmpl", MG_ClusAmpl, &b_MG_ClusAmpl);
      MG_ClusSize = new Double_t[det_N[Tomography::MG]][300];
      MG_ClusPos = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusPos", MG_ClusPos, &b_MG_ClusPos);
      fChain->SetBranchAddress("MG_ClusSize", MG_ClusSize, &b_MG_ClusSize);
      MG_ClusMaxStripAmpl = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusMaxStripAmpl", MG_ClusMaxStripAmpl, &b_MG_ClusMaxStripAmpl);
      MG_ClusMaxSample = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusMaxSample", MG_ClusMaxSample, &b_MG_ClusMaxSample);
      MG_ClusTOT = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusTOT", MG_ClusTOT, &b_MG_ClusTOT);
      MG_ClusT = new Double_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusT", MG_ClusT, &b_MG_ClusT);
      MG_ClusMaxStrip = new Int_t[det_N[Tomography::MG]][300];
      fChain->SetBranchAddress("MG_ClusMaxStrip", MG_ClusMaxStrip, &b_MG_ClusMaxStrip);
      MG_StripMaxAmpl = new Double_t[det_N[Tomography::MG]][MG_Detector::Nchannel];
      fChain->SetBranchAddress("MG_StripMaxAmpl", MG_StripMaxAmpl, &b_MG_StripMaxAmpl);
   }
   if(det_N.count(Tomography::MGv2)>0){
      MGv2_NClus = new int[det_N[Tomography::MGv2]];
      fChain->SetBranchAddress("MGv2_NClus", MGv2_NClus, &b_MGv2_NClus);
      MGv2_Spark = new int[det_N[Tomography::MGv2]];
      for(int i = 0;i<det_N[Tomography::MGv2];i++){
         MGv2_Spark[i] = 0;
      }
      //fChain->SetBranchAddress("MGv2_Spark", MGv2_Spark, &b_MGv2_Spark);
      MGv2_ClusAmpl = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusAmpl", MGv2_ClusAmpl, &b_MGv2_ClusAmpl);
      MGv2_ClusSize = new Double_t[det_N[Tomography::MGv2]][300];
      MGv2_ClusPos = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusPos", MGv2_ClusPos, &b_MGv2_ClusPos);
      fChain->SetBranchAddress("MGv2_ClusSize", MGv2_ClusSize, &b_MGv2_ClusSize);
      MGv2_ClusMaxStripAmpl = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusMaxStripAmpl", MGv2_ClusMaxStripAmpl, &b_MGv2_ClusMaxStripAmpl);
      MGv2_ClusMaxSample = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusMaxSample", MGv2_ClusMaxSample, &b_MGv2_ClusMaxSample);
      MGv2_ClusTOT = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusTOT", MGv2_ClusTOT, &b_MGv2_ClusTOT);
      MGv2_ClusT = new Double_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusT", MGv2_ClusT, &b_MGv2_ClusT);
      MGv2_ClusMaxStrip = new Int_t[det_N[Tomography::MGv2]][300];
      fChain->SetBranchAddress("MGv2_ClusMaxStrip", MGv2_ClusMaxStrip, &b_MGv2_ClusMaxStrip);
      MGv2_StripMaxAmpl = new Double_t[det_N[Tomography::MGv2]][MGv2_Detector::Nchannel];
      fChain->SetBranchAddress("MGv2_StripMaxAmpl", MGv2_StripMaxAmpl, &b_MGv2_StripMaxAmpl);
   }
   */
   Notify();
}

Bool_t Tanalyse_R::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tanalyse_R::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
